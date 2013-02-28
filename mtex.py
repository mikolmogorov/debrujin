import argparse
import Queue

######################################################

class Vertex(object):
    def __init__(self):
        self.edges = []
        self.inDegree  = 0
        self.outDegree = 0

class Edge(object):
    def __init__(self, seq, kmer1, kmer2):
        self.seq = seq
        self.kmer1 = kmer1
        self.kmer2 = kmer2
        self.avg_cover = 1

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        first, last = self.seq[:3], self.seq[-3:]
        return "{0}..{1} ({2}) {3}".format(first, last,
                                           len(self), self.avg_cover)

######################################################

class FastaSeqGenerator(object):
    def __init__(self, fileName):
        self.fileName = fileName

    def next(self):
        fd = open(self.fileName, "r")
        seq = ""
        for line in fd:
            l = line.strip("\n")
            if l.startswith('>'):
                yield seq
                seq = ""
            else:
                seq += l
        if seq != "":
            yield seq
        fd.close()

class FastqSeqGenerator(object):
    def __init__(self, fileName):
        self.fileName = fileName

    def next(self):
        fd = open(self.fileName, "r")
        i = 0
        for line in fd:
            l = line.strip("\n")
            if l != "":
                if i % 4 == 1:
                    yield l
                i += 1
        fd.close()

def getGenerator( fileName, fastq ):
    if fastq:
        return FastqSeqGenerator(fileName)
    else:
        return FastaSeqGenerator(fileName)

######################################################

def getReverseCompliment( seq ):
    COMPLIMENT = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    s = []
    for i in seq:
        s.append(COMPLIMENT[i])
    return "".join(s[::-1])

######################################################

def kmerGenerator( seq, k ):
    n = len(seq)
    for i in xrange(n - k + 1):
        yield i,seq[i:i+k]
    for i in xrange(n - k + 1):
        yield i,getReverseCompliment(seq[i:i+k])

def generateHash( seqGenerator, k ):
    h = {}
    v = {}
    for seq in seqGenerator.next():
        n = len(seq)
        for i in xrange(n - k + 1):
            kmer = seq[i:i+k]
            kmer_comp = getReverseCompliment(kmer)
            h[kmer] = h.get(kmer,0) + 1
            h[kmer_comp] = h.get(kmer_comp,0) + 1

            v[kmer] = False
            v[kmer_comp] = False
    return h, v

######################################################

def getKppmerNextList( kppmer, kppmerHash ):
    NUCLEO = ['A', 'C', 'G', 'T']

    l = []
    base = kppmer[1:]
    for n in NUCLEO:
        if (base + n) in kppmerHash:
            l.append(n)
    return l,base

def getKppmerPrevList( kppmer, kppmerHash ):
    NUCLEO = ['A', 'C', 'G', 'T']

    l = []
    base = kppmer[:-1]
    for n in NUCLEO:
        if (n + base) in kppmerHash:
            l.append(n)
    return l,base

######################################################

def go( kppmer, kppmerHash, visitedHash, isLeft = False ):
    seq = "" 
    cover = 0
    l,base = None,None
    if isLeft:
        l,base = getKppmerPrevList(kppmer, kppmerHash)
    else:
        l,base = getKppmerNextList(kppmer, kppmerHash)
    while len(l) == 1:
        n = l[0]
        kppmer_next = None
        if isLeft:
            seq = n + seq
            kppmer_next = n + base
            visitedHash[kppmer_next] = True
            l,base = getKppmerPrevList(kppmer_next, kppmerHash)
        else:
            seq = seq + n
            kppmer_next = base + n
            visitedHash[kppmer_next] = True
            l,base = getKppmerNextList(kppmer_next, kppmerHash)
        cover += kppmerHash[kppmer_next]
    return seq,base,cover

def buildGraph( seqGenerator, k ):
    print "Building (k+1)-mer hash",
    kppmerHash,visitedHash = generateHash(seqGenerator, k+1)
    print "OK ({0})".format(len(kppmerHash))
    avg_cover  = 0
    edge_count = 0
    graph = {}
    for kppmer in kppmerHash:
        if visitedHash[kppmer]:
            continue
        # to check if kppmer contains left of right vertex
        p,pkmer = getKppmerPrevList(kppmer, kppmerHash) 
        n,nkmer = getKppmerNextList(kppmer, kppmerHash)
        # not left or right
        if len(p) == 1 and len(n) == 1:
            seqLeft,kmerLeft,coverLeft    = go(kppmer, kppmerHash, visitedHash, True)
            seqRight,kmerRight,coverRight = go(kppmer, kppmerHash, visitedHash, False)
            if kmerLeft not in graph:
                graph[kmerLeft]  = Vertex()
            if kmerRight not in graph:
                graph[kmerRight] = Vertex()

            # Setup out edge to left k-mer
            graph[kmerLeft].edges.append(Edge(seqLeft + kppmer + seqRight, kmerLeft, kmerRight))
            # Setup edge average coverage
            graph[kmerLeft].edges[-1].avg_cover = (coverLeft + kppmerHash[kppmer] + coverRight) 
            graph[kmerLeft].edges[-1].avg_cover /= (len(seqLeft) + len(seqRight) + 1)
            # Setup vertices degrees
            graph[kmerLeft].outDegree += 1
            graph[kmerRight].inDegree += 1
            # Append edge average coverage to count total average coverage
            avg_cover += graph[kmerLeft].edges[-1].avg_cover
            edge_count += 1
        else: 
            if len(p) != 1 and kppmer[:-1] not in graph:
                graph[kppmer[:-1]] = Vertex()
            if len(n) != 1 and kppmer[1:] not in graph:
                graph[kppmer[1:]]  = Vertex()
        visitedHash[kppmer] = True
    del kppmerHash

    return graph, avg_cover/edge_count

######################################################

def filterGraph( graph, avg_cover, t_cover_part ):
    print "Filtering graph"
    threshold = avg_cover * t_cover_part
    graph_new = graph.copy()
    # Remove tails
    for vertex in graph:
        for edge in graph[vertex].edges:
            if edge.avg_cover < threshold and not graph[edge.kmer2].edges and graph[edge.kmer2].inDegree == 1:
                graph_new[vertex].edges.remove(edge)
                del graph_new[edge.kmer2]
    return graph_new

######################################################

def drawGraph( graph, fileName ):
    print "Saving draph to .dot file"
    fd = open(fileName, "w")
    fd.write("digraph {\n")
    for kmer in graph:
        fd.write("\t" + kmer + ";\n")

    for kmer in graph:
        for edge in graph[kmer].edges:
            fd.write("\t" + edge.kmer1 + " -> " + edge.kmer2 + " [label=\"{0}\"];\n".format(edge) )
    fd.write("}\n")
    fd.close()

######################################################

def writeFasta( graph, fileName ):
    print "Exporting edges to FASTA file"
    fd = open(fileName, "w")

    counter = 1
    for kmer in graph:
        for edge in graph[kmer].edges:
            fd.write("> edge" + str(counter) + "\n")
            seq = edge.seq
            while seq != "":
                fd.write(seq[:20] + "\n")
                seq = seq[20:]
            counter += 1
            fd.write("\n")
    fd.close()

######################################################

def main():
    parser = argparse.ArgumentParser(description="DeBrujn Graph Creator")
    parser.add_argument('source', action="store", help="input source (default: FASTA)")
    parser.add_argument('-q', action="store_const", metavar="fastq", dest="fastq", default=False,
        const=True, help="use FASTQ instead of FASTA")
    parser.add_argument('-k', action="store", metavar="k", dest="k", default=55, type=int,
        help="length of k-mers (default: 55)")
    parser.add_argument('-f', action="store", metavar="part", dest="f", default=None, type=int,
        help="filter graph edges with provided part of average cover")
    parser.add_argument('-o', metavar="out", action="store", dest="out", default="out", 
        help="output file prefix (default: out)")

    args = parser.parse_args()

    graph,avg_cover = buildGraph(getGenerator(args.source, args.fastq), args.k)
    if args.f != None:
        graph = filterGraph(graph,avg_cover,0.01 * args.f)
    drawGraph(graph, args.out + ".dot")
    writeFasta(graph, args.out + ".fasta")

if __name__ == "__main__":
    main()