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
    for seq in seqGenerator.next():
        n = len(seq)
        for i in xrange(n - k + 1):
            kmer = seq[i:i+k]
            kmer_comp = getReverseCompliment(kmer)
            h[kmer] = h.get(kmer,0) + 1
            h[kmer_comp] = h.get(kmer_comp,0) + 1
    return h

######################################################

def getKppmerNextList( kmer, kppmerHash ):
    NUCLEO = ['A', 'C', 'G', 'T']

    l = []
    for n in NUCLEO:
        if (kmer + n) in kppmerHash:
            l.append(kmer + n)
    return l

def getKppmerPrevList( kmer, kppmerHash ):
    NUCLEO = ['A', 'C', 'G', 'T']

    l = []
    for n in NUCLEO:
        if (n + kmer) in kppmerHash:
            l.append(n + kmer)
    return l

######################################################

def buildGraph( seqGenerator, k ):
    print "Building k-mer hash",
    kmerHash = generateHash(seqGenerator, k)
    print "OK ({0})".format(len(kmerHash))
    print "Building (k+1)-mer hash",
    kppmerHash = generateHash(seqGenerator, k+1)
    print "OK ({0})".format(len(kppmerHash))

    graph = {}

    # Build graph vertices
    print "Detecting vertices k-mers"
    for seq in seqGenerator.next():
        if seq == "":
            continue
        kmer_i_last = len(seq) - k
        for kmer_i, kmer in kmerGenerator(seq, k):
            # Get k+1-mers
            n = getKppmerNextList(kmer, kppmerHash)
            p = getKppmerPrevList(kmer, kppmerHash)

            # Start or end of read or contig
            if kmer_i == 0 or kmer_i == kmer_i_last:
                if len(n) == 0 or len(p) == 0:
                    graph[kmer] = Vertex()
                elif len(n) > 1 or len(p) > 1:
                    graph[kmer] = Vertex()
            else:
                if kmerHash[kmer] < 2:
                    continue
                elif len(n) < 2 and len(p) < 2:
                    continue

                graph[kmer] = Vertex()
            if kmer in graph:
                graph[kmer].inDegree  = len(p)
                graph[kmer].outDegree = len(n)

    # for i in graph:
    #     print i

    # Build graph edges
    print "Creating edges between vertices"
    kmers_queue = Queue.Queue()
    for vertex in graph:
        kmers_queue.put( (vertex, (vertex, kmerHash[vertex], 1)) )
        while not kmers_queue.empty():
            kmer_prev,(kmer_edge,kmer_edge_cover,kmer_edge_count) = kmers_queue.get()
            kppmersNext = getKppmerNextList(kmer_prev, kppmerHash)
            for kppmer in kppmersNext:
                kmer_next   = kppmer[1:]
                edge_next   = kmer_edge + kppmer[-1]
                edge_cover = kmer_edge_cover + kppmerHash[kppmer]
                if kmer_next in graph:
                    graph[vertex].edges.append(Edge(edge_next, vertex, kmer_next))
                    graph[vertex].edges[-1].avg_cover = edge_cover / kmer_edge_count
                else:
                    kmers_queue.put( (kmer_next, (edge_next, edge_cover, kmer_edge_count + 1)) )

    del kmerHash
    del kppmerHash

    return graph

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
    parser.add_argument('-o', metavar="out", action="store", dest="out", default="out", 
        help="output file prefix (default: out)")

    args = parser.parse_args()

    graph = buildGraph(getGenerator(args.source, args.fastq), args.k)
    drawGraph(graph, args.out + ".dot")
    writeFasta(graph, args.out + ".fasta")

if __name__ == "__main__":
    main()