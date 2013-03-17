#!/usr/bin/env python2

import sys
from Queue import Queue



class Vertex:
    def __init__(self):
        self.in_edges = []
        self.out_edges = []


class Edge:
    def __init__(self, seq):
        self.seq = seq
        self.vertex = None
        self.coverage = 0.0

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return "L={0} C={1:5.2f}".format(len(self), self.coverage)

class FastqSource:
    def __init__(self, filename):
        self.fastq = filename
    def reads(self):
        fin = open(self.fastq, "r")
        counter = 0
        for line in fin:
            counter = (counter + 1) % 4
            if counter == 2:
                yield line.strip()
                yield getReverseComplement(line.strip())
        fin.close()


def getReverseComplement( seq ):
    COMPLEMENT = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    s = []
    for i in seq:
        s.append(COMPLEMENT[i])
    return "".join(s[::-1])


def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def getKmerNextList(kmer, kmerHash):
    l = []
    base = kmer[1:]
    for n in NUCLEOTIDES:
        if (base + n) in kmerHash:
            l.append(n)
    return l,base

def getKmerPrevList(kmer, kmerHash):
    l = []
    base = kmer[:-1]
    for n in NUCLEOTIDES:
        if (n + base) in kmerHash:
            l.append(n)
    return l,base

def build_hash(fastqSrc, k):
    print "Buliding kmer hash"
    h = {}
    for read in fastqSrc.reads():
        for i, kmer in iter_kmers(read, k):
            h[kmer] = h.get(kmer, 0) + 1
    return h


def count_vertex(kmerHash):
    count = 0
    for kmer in kmerHash:
        r, rbase = getKmerNextList(kmer, kmerHash)
        l, lbase = getKmerPrevList(kmer, kmerHash)
        if len(r) != 1 or len(l) != 1:
            count += 1
    return count


def goBranch(kmer, kmerHash, visited):
    coverage = 0.0
    nVertex = 0
    seq = kmer[:-1]
    while True:
        coverage += kmerHash[kmer]
        nVertex += 1
        seq += kmer[-1]

        r, rbase = getKmerNextList(kmer, kmerHash)
        l, lbase = getKmerPrevList(kmer, kmerHash)

        if len(r) == 1 and len(l) == 1:
            visited.add(kmer)
            kmer = rbase + r[0]
        else:
            return kmer, coverage / nVertex, seq, r


def build(fastqSrc, kmerLen):
    graph = {}
    tot_cover = 0
    num_edges = 0
    kmerHash = build_hash(fastqSrc, kmerLen)
    print "Building graph"
    #print count_vertex(kmerHash)
    visited = set()
    for newKmer in kmerHash:
        if newKmer in visited:
            continue

        bfsqueue = Queue()
        bfsqueue.put((newKmer, None))
        while not bfsqueue.empty():
            kmer, prevVertex = bfsqueue.get()

            vertex, cover, seq, outBranch = goBranch(kmer, kmerHash, visited)

            if not vertex in graph:
                graph[vertex] = Vertex()
            if prevVertex:
                graph[prevVertex].out_edges.append(Edge(seq))
                graph[prevVertex].out_edges[-1].vertex = vertex
                graph[prevVertex].out_edges[-1].coverage = cover
                graph[vertex].in_edges.append(prevVertex)
                tot_cover += cover
                num_edges += 1

            if vertex in visited:
                continue
            visited.add(vertex)

            for nr in outBranch:
                newMer = vertex[1:] + nr
                bfsqueue.put((newMer, vertex))

    #print len(visited), len(kmerHash)
    return graph, tot_cover / num_edges


def cut_graph(graph, threshold):
    print "Removing tips",
    visited = set()
    toDelete = []
    marked = []
    for v in graph:
        #cut vertex with in_degree = 1, out_degree = 0
        if (len(graph[v].in_edges) == 1 and
                len(graph[v].out_edges) == 0):
            prev = graph[v].in_edges[0]
            inEdge = [e for e in graph[prev].out_edges if e.vertex == v][0]
            if inEdge.coverage < threshold:
                toDelete.append(v)
                marked.append(prev)
                graph[prev].out_edges = [e for e in graph[prev].out_edges if e.vertex != v]

        if (len(graph[v].in_edges) == 0 and
                len(graph[v].out_edges) == 1):
            nextV = graph[v].out_edges[0].vertex
            if graph[v].out_edges[0].coverage < threshold:
                toDelete.append(v)
                marked.append(nextV)
                graph[nextV].in_edges.remove(v)

    print "OK", len(toDelete)
    for elem in toDelete:
        if elem in graph:
            del graph[elem]

    for vertex in marked:
        #print len(graph[vertex].in_edges), len(graph[vertex].out_edges)
        if not vertex in graph:
            continue
        if len(graph[vertex].in_edges) != 1 or len(graph[vertex].out_edges) != 1:
            continue

        prev_v = graph[vertex].in_edges[0]
        next_v = graph[vertex].out_edges[0].vertex
        for edge in graph[prev_v].out_edges:
            if edge.vertex == vertex:
                edge.vertex = next_v
                edge.seq += graph[vertex].out_edges[0].seq
                edge.coverage = (edge.coverage + graph[vertex].out_edges[0].coverage) / 2
                break
        graph[next_v].in_edges.remove(vertex)
        graph[next_v].in_edges.append(prev_v)
        del graph[vertex]
    return graph


def read_fasta(filename):
    fasta = open(filename, "r")
    seq = ""
    for line in fasta:
        if line[0] == '>':
            continue
        line = line.strip('\n')
        seq += line
    return seq


def write_dot(graph, dot_file):
    print "Writing dot file"
    dot_file.write("digraph {\n")
    v_ids = {}
    id_count = 0
    for v in graph:
        if v not in v_ids:
            v_ids[v] = id_count
            id_count += 1
        for edge in graph[v].out_edges:
            if edge.vertex not in v_ids:
                v_ids[edge.vertex] = id_count
                id_count += 1
            dot_file.write("""{0} -> {1} [label = "{2}"];\n""".format(v_ids[v], v_ids[edge.vertex], edge))
    dot_file.write("}")

def write_fasta(graph, fileName):
    print "Exporting edges to FASTA file"
    fd = open(fileName, "w")

    counter = 1
    for kmer in graph:
        for edge in graph[kmer].out_edges:
            fd.write("> edge" + str(counter) + "\n")
            seq = edge.seq
            while seq != "":
                fd.write(seq[:80] + "\n")
                seq = seq[80:]
            counter += 1
            fd.write("\n")
    fd.close()

if __name__ == "__main__":
    try:
        seq_path, k, threshold = sys.argv[1:]
        k = int(k)
    except ValueError:
        print("Usage: build.py SEQ_PATH K THRESHOLD")
    else:
        graph, avg_cover = build(FastqSource(seq_path), k)
        print "Average coverage:", avg_cover
        cut_graph(graph, avg_cover * float(threshold) / 100)
        prefix = seq_path.split(".")[0]
        write_dot(graph, open(prefix + ".dot", "w"))
        write_fasta(graph, prefix + ".fasta")
