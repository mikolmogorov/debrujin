#!/usr/bin/env python2

import sys
from Queue import Queue

NUCLEOTIDES = ['A', 'C', 'G', 'T']


class Vertex:
    def __init__(self):
        self.edges = []
        self.in_degree = 0
        #self.l = 0
        #self.r = 0


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
                yield getReverseCompliment(line.strip())
        fin.close()


def getReverseCompliment( seq ):
    COMPLIMENT = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    s = []
    for i in seq:
        s.append(COMPLIMENT[i])
    return "".join(s[::-1])


def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]

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

        #assert(len(l) != 0)
        if len(r) == 1 and len(l) == 1:
            visited.add(kmer)
            kmer = rbase + r[0]
        else:
            return kmer, coverage / nVertex, seq, r


def build(fastqSrc, kmerLen):
    graph = {}
    kmerHash = build_hash(fastqSrc, kmerLen)
    print count_vertex(kmerHash)
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
                #graph[vertex].l = len(inBranch)
                #graph[vertex].r = len(outBranch)
            if prevVertex:
                graph[prevVertex].edges.append(Edge(seq))
                graph[prevVertex].edges[-1].vertex = vertex
                graph[prevVertex].edges[-1].coverage = cover
                graph[vertex].in_degree += 1

            if vertex in visited:
                continue
            visited.add(vertex)

            for nr in outBranch:
                newMer = vertex[1:] + nr
                bfsqueue.put((newMer, vertex))

    print len(visited), len(kmerHash)
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
    dot_file.write("digraph {\n")
    for v in graph:
        for edge in graph[v].edges:
            dot_file.write("""{0} -> {1.vertex} [label = "{1}"];\n""".format(v, edge))
    dot_file.write("}")


if __name__ == "__main__":
    try:
        seq_path, dot_path, k = sys.argv[1:]
        k = int(k)
    except ValueError:
        print("Usage: build.py SEQ_PATH DOT_PATH K")
    else:
        graph = build(FastqSource(seq_path), k)
        write_dot(graph, open(dot_path, "w"))
