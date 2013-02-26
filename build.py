#!/usr/bin/env python

import sys
from collections import defaultdict
from Queue import Queue


class Vertex(object):
    def __init__(self):
        self.edges = []
        self.in_degree = 0

    def __str__(self):
        return ",".join(map(str, self.edges))

    __repr__ = __str__

    def __eq__(self, other):
        return self.edges == other.edges and self.in_degree == other.in_degree


class Edge(object):
    def __init__(self, seq):
        self.seq = seq
        self.vertex = None
        self.degree = 1

    def __str__(self):
        first, last = self.seq[:3], self.seq[-3:]
        return "{0}..{1} ({2}) {3}".format(first, last,
                                           len(self), self.degree)

    def __len__(self):
        return len(self.seq)

    def __eq__(self, other):
        return all([self.seq == other.seq,
                    self.vertex == other.vertex,
                    self.degree == other.degree])



def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]


def build_uncompressed_graph(seq, k):
    n = len(seq)
    g, seen = {}, defaultdict(int)
    for i, kmer in iter_kmers(seq, k):
        seen[kmer] += 1
        if seen[kmer] > 1 or i in [0, n - k]:
            g[kmer] = Vertex()

    return g


def build_compressed_graph(inputSeq, kmerLen):
    graph = build_uncompressed_graph(inputSeq, kmerLen)

    #build de-brujin graph
    #contained as set of vertexes
    curEdge = None
    for i, kmer in iter_kmers(inputSeq, kmerLen):
        if curEdge is not None:
            graph[curKmer].edges[curEdge].seq += kmer[-1]

        if kmer in graph:
            if curEdge != None:
                #seerch for the same edges
                edge = graph[curKmer].edges[curEdge]
                deleted = False
                for e in graph[curKmer].edges:
                    if e == edge:
                        #print "bbb"
                        continue
                    if e.seq == edge.seq:
                        e.degree += 1
                        del graph[curKmer].edges[curEdge]
                        #print "aaa"
                        deleted = True
                        break
                if not deleted:
                    graph[curKmer].edges[curEdge].vertex = kmer
                    graph[kmer].in_degree += 1
            #print kmer
            curKmer = kmer
            graph[curKmer].edges.append(Edge(kmer))
            curEdge = len(graph[curKmer].edges) - 1

    #TODO: ugly hack
    del graph[curKmer].edges[curEdge]

    #simplify graph
    print "Simplifying graph"
    wfsqueue = Queue()
    visited = set()
    firstKmer = inputSeq[0 : kmerLen]
    visited.add(firstKmer)
    wfsqueue.put(firstKmer)

    counter = 1
    while not wfsqueue.empty():
        kmer = wfsqueue.get()
        #print graph[kmer].edges
        for edge in graph[kmer].edges:
            curVertex = edge.vertex

            if curVertex in visited:
                continue

            while len(graph[curVertex].edges) == 1 and graph[curVertex].in_degree == 1:
                edge.seq += graph[curVertex].edges[0].seq[kmerLen:]
                edge.vertex = graph[curVertex].edges[0].vertex
                del graph[curVertex]
                curVertex = edge.vertex

            print counter, "/", len(graph)
            counter += 1

            #TODO: assert fails
            #assert(curVertex not in visited)
            visited.add(curVertex)
            wfsqueue.put(curVertex)

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
            dot_file.write("""{0} -> {1.vertex} [label = "{1}"];\n"""
                           .format(v, edge))
    dot_file.write("}")


if __name__ == "__main__":
    try:
        seq_path, dot_path, k = sys.argv[1:]
        k = int(k)
    except ValueError:
        print("Usage: build.py SEQ_PATH DOT_PATH K")
    else:
        graph = build_compressed_graph(read_fasta(seq_path), k)
        write_dot(graph, open(dot_path, "w"))
