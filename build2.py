#!/usr/bin/env python2

import sys
from collections import defaultdict
from Queue import Queue

NUCLEOTIDES = ['A', 'C', 'G', 'T']


class Vertex:
    def __init__(self):
        self.edges = []
        self.inDegree = 0


class Edge:
    def __init__(self, seq):
        self.seq = seq
        self.vertex = None
        self.degree = 1

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        first, last = self.seq[:3], self.seq[-3:]
        return "{0}..{1} ({2}) {3}".format(first, last,
                                           len(self), self.degree)

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
        fin.close()


def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]


def build_hash(fastqSrc, k):
    h = {}
    for read in fastqSrc.reads():
        for _, kmer in iter_kmers(read, k):
            h[kmer] = h.get(kmer, 0) + 1
    return h


def build_compressed_graph(fastqSrc, kmerLen):
    graph = {}
    #toCleanup = set("")
    kmerHash = build_hash(fastqSrc, kmerLen)
    kPlusOneHash = build_hash(fastqSrc, kmerLen + 1)

    prevVertex = None
    edgeSeq = ""


    for read in fastqSrc.reads():
        for i, kmer in iter_kmers(read, kmerLen):
            edgeSeq += kmer[-1]
            #check, if kmer is a vertex in compressed graph
            isVertex = False

            if kmer in graph:
                isVertex = True
            elif i == 0:
                prevVertex = kmer
                edgeSeq = kmer
                graph[kmer] = Vertex()
                #toCleanup.add(kmer)
            elif i == len(read) - kmerLen:
                isVertex = True
                #toCleanup.add(kmer)
            elif kmerHash[kmer] > 1:
                #search for different k+1-mers
                nRight, nLeft = 0, 0
                for n in NUCLEOTIDES:
                    if kmer + n in kPlusOneHash:
                        nRight += 1

                    if n + kmer in kPlusOneHash:
                        nLeft += 1
                assert(nRight > 0 and nLeft > 0)
                if nRight > 1 or nLeft > 1:
                    isVertex = True

            if isVertex:
                if kmer not in graph:
                    graph[kmer] = Vertex()
                #there is edge between previous vertex and this one
                found = False
                for e in graph[prevVertex].edges:
                    if e.seq == edgeSeq:
                        e.degree += 1
                        found = True
                        break
                if not found:
                    graph[prevVertex].edges.append(Edge(edgeSeq))
                    graph[prevVertex].edges[-1].vertex = kmer
                    graph[kmer].inDegree += 1

                prevVertex = kmer
                edgeSeq = kmer

    #simplify_graph(graph, kmerLen)
    for vertex in graph.items():
        #TODO: do something!
        if vertex in graph and vertex not in visited:
            dfs(graph, kmerLen, vertex, visited)
    return graph


def dfs(graph, kmerLen, vertex, visited):
    visited.add(vertex)

    for edge in graph[vertex].edges:
        assert(graph[vertex].inDegree > 0)
        itVertex = edge.vertex
        if itVertex in visited:
            continue

        while len(graph[itVertex].edges) == 1 and graph[itVertex].inDegree == 1:
            #print curVertex
            edge.seq += graph[itVertex].edges[0].seq[kmerLen:]
            edge.vertex = graph[itVertex].edges[0].vertex
            del graph[itVertex]
            itVertex = edge.vertex

        dfs(graph, kmerLen, itVertex, visited)


def simplify_graph(graph, kmerLen):
    print "Simplifying graph"
    bfsqueue = Queue()
    visited = set()
    firstKmer = graph.keys()[0]
    #print firstKmer
    visited.add(firstKmer)
    bfsqueue.put(firstKmer)

    counter = 1
    deleted = set()
    while not bfsqueue.empty():
        kmer = bfsqueue.get()
        #print graph[kmer].edges
        for edge in graph[kmer].edges:
            #print graph[edge.vertex].inDegree
            #print edge.vertex
            #print graph[edge.vertex].end
            assert(graph[edge.vertex].inDegree > 0)
            curVertex = edge.vertex

            if curVertex in visited:
                continue

            while len(graph[curVertex].edges) == 1 and graph[curVertex].inDegree == 1:
                #print curVertex
                edge.seq += graph[curVertex].edges[0].seq[kmerLen:]
                edge.vertex = graph[curVertex].edges[0].vertex
                del graph[curVertex]

                deleted.add(curVertex)

                curVertex = edge.vertex
                assert(curVertex not in deleted)

                #assert(curVertex in graph)
                #print curVertex

            print counter, "/", len(graph)
            counter += 1

            #TODO: assert fails
            #assert(curVertex not in visited)
            visited.add(curVertex)
            bfsqueue.put(curVertex)


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
        graph = build_compressed_graph(FastqSource(seq_path), k)
        write_dot(graph, open(dot_path, "w"))
