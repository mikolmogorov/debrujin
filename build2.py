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


def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]


def build_hash(seq, k):
	h = {}
	for _, kmer in iter_kmers(seq, k):
		h[kmer] = h.get(kmer, 0) + 1
	return h


def build_compressed_graph(inputSeq, kmerLen):
	graph = {}
	kmerHash = build_hash(inputSeq, kmerLen)
	kPlusOneHash = build_hash(inputSeq, kmerLen + 1)

	prevVertex = None
	edgeSeq = ""

	for i, kmer in iter_kmers(inputSeq, kmerLen):
		edgeSeq += kmer[-1]
		#check, if kmer is a vertex in compressed graph
		isVertex = False

		if i == 0:
			prevVertex = kmer
			edgeSeq = kmer
			graph[kmer] = Vertex()
		elif i == len(inputSeq) - kmerLen:
			isVertex = True
		elif kmer in graph:
			isVertex = True
		elif kmerHash[kmer] > 1:
			#search for different k+1-mers
			nRight, nLeft = 0, 0
			for n in NUCLEOTIDES:
				if kmer + n in kPlusOneHash:
					nRight += 1

				if n + kmer in kPlusOneHash:
					nLeft += 1
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
				graph[prevVertex].inDegree += 1

			prevVertex = kmer
			edgeSeq = kmer

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
