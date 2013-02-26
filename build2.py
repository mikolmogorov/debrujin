#!/usr/bin/env python2

import sys
from collections import defaultdict, namedtuple
from Queue import Queue

from build import Vertex, Edge, iter_kmers, read_fasta, write_dot

NUCLEOTIDES = ['A', 'C', 'G', 'T']


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
				graph[prevVertex].in_degree += 1

			prevVertex = kmer
			edgeSeq = kmer

	return graph


if __name__ == "__main__":
    try:
        seq_path, dot_path, k = sys.argv[1:]
        k = int(k)
    except ValueError:
        print("Usage: build.py SEQ_PATH DOT_PATH K")
    else:
        graph = build_compressed_graph(read_fasta(seq_path), k)
        write_dot(graph, open(dot_path, "w"))
