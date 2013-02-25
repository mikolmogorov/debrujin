#!/usr/bin/env python2

import sys
import Queue

class Vertex:
	def __init__(self):
		self.edges = []
		self.inDegree = 0

class Edge:
	def __init__(self, seq):
		self.seq = seq
		self.vertex = None
		self.degree = 1

def buildGraph(inputSeq, kmerLen):
	graph = {}

	#count kmers
	kmerHash = {}

	for i in xrange(0, len(inputSeq) - kmerLen + 1):
		kmer = inputSeq[i : i + kmerLen]
		if not kmer in kmerHash:
			kmerHash[kmer] = 1
		else:
			kmerHash[kmer] += 1
		
		if kmerHash[kmer] > 1 or i == 0 or i == len(inputSeq) - kmerLen:
			graph[kmer] = Vertex()

	del kmerHash

	#build de-brujin graph
	#contained as set of vertexes
	curEdge = None
	
	for i in xrange(0, len(inputSeq) - kmerLen + 1):
		kmer = inputSeq[i : i + kmerLen]
		if curEdge != None:
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
					graph[kmer].inDegree += 1
			#print kmer
			curKmer = kmer
			graph[curKmer].edges.append(Edge(kmer))
			curEdge = len(graph[curKmer].edges) - 1
	#TODO: ugly hack
	del graph[curKmer].edges[curEdge]

	#simplify graph
	print "Simplifying graph"
	wfsqueue = Queue.Queue()
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

			while len(graph[curVertex].edges) == 1 and graph[curVertex].inDegree == 1:
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

def readFasta(filename):
	fasta = open(filename, "r")
	seq = ""
	for line in fasta:
		if line[0] == '>':
			continue
		line = line.strip('\n')
		seq += line
	return seq

def outputDot(graph, outFile):
	outFile.write("digraph {\n")
	for ver in graph:
		for edge in graph[ver].edges:
			outFile.write(ver + " -> " + str(edge.vertex))
			outFile.write(" [label=\"" + edge.seq[0:3] + "..." + edge.seq[-4:-1])
			outFile.write("(" + str(len(edge.seq)) + ") " + str(edge.degree) + "\"];\n")
	outFile.write("}")

#seq = "AACGTGCGCTAGCTGGCTAGCTAGCCGATAGCTCTAGGCTAGGCGATCGCTACAT"
seq = readFasta(sys.argv[1])
k = 20

graph = buildGraph(seq, k)
outputDot(graph, open(sys.argv[2], "w"))
