"""
Problem
A graph whose nodes have all been labeled can be represented by an adjacency list, 
in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. 
That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. 
The directed edge with tail v and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v).

For a collection of strings and a positive integer k, 
the overlap graph for the strings is a directed graph Ok in which each string is represented by a node, 
and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, 
as long as s≠t; we demand s≠t to prevent directed loops in the overlap graph (although directed cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
Return: The adjacency list corresponding to O3. You may return edges in any order.
"""
# import libraries
import networkx as nx
from Bio import SeqIO

# Parse the sequences
records = list(SeqIO.parse("data/rosalind_grph.txt", "fasta"))

# Create graph
G = nx.DiGraph()

# Add nodes where DNA sequences are the nodes
for record in records:
    G.add_node(record.id, sequence=str(record.seq))

# Add edge where an edge u → v exists
# If the suffix of length k of u matches the prefix of length k of v
k = 3
for u in G.nodes:
    for v in G.nodes:
        if u != v:
            seq_u = G.nodes[u]["sequence"]
            seq_v = G.nodes[v]["sequence"]
            if seq_u[-k:] == seq_v[:k]:
                G.add_edge(u, v)

with open("output/rosalind_grph_output.txt", 'w') as f:
    for u, v in G.edges():
        f.write(f"{u} {v}\n")
