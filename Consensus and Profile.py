"""
PROBLEM
A matrix is a rectangular table of values divided into rows and columns. An m×n
 matrix has m
 rows and n
 columns. Given a matrix A
, we write Ai,j
 to indicate the value found at the intersection of row i
 and column j
.

Say that we have a collection of DNA strings, all having the same length n
. Their profile matrix is a 4×n
 matrix P
 in which P1,j
 represents the number of times that 'A' occurs in the j
th position of one of the strings, P2,j
 represents the number of times that C occurs in the j
th position, and so on (see below).

A consensus string c
 is a string of length n
 formed from our collection by taking the most common symbol at each position; the j
th symbol of c
 therefore corresponds to the symbol having the maximum value in the j
-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
"""
from Bio import SeqIO
import pandas as pd

# Parse the sequences
records = list(SeqIO.parse("data/rosalind_cons.txt", "fasta"))
sequences = [str(record.seq) for record in records]

# Convert sequences into a table
df = pd.DataFrame([list(seq) for seq in sequences])

# Get frequencies of bases
counts = df.apply(pd.Series.value_counts).fillna(0).astype(int)

# Build consensus string
consensus = "".join(df.apply(lambda col: col.value_counts().idxmax()))

with open("output/consensus.txt", "w") as f:
    f.write(consensus + "\n")
    for base in ["A", "C", "G", "T"]:
        row = counts.loc[base].tolist()  # e.g. [5, 1, 0, 0, 5, 5, 0, 0]
        
        numbers = ""
        for count in row:
            numbers = numbers + str(count) + " "
        
        numbers = numbers.strip()  # remove trailing space
        line = base + ": " + numbers
        f.write(line + "\n")