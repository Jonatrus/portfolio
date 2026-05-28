"""
Problem
A prefix of a length n string s is a substring s[1:j]; a suffix of s is a substring s[k:n].

The failure array of s is an array P of length n for which P[k] is the length of the longest substring s[j:k] that is equal to some prefix s[1:k−j+1], where j cannot equal 1 (otherwise, P[k] would always equal k). By convention, P[1] = 0.

Given:
A DNA string s (of length at most 100 kbp) in FASTA format.

Return:
The failure array of s.
"""
# Load library
from Bio import SeqIO

# Parse data (expects exactly one sequence)
record = SeqIO.read("data/rosalind_kmp.txt", "fasta")
s = str(record.seq)
m = len(s)

# Initialise the failure array (also called prefix table) with zeros
# failureArr[k] will store the length of the longest prefix that is also a suffix
# for the substring s[0:k+1]
failureArr = [0]*m

# j tracks the length of the current matching prefix
# (i.e. how many characters we've matched so far)
j = 0

# i is the current position in the string we are computing the value for
# we start at 1 because failureArr[0] is always 0 by definition
i = 1

# Build failure array
while i < m:

    # Case 1: characters match
    # This means we can extend the current prefix-suffix match
    if s[i] == s[j]:
        j += 1                     # increase length of current match
        failureArr[i] = j          # store this length at position i
        i += 1                     # move to next character

    else:
        # Case 2: characters do not match

        if j != 0:
            # Instead of restarting from scratch,
            # fall back to the previous longest prefix length
            # This is the key optimisation of KMP
            j = failureArr[j - 1]

        else:
            # If no prefix matched so far (j == 0),
            # then failureArr[i] stays 0
            failureArr[i] = 0
            i += 1                 # move to next character


with open("output/rosalind_kmp_output.txt",'w') as f:
    f.write(' '.join(map(str, failureArr)) + '\n')
