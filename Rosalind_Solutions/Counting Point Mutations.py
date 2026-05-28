"""
PROBLEM
Given two strings s
 and t
 of equal length, the Hamming distance between s
 and t
, denoted dH(s,t)
, is the number of corresponding symbols that differ in s
 and t
. See Figure 2.

Given: Two DNA strings s
 and t
 of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t)
.
"""
with open("data/rosalind_hamm.txt",'r') as f:
    s,t = f.readlines()

hamming_distance = 0

# Using zip as it will enable parralel comparison of the strings. Avoids using nested for loops.
for a,b in zip(s,t):
    if a != b:
        hamming_distance += 1

print(hamming_distance)
