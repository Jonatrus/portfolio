"""
An array is a structure containing an ordered collection of objects (numbers, strings, other arrays, etc.). 
We let A[k] denote the k-th value in array A. You may like to think of an array as simply a matrix having only one row.

A random string is constructed so that the probability of choosing each subsequent symbol is based on a fixed underlying symbol frequency.

GC-content offers us natural symbol frequencies for constructing random DNA strings. 
If the GC-content is x, then we set the symbol frequencies of C and G equal to x/2 and the symbol frequencies of A and T equal to (1−x)/2. 
For example, if the GC-content is 40%, then as we construct the string, the next symbol is 'G'/'C' with probability 0.2, and the next symbol is 'A'/'T' with probability 0.3.

In practice, many probabilities wind up being very small. In order to work with small probabilities, 
we may plug them into a function that "blows them up" for the sake of comparison. Specifically, 
the common logarithm of x (defined for x>0 and denoted log10(x)) is the exponent to which we must raise 10 to obtain x.

See Figure 1 for a graph of the common logarithm function y=log10(x). 
In this graph, we can see that the logarithm of x-values between 0 and 1 always winds up mapping to y-values between −∞ and 0: x-values near 0 have logarithms close to −∞, 
and x-values close to 1 have logarithms close to 0. Thus, we will select the common logarithm as our function to "blow up" small probability values for comparison.

Given: A DNA string s of length at most 100 bp and an array A containing at most 20 numbers between 0 and 1.

Return: An array B having the same length as A in which B[k] represents the common logarithm of the probability 
        that a random string constructed with the GC-content found in A[k] will match s exactly.
"""
# Import libraries
import math

# If GC-content = x, then:

# P(G) = x / 2
# P(C) = x / 2
# P(A) = (1 − x) / 2
# P(T) = (1 − x) / 2

with open("data/rosalind_prob.txt") as f:
    s, A = f.readlines()
    s = s.strip()
    A = list(map(float, A.split()))
    


B = []
base_probabilities = {}

# Compute the probability of ATCG
for i in range(len(A)):
    x = A[i] # The GC content as represented by x in the problem description.
    base_probabilities["G"] = x / 2
    base_probabilities["C"] = x / 2
    base_probabilities["A"] = (1 - x) / 2
    base_probabilities["T"] = (1 - x) / 2

    # Compute the probability of generating the string s
    perCharacterProbabilities = 0

    for base in s:
        perCharacterProbabilities += math.log10(base_probabilities[base])
    
    B.append(perCharacterProbabilities)
        
with open("output/rosalind_prob_out.txt",'w') as f:
    B = ' '.join(map(str, B))
    f.write(B)

