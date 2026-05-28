"""Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: 
k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). 
Assume that any two organisms can mate."""


#mendel first law ratios
# heterozygous (m*m) = 3 dominant :1 recessive or 75%. 
# homozygous (k*k, n*n) = 100% homozygous
# homozygous dominant and homozygous recessive (k*n)= 2:2 or 50%
# homozygous and heterzygous (k*m, m*n) = 1:1 or 50%

with open("rosalind_iprb.txt") as f:
    k, m, n = map(int, f.readline().split())

p=k+m+n

#here we will subtract the cases where dominance is NOT produced as these cases are fewer and therfore easier to then caluclate the probability of dominance. that why k*k is not calculated.
Hetero = (1/4) * (m/p)*((m-1)/(p-1)) # probability of two heterozygous producing a homozygous recessive.
homoRec = (n/p)*((n-1)/(p-1)) # probability of two homozygous reccesive producing homozygous reccesive
HeteroAndHomoRecOne = (1/2) * (m/p)*(n/(p-1)) 
HeteroAndHomoRecTwo = (1/2) * (n/p)*(m/(p-1))

probabilityDominant = 1 - (Hetero + homoRec + HeteroAndHomoRecOne + HeteroAndHomoRecTwo)
print(probabilityDominant)