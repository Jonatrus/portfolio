"""
In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.

The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

Given: A protein string P of length at most 1000 aa.

Return: The total weight of P. Consult the monoisotopic mass table.
"""

with open("rosalind_prtm.txt") as f:
    oneLetterProtein = f.read().strip()

print(oneLetterProtein)
mass_table = {}

with open("mass_table.txt") as f:
    for line in f:
        aa,mass = line.split()
        mass_table[aa] = float(mass)


mass = 0

for aa in oneLetterProtein:
    mass += mass_table[aa]


print(mass)
