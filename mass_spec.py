import math

mass_table = {}
result = ""

with open("mass_table.txt",'r') as f:
    for line in f:
        key, value = line.split()
        mass_table[key] = float(value) 

with open("rosalind_spec.txt", 'r') as f:
    L = [float(line.strip()) for line in f if line.strip()]

for i in range(1,len(L)):
    diff = L[i] - L[i-1]
    
    best_aa = None
    best_error = float("inf")

    for aa, mass in mass_table.items():
        error = abs(diff - mass)

        if error < best_error:
            best_error = error
            best_aa = aa

    result += best_aa

with open("rosalind_spec_output.txt",'w') as f:
    f.write(result)