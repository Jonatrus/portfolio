import itertools
import math

n = 5

with open("Rosalind coding problems/output/rosalind_perm_output.txt",'w') as f:
    # Due to the simplicity of the problem, i will just calcualte the number of permutaions instead of counting them.
    f.write(str(math.factorial(n))+'\n')
    for j in itertools.permutations(range(1,n+1)):
        f.write(' '.join(str(val) for val in j)+'\n')