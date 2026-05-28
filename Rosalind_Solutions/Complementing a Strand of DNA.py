with open("rosalind_revc.txt") as f:
    seq = f.read().strip()

myTable = str.maketrans("ACGT", "TGCA")
seq = seq.translate(myTable)
seq = seq[::-1]

with open("output.txt", "w") as f:
    f.write(seq)
