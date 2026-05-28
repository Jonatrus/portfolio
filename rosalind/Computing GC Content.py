"""
Problem
The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
"""
# Import libraries
from Bio import SeqIO
import collections

# Instantiate variables
record_id = ""
gc_max = 0

# Loop through the sequences
for record in SeqIO.parse("data/rosalind_gc.txt", "fasta"):
    count = collections.Counter(record.seq) # counter object thtat stoes all the counts of ATGC
    GC = count['G'] + count["C"]
    gc_content = 100 * (GC/len(record.seq))
    # Single pass max tracking so only the maximum is stored in the end
    if gc_content > gc_max:
        record_id = record.id
        gc_max = gc_content

# Save data for easier answer submitting on Rosalind webstie
with open("output/rosalind_gc_output.txt",'w+') as f:
    f.write(f"{record_id}\n{gc_max}")
    
