from Bio import Entrez
from Bio import SeqIO

# NCBI requires email
Entrez.email = "jonny.kobrus@outlook.com"

accession = "PRJNA891285"

handle = Entrez.efetch(
    db="nucleotide",
    id=accession,
    rettype="gb",
    retmode="text"
)

record = SeqIO.read(handle, "genbank")