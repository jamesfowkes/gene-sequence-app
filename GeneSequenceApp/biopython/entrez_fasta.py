from Bio import Entrez, SeqIO

Entrez.email = "stxms13@nottingham.ac.uk"

def entrez_get_fasta(gid):
    handle = Entrez.efetch(db="nucleotide", id=gid, rettype="gb", retmode="text")
    return SeqIO.read(handle,"genbank")
    