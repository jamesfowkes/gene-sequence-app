import sys
import argparse
import pickle
import os

from Bio import Entrez, SeqIO

Entrez.email = "stxms13@nottingham.ac.uk"

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def entrez_get_fasta_seqrecord(gid):
    try:
        handle = Entrez.efetch(db="nucleotide", id=gid, rettype="gb", retmode="text")
        return SeqIO.read(handle,"genbank")
    except:
        return None

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("gid", help="The GenBank ID to fetch")
    parser.add_argument("-l", "--load_pickle", nargs='?', const="use_default", help="Specify to load from pickle file")
    parser.add_argument("-p", "--pickle", nargs='?', const="use_default", help="Output pickled object to this file (defaults to {gid}.pickle)")

    args = parser.parse_args()

    if args.load_pickle:
        if args.load_pickle == "use_default":
            args.load_pickle = "{}.pickle".format(args.gid)

    if args.pickle == "use_default":
        args.pickle = "{}.pickle".format(args.gid) 

    return args

if __name__ == "__main__":

    args = get_args()

    if args.load_pickle:
        with open(__location__+ "/" + args.load_pickle, 'rb') as f:
            fasta_data = pickle.load(f)
    else:
        fasta_data = entrez_get_fasta_seqrecord(args.gid)

    if args.pickle:
        with open(__location__+ "/" + args.pickle, 'wb') as f:
            pickle.dump(fasta_data, f)
    else:
        sys.stdout.write(str(fasta_data))
