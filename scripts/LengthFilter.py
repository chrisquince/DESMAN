#!/usr/bin/env python3

from Bio import SeqIO
import sys
import argparse

def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("inputfile", help="fasta file",metavar="FILE")

    parser.add_argument('-m','--minlength', type=float, default=1000.0, help=("minimum coverage for sample to be included"))

    args = parser.parse_args()    

    filtered = []
    handle = open(args.inputfile, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        seq = record.seq
    
        if len(seq) > args.minlength:
            filtered.append(record)
 
    SeqIO.write(filtered, sys.stdout, "fasta")

if __name__ == "__main__":
    main(sys.argv[1:])
