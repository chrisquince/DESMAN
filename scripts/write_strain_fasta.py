#!/usr/bin/env python
import sys
from Bio import SeqIO
from argparse import ArgumentParser


def extract_strain_seqs(etaS_file, seq_index):
    strainseqs = []
    with open(etaS_file, 'r') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                hl = line.split(",")
                strainseqs = [{} for x in range(len(hl) - 1)]
            d = line.rstrip("\n").split(",")
            if (len(d[0]) == 0): continue
            if not d[0] in seq_index: continue
            seq = str(seq_index[d[0]].seq)
            for k in range(0, len(d) - 1):
                if (float(d[k + 1]) >= 0.5):
                    strainseqs[k][d[0]] = seq
    return strainseqs


def update_strain_seqs(tau_star_file, strainseqs):
    acgt = "ACGT"
    with open(tau_star_file, 'r') as fh:
        for line in fh:
            d = line.rstrip("\n").split(",")
            # Iterate all strains
            for strain_index in range(0, len(strainseqs)):
                # Skip if this gene is not predicted to be in the strain
                gene_id = d[0]
                if not gene_id in strainseqs[strain_index]: continue
                for n in range(0, 3):
                    val = float(d[strain_index * 4 + n + 2])
                    if val > 0.5:
                        # get position in gene
                        p = int(d[1])
                        # Update nucleotide position for strain
                        strainseqs[strain_index][gene_id] = strainseqs[strain_index][gene_id][:p] + acgt[n] + strainseqs[strain_index][gene_id][(p + 1):]
    return strainseqs


def write_strains(strainseqs, outbase):
    for k in range(0, len(strainseqs)):
        with open("{}.strain.{}.fa".format(outbase, k), 'w') as fh:
            for s in sorted(strainseqs[k]):
                seq = strainseqs[k][s]
                fh.write(">{s}\n{seq}\n".format(s=s, seq=seq))


def main(args):
    seq_index = SeqIO.index(args.fasta, 'fasta')
    tau_star_file = args.tau_star_file
    etaS_file = args.etaS_file
    outbase = args.outbase
    strainseqs = extract_strain_seqs(etaS_file, seq_index)
    # Update strains
    strainseqs = update_strain_seqs(tau_star_file, strainseqs)
    # Write strains
    write_strains(strainseqs, outbase)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("fasta", type=str,
                        help="Gene nucleotide fasta file")
    parser.add_argument("tau_star_file", type=str,
                        help="Prediction for strain haplotypes ('_tau_star.csv' output from GeneAssign.py)")
    parser.add_argument("etaS_file", type=str,
                        help="Prediction for error transition matrix ('etaD_df.csv' output from GeneAssign.py")
    parser.add_argument("outbase", type=str,
                        help="Basename for output files")
    args = parser.parse_args()
    main(args)

