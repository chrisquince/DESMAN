#!/usr/bin/env python3
import sys
from Bio import SeqIO

seq_index = SeqIO.index(sys.argv[1], 'fasta')
tau_star_file = open(sys.argv[2])
etaS_file = open(sys.argv[3])
outbase = sys.argv[4]

head = etaS_file.readline()
hl = head.split(",")
strainseqs = [{} for x in range(len(hl)-1)]
for line in etaS_file:
    d = line.rstrip("\n").split(",")
    if(len(d[0])==0): continue
    if not d[0] in seq_index: continue
    seq = str(seq_index[d[0]].seq)

    for k in range(0,len(d)-1):
        if(float(d[k+1])>=0.5):
            strainseqs[k][d[0]] = seq

acgt = "ACGT"
for line in tau_star_file:
    d = line.rstrip("\n").split(",")
    for k in range(0,len(strainseqs)):
        if not d[0] in strainseqs[k]: continue
        for n in range(0,3):
            if(float(d[k*4+n+2])>0.5):  
                p = int(d[1])
                strainseqs[k][d[0]] = strainseqs[k][d[0]][:p] + acgt[n] + strainseqs[k][d[0]][(p+1):]

for k in range(0,len(strainseqs)):
    strainout = open(outbase + ".strain." + str(k) + ".fa", "w")
    for s in strainseqs[k]:
        strainout.write(">" + s + "\n" + strainseqs[k][s] + "\n")
