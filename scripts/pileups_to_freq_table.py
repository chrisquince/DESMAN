#!/usr/bin/env python3
import sys
import re
from Bio import SeqIO


if len(sys.argv)<4:
    print("Error: need at least two pileup files as input\n")
    sys.exit(-1)


nucleotides = ["A","T","C","G","a","t","c","g"]
revcomp = {'a':'T','g':'C','t':'A','c':'G'}

min_base_q = 20
min_base_q += 33 # assume phred33 scale


# get ordered identifiers for contigs
tig_ids = {}
tig_i = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
    tig_ids[record.id] = tig_i
    tig_i += 1

# open up the output file and write header info for the variant sites
desmanout = open(sys.argv[-1],'w')
header = "Contig,Position"
for f in range(1,len(sys.argv)-2):
    header += ","+str(f)+"-A,"+str(f)+"-C,"+str(f)+"-G,"+str(f)+"-T"
desmanout.write(header+"\n")

# open all of the pileup files
pileup_ins = {}
for f in range(1,len(sys.argv)-2):
    pileup_ins[f]=open(sys.argv[f+1])

# write out variant site counts one contig at a time
curtig = -1
inspattern = re.compile("\d+")
while curtig < tig_i:
    curtig += 1
    allsites = {}
    counts = {}
    for f in range(1,len(sys.argv)-2):
        counts[f] = {'a':{},'c':{},'g':{},'t':{}}
        while 1:
            curpos = pileup_ins[f].tell()
            line = pileup_ins[f].readline()
            if line == '':
                break
            if line[0] == "#":
                continue
            dd = line.rstrip("\n").split('\t')

            # if we've scanned past the current contig then reset the file pointer
            if tig_ids[dd[0]] > curtig:
#                print "hit contig " + dd[0] + " seeking to " + str(curpos)
                pileup_ins[f].seek(curpos)
                break
            curpos = pileup_ins[f].tell() # otherwise update the file pointer

            dd[2] = dd[2].upper()
            freqs = {"A":0,"T":0,"C":0,"G":0}

            index = -1
            bindex = -1
            while index < len(dd[4])-1:
                index += 1
                i=dd[4][index]
                if i == "." or i == ",":
                    bindex += 1
                    if ord(dd[5][bindex]) >= min_base_q:
                        freqs[dd[2]] = freqs[dd[2]] + 1
                    continue

                if i == '^':
                    index += 1 # skip the next character
                    continue

                if i == '+' or i == '-':
                    m = inspattern.match(dd[4],index+1)
                    index += int(m.group(0)) + len(m.group(0)) + 1

                if i in nucleotides:
                    bindex += 1
                    if ord(dd[5][bindex]) >= min_base_q:
                         freqs[i.upper()] = freqs[i.upper()]+1

            sitename = dd[0]+","+dd[1]
            allsites[sitename]=1
            counts[f]['a'][sitename] = freqs['A']
            counts[f]['c'][sitename] = freqs['C']
            counts[f]['g'][sitename] = freqs['G']
            counts[f]['t'][sitename] = freqs['T']
               
    vsites = [] # list of variant sites & counts
    for c in allsites:
        nonzero=0    
        for nt in counts[1]:
            for f in counts:
                if(c in counts[f][nt] and counts[f][nt][c] > 0): 
                    nonzero+=1
                    break
        # ensure at least two different nucleotides have been observed in the data            
        if(nonzero<2): continue
        vsites.append(c)

    # BUG WORKAROUND: only write out sites if the contig has at least 2
    if len(vsites) < 2: continue
    for v in vsites:
        line=v
        for f in counts:
            if not v in counts[f]['a']: 
                counts[f]['a'][v] = 0
                counts[f]['c'][v] = 0
                counts[f]['g'][v] = 0
                counts[f]['t'][v] = 0
            line += ","+str(counts[f]['a'][v])
            line += ","+str(counts[f]['c'][v])
            line += ","+str(counts[f]['g'][v])
            line += ","+str(counts[f]['t'][v])
        desmanout.write(line+"\n")


