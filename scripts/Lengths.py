#!/usr/bin/env python3

from Bio import SeqIO

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="ifilename",
                  help="fasta file", metavar="FILE")



(options, args) = parser.parse_args()

minqual = 20

handle = open(options.ifilename, "rU")
for record in SeqIO.parse(handle, "fasta"):
	seq = record.seq
	
	print(record.id + "\t" + repr(len(seq)))
handle.close()
