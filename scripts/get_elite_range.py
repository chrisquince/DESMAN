#!/usr/bin/env python
import sys

# read the sequence name to ID map
psids = open(sys.argv[1])
ids = {}
for line in psids:
    d = line.rstrip("\n").split("\t")
    d2 = d[1].split("/")
    ids[d2[0]]=d[0]

# read the list of elite gene regions and write out in bed format
for f in sys.argv[2:]:
    psaln = open(f)
    line = psaln.readline()
    d = line[1:-1].split(".")
    chrname = ids[d[0]].split(" ")
    print(chrname[0] + "\t" + d[2] + "\t" + d[3])
