#!/usr/bin/env python3
import sys

slistfile = open(sys.argv[2])
tigs = {}
for line in slistfile:
    line = line.rstrip("\n")
    tigs[line]=1

asmfile = open(sys.argv[1])
printing=0
for line in asmfile:
    line = line.rstrip("\n")
    if line[0] == '>':
        if line[1:] in tigs:
            printing=1
        else:
            printing=0
    if printing==1: 
        print(line)
