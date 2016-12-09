#from ete2 import NCBITaxa
import argparse
import sys, getopt
import glob, os
import pandas as p
import numpy as np
import re
import operator
import gzip

from collections import defaultdict
from collections import Counter
import logging

#MIN_IDENTITY_TAXA = (0.40,0.40,0.50,0.50,0.50,0.60,0.85)
#50,60,70,80,90,95
MIN_IDENTITY_TAXA = (0.40,0.50,0.60,0.70,0.80,0.90,0.95)
MIN_IDENTITY = 0.40

MAX_MATCHES = 1e100
MIN_FRACTION = 0.9

DEF_DMP_FILE = "/home/chris/native/Databases/nr/FASTA/gi_taxid_prot.dmp"

DEF_LINE_FILE = "/home/chris/native/Databases/nr/FASTA/all_taxa_lineage_notnone.tsv"

def read_blast_input(blastinputfile,lengths): 
    #k191_83_2       gi|973180054|gb|KUL19018.1|     71.2    73      21      0       9       81      337     409     6.6e-24 118.2
    
    #queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore
    
    matches = defaultdict(list)
    gids = Counter()
    nmatches = Counter();
    
    for line in open(blastinputfile):
        line = line.rstrip()
        
        (queryId, subjectId, percIdentity, alnLength, mismatchCount, 
        gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore) = line.split("\t")
        
        m = re.search(r"gi\|(.*?)\|.*", subjectId)
        gid = m.group(1)
        qLength = lengths[queryId]
        alnLength_in_query = abs(int(queryEnd) - int(queryStart)) + 1
        fHit = float(alnLength_in_query)/qLength
        fHit *= float(percIdentity)/100.0
        fHit = min(1.0,fHit)
        #hits[queryId] = hits[queryId] + 1
        if percIdentity > MIN_IDENTITY and nmatches[queryId] < MAX_MATCHES:
            matches[queryId].append((m.group(1),fHit))
            nmatches[queryId] += 1
            gids[gid] +=1


    return (matches, gids.keys())
    
def read_lineage_file(lineage_file): 
    
    mapping = {}
    mapBack = defaultdict(lambda: defaultdict(list))
    for line in open(lineage_file):
        line = line.rstrip()
    
        tokens = line.split("\t")
        (taxaid, domainid, phylumid, classid, orderid, familyid, genusid, speciesid) = tokens
    
        mapping[int(taxaid)]=[domainid, phylumid, classid, orderid, familyid, genusid, speciesid];
        tokens.pop(0)
        for depth in range(6,0,-1):
            if tokens[depth] not in mapBack[depth] and tokens[depth] != 'None':
                for depth2 in range(depth - 1,-1,-1):
                    mapBack[depth][tokens[depth]].append(tokens[depth2])
    return (mapping,mapBack)

def read_query_length_file(query_length_file): 
    
    lengths = {}
    
    for line in open(query_length_file):
        line = line.rstrip()
    
        (queryid, length) = line.split("\t")
    
        lengths[queryid] = float(length)

    return lengths

def map_gids_binary(gids, mapping_file):
    
    t = open(mapping_file, 'r')
    c = 0
    size = long(os.path.getsize(mapping_file))
    mapping = {}
    
    for gis in gids: #binary search #########################################
        c=c+1
        #print gis
        found=False
        offset=0
        chunk=size
        pos=chunk/2
        while found == False and chunk>0:
            chunk = chunk/2
            t.seek(pos)
            t.readline()
            entry = t.readline().split("\t")
            if entry[0]:
                filegi = int(entry[0])
                filetax = int(entry[1].rstrip("\n"))
                
                if filegi == int(gis):
                    answer = filetax
                    found = True
                elif filegi > int(gis):
                    pos = offset +(chunk/2)
                elif filegi < int(gis):
                    offset = offset+chunk
                    pos = pos + (chunk/2)

        if found == False:
            answer = -1

        mapping[gis] = answer

    return mapping

def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("blast_input_file", help="directory with blast 6 matches to taxaid database *.b6")

    parser.add_argument("query_length_file", help="tab delimited file of query lengths")

    parser.add_argument('-g','--gid_taxaid_mapping_file',default=DEF_DMP_FILE,help="mapping from gid to taxaid gzipped")

    parser.add_argument('-l','--lineage_file',default=DEF_LINE_FILE,help="text taxaid to lineage mapping")

    parser.add_argument('-o','--output_dir', type=str, default="output",
        help=("string specifying output directory and file stubs"))

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()
    
    lengths = read_query_length_file(args.query_length_file)
    logging.info("Finished reading lengths file")

    (matches,gids) = read_blast_input(args.blast_input_file,lengths)
    logging.info("Finished reading in blast results file")

    (lineages,mapBack) = read_lineage_file(args.lineage_file)
    logging.info("Finished reading in lineage file")

    mapping = map_gids_binary(gids, args.gid_taxaid_mapping_file)
    logging.info("Finished loading map gids file")

    geneAssign = defaultdict(dict)
    contigLengths = Counter()
    contigAssignDepth = list()
    for depth in range(7):
            contigAssignDepth.append(defaultdict(lambda: Counter()))

    contigGenes = defaultdict(list)
    for gene, matchs in matches.iteritems(): 
        #print str(gene)
        m = re.search(r"(.*)_\d+", gene)
        contig = m.group(1)
        contigLengths[contig] += lengths[gene]
        contigGenes[contig].append(gene)
        
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
        
        
        added_matches = set()
        for (match,fHit) in sorted(matchs, key=lambda x: x[1], reverse=True):

            if mapping[match] > -1:
                tax_id = mapping[match]
                if tax_id not in added_matches:     # Only add the best hit per species
                    added_matches.add(tax_id)
                    hits = lineages[tax_id]
                    for depth in range(7):
                        if hits[depth] != "None":
                            weight = (fHit - MIN_IDENTITY_TAXA[depth])/(1.0 - MIN_IDENTITY_TAXA[depth])
                            weight = max(weight,0.0)
                            if weight > 0.0:
                                collate_hits[depth][hits[depth]] += weight #could put a transform in here


        #import ipdb; ipdb.set_trace()        
        for depth in range(6,-1,-1):
            collate = collate_hits[depth]
            dWeight = sum(collate.values())
        
            sortCollate = sorted(collate.items(), key=operator.itemgetter(1),reverse=True)
            nL = len(collate)
            if nL > 0:
                dP = 0.0
                if dWeight > 0.0:
                    dP = float(sortCollate[0][1])/dWeight
                    
                    if dP > MIN_FRACTION:
                        geneAssign[gene][depth] = (sortCollate[0][0],dP)
                        assignBack = mapBack[depth][sortCollate[0][0]]
                        depth2 = depth -1
                        for assignB in assignBack:
                            geneAssign[gene][depth2] = (assignB,1.0)
                            depth2 -= 1
                        
                        break
                    else:
                        geneAssign[gene][depth] = ('No hits', -1.0)
                else:
                    geneAssign[gene][depth] = ('No hits', -1.0)
            else:
                geneAssign[gene][depth] = ('No hits',-1.0)

    
    with open(args.output_dir+"_genes.csv", "w") as text_file:
        for gene in geneAssign.keys():
            text_file.write('%s'%gene)
            for depth in range(7):
                (assign,p) = geneAssign[gene][depth]
                text_file.write(',%s->%.3f'%(assign,p))
            text_file.write('\n')
            text_file.flush()

    #collate contig predictions 
    #import ipdb; ipdb.set_trace()
    
    contigAssign = defaultdict(dict)
    
    for contig, genes in contigGenes.iteritems():
    
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
    
        for gene in genes:
        
            for depth in range(7):
                (assignhit, genef) = geneAssign[gene][depth] 
            
                if assignhit != 'No hits':
                    collate_hits[depth][assignhit] += lengths[gene]#*genef
        
        for depth in range(7):
            collate = collate_hits[depth]
            dWeight = sum(collate.values())
            sortCollate = sorted(collate.items(), key=operator.itemgetter(1),reverse=True)
            nL = len(collate)
            if nL > 0:
                dP = 0.0
                if dWeight > 0.0:
                    dP = float(sortCollate[0][1])/dWeight
                    contigAssign[contig][depth] = (sortCollate[0][0],dP,sortCollate[0][1])   
                else:
                    contigAssign[contig][depth] = ('No hits',0.,0.)
            else:    
                contigAssign[contig][depth] = ('No hits',0.,0.)

    with open(args.output_dir+"_contigs.csv", "w") as text_file:
        for contig in contigAssign.keys():
            text_file.write('%s,%f'%(contig,contigLengths[contig]))
            for depth in range(7):
                (assign,p,dF) = contigAssign[contig][depth]
                dFN = dF/contigLengths[contig]
                text_file.write(',%s->%.3f->%.3f'%(assign,p,dFN))
            text_file.write('\n')
            text_file.flush()
    
if __name__ == "__main__":
    main(sys.argv[1:])





