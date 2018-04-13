#!/usr/bin/env python3
# ***************************************************************
# Name:      COG_table.py
# Purpose:   This script integrates with PROKKA or Prodigal and generates Cogs assignments for the protein sequences.
# Version:   0.1
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Created:   2014-01-11
# License:   Copyright (c) 2014 Computational Microbial Genomics Group, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/
import sys
from BCBio import GFF
import argparse
from Bio import Entrez
from collections import defaultdict
import math
import re
def get_records_from_cdd(queries, email):
    # We need CDD accession to COG accession mapping. For this we will use NCBI eutils and parse the returned XML
    # file. For example,
    #
    #   http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=cdd&id=223855
    #
    # returns the following XML record
    #
    #   <eSummaryResult>
    #       <DocSum>
    #           <Id>223855</Id>
    #           <Item Name="Accession" Type="String">COG0784</Item>
    #           <Item Name="Title" Type="String">CheY</Item>
    #           <Item Name="Abstract" Type="String">FOG: CheY-like receiver [Signal transduction mechanisms]</Item>
    #           <Item Name="Status" Type="Integer">0</Item>
    #           <Item Name="LivePssmID" Type="Integer">0</Item>
    #       </DocSum>
    #   </eSummaryResult>
    Entrez.email = email # Always tell ncbi who you are.
    search_result = Entrez.read(Entrez.epost("cdd", id=",".join(queries)))
    records = Entrez.read(Entrez.efetch(db="cdd",
            rettype='docsum',
            webenv=search_result['WebEnv'],
            query_key=search_result['QueryKey']))
    return records

def get_records_from_file(queries, cog_file):
    # Read a simple tsv file with two columns, cddid and cogid respectively.
    with open(cog_file, 'r') as cf:
        records = dict([(row.split('\t')[0], {'Accession': row.split('\t')[1].strip()}) for row in cf.readlines()])
    return records

def usage():
    return '\n'.join([
           'Example usage:',
       '',              
           '\tStep 1: Run PROKKA_XXXXXXXX.faa with rpsblast against the  Cog database',
       '\twith following format:',    
           '\t\t\trpsblast -query PROKKA_XXXXXXXX.faa -db Cog -evalue 0.00001', 
           '\t\t\t-outfmt \"6 qseqid sseqid evalue pident score qstart qend', 
           '\t\t\tsstart send length slen\" -out blast_output.out',
           '',
           '\tStep 2: Run this script to generate the table with marker gene abundance per cluster.:',
           '\t\t\t./COG_table.py -g PROKKA_XXXXXXXX.gff -b blast_output.out -e mail@example.com',
           '\t\t\t -c clustering_gt1000.csv -m marker_genes.txt > scg_table.tsv',
       '',  
           'Refer to rpsblast tutorial: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/',
           ''])

def read_blast_output(blastoutfile): 
    sseq_ids = []
    records = []
    with open(blastoutfile) as in_handle:
        for line in in_handle:
            line_items = line.split("\t")
            qseq = line_items[0]
            sseq = line_items[1]
            pident = line_items[3]
            send = line_items[7]
            sstart = line_items[8]
            slen = line_items[10]

            records.append({'qseqid': qseq,
                            'sseqid': sseq,
                            'pident': float(pident),
                            'send': float(send),
                            'sstart': float(sstart),
                            'slen': float(slen)})

            sseq_ids.append(sseq.split('|')[2])
    return records, sseq_ids

def read_gff_file(gfffile):
    featureid_locations={}
    limits=dict(gff_type=["gene","mRNA","CDS"])
    with open(gfffile) as in_handle:
        for rec in GFF.parse(in_handle, limit_info=limits):
            for feature in rec.features:
                featureid_locations[feature.id] = rec.id
    return featureid_locations

def read_gff_file2(gfffile):
    featureid_locations={}
    limits=dict(gff_type=["gene","mRNA","CDS"])
    with open(gfffile) as in_handle:
        for rec in GFF.parse(in_handle, limit_info=limits):
            for feature in rec.features:
                stoks = feature.id.split("_")
                cid = rec.id + "_" + stoks[1]
                featureid_locations[cid] = feature
    return featureid_locations


def read_markers_file(marker_file):
    # Stores each line of marker_file as an item in a list
    with open(marker_file) as mf:
        return [l.strip() for l in mf.readlines()]

def read_clustering_file(cluster_file):
    # Returns the cluster names and the contig names per cluster
    contigs_per_cluster = defaultdict(list)
    clusters = set()
    with open(cluster_file) as cf:
        for line in cf.readlines():
            line_items = line.strip().split(',')
            cluster = line_items[1]
            contig = line_items[0]
            clusters.add(cluster)
            contigs_per_cluster[cluster].append(contig)
    return list(clusters), contigs_per_cluster


def main(args):
    
    minLength = 100
    #import ipdb; ipdb.set_trace()
    if args.gfffile:
        limits=dict(gff_type=["gene","mRNA","CDS"])
        with open(args.gfffile) as in_handle:
            for rec in GFF.parse(in_handle, limit_info=limits):
                for feature in rec.features:
                    if math.fabs(feature.location.end - feature.location.start) > minLength :                
                        m = re.search(".*_(\d+)", feature.id)
                        idx = m.groups()[0]                        

                        gene_id = rec.id + "_" + idx                        

                        print(gene_id + "," + rec.id + "," + str(feature.location.start) + "," +  str(feature.location.end) + "," + str(feature.strand))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=usage())
   
    parser.add_argument('-g', '--gfffile',
           help=('GFF file generated by e.g. prodigal '
           'only needed if the contig names are not recoverable from the '
           'blast output file.'))

    args = parser.parse_args()

    main(args)
