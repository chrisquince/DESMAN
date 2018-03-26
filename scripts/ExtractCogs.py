#!/usr/bin/env python
import sys
from BCBio import GFF
import argparse
from Bio import Entrez
from collections import defaultdict

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
     
    RPSBLAST_SCOVS_THRESHOLD = args.scovs_threshold
    RPSBLAST_PIDENT_THRESHOLD = args.pident_threshold

    records, sseq_ids = read_blast_output(args.blastoutfile)

    #import ipdb; ipdb.set_trace()

    if args.cdd_cog_file:
        # Retrieve the cog accession number from file
        cogrecords = get_records_from_file(sseq_ids, args.cdd_cog_file)
    else:
        sys.exit("Need CDD to COG file")   
    
    if args.gfffile:
        featureid_locations = read_gff_file2(args.gfffile)
    
    features_per_contig = defaultdict(list)
    for record_d in records:
        pident_above_threshold = record_d['pident'] >= RPSBLAST_PIDENT_THRESHOLD

        # A certain fraction of the cog should be covered to avoid the same cog 
        # to be counted twice in the case when a cog is split across two or more contigs.
        alignment_length_in_subject = abs(record_d['send'] - record_d['sstart']) + 1
        percent_seq_covered = (alignment_length_in_subject / record_d['slen']) * 100.0
        seq_cov_above_threshold =  percent_seq_covered >= RPSBLAST_SCOVS_THRESHOLD
        
        if pident_above_threshold and seq_cov_above_threshold:
            cog_accession = cogrecords[record_d['sseqid'].split('|')[2]]['Accession']
            queryloc = featureid_locations[record_d['qseqid']]
            print(record_d['qseqid']+","+cog_accession + "," + str(queryloc.location.start) + "," +  str(queryloc.location.end) + "," + str(record_d['send']) + "," + str(record_d['sstart']) + "," + str(
queryloc.strand))


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('-b', '--blastoutfile', required=True,
           help=('Output of rpsblast run, assumed to be in tabular format whith '
               'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen. '
               'The contigs ids are assumed to be recoverable by removing the last underscore '
               'and the characters following it from the qseqid column.' ))
   parser.add_argument('-s', '--scovs-threshold', type=float, default=50.0,
           help='Threshold covered in percent, default=50.0')
   parser.add_argument('-p', '--pident-threshold', type=float, default=0.0,
           help='Threshold identity in percent, default=0.0')
   parser.add_argument('--cdd_cog_file',
           help = ('Supply a cdd to cog mapping file in a tsv format '
           'to take precedence over eutils fetching of name. '
           'Useful if running this script in parallel, since '
           'NCBI eutils has a limit on the number of requests per '
           'time unit you can make.'))
   parser.add_argument('-g', '--gfffile',
           help=('GFF file generated by e.g. prodigal '
           'only needed if the contig names are not recoverable from the '
           'blast output file.'))

   args = parser.parse_args()

   main(args)
