#!/usr/bin/env nextflow
/**
 * Usage: desmanflow.nf --assembly=/path/to/assembly.fa --speciescontigs=contig.list --inputreads=/path/to/read_directory 
 *
 * Optional parameters:
 * --straincount=<number>  sets the number of strains
 * --output=<path> sets the output path
 * Reads are expected to be paired gzip'd FastQ with extensions .r1.fq.gz and .r2.fq.gz
 */

params.straincount = 0
params.output = "out"
params.r12extglob = ".r{1,2}.fq.gz"

speciescontigslist = file(params.speciescontigs)
assembly4elite = file(params.assembly)
assembly4map = file(params.assembly)
assembly4pileup = file(params.assembly)
assembly4elvar = file(params.assembly)
assembly4straintigs = file(params.assembly)

// check for input files
//Channel.fromFilePairs(params.inputreads+"/*."+params.r12extglob).
                
// get the input files
readfiles = Channel.fromFilePairs(params.inputreads+"/*"+params.r12extglob)
Channel.fromFilePairs(params.inputreads+"/*"+params.r12extglob).ifEmpty('\n==\n==> ERROR: no paired FastQ files were found in the directory '+params.inputreads+' matching the file extension glob '+params.r12extglob+'\n==\n').println()


// map reads to coassembly reference. 
// this might be skipped in a future workflow since it will likely already be done to make the species bins
process mapReads {
    input:
    set x, file(x1), file(x2) from readfiles
    file 'assembly.fa' from assembly4map

    output:
    file '*.sort.bam' into bamfiles

    """
    bwa index assembly.fa
    bwa mem assembly.fa x1 x2 | samtools view -q40 -S -b - | samtools sort -o - - > ${x}.sort.bam
    """
}

// uses phylosift to find marker genes in the assembly
process findEliteGenes {
    input:
    file 'species.contigs' from speciescontigslist
    file 'assembly.fa' from assembly4elite

    output:
    file 'elites.bed' into elitebed, elitebed2
    file 'species_contigs.fa' into speciesfa

    """
    ${DESMANHOME}/scripts/extract_species_contigs.py assembly.fa species.contigs > species_contigs.fa
    WORKDIR=`pwd`
    cd ${DESMANHOME}/external/
    tar xjf phylosift_v1.0.1.tar.bz2
    cd \$WORKDIR
    ${DESMANHOME}/external/phylosift_v1.0.1/phylosift search --besthit --isolate species_contigs.fa
    ${DESMANHOME}/external/phylosift_v1.0.1/phylosift align --besthit --isolate species_contigs.fa
    ${DESMANHOME}/scripts/get_elite_range.py PS_temp/species_contigs.fa/blastDir/lookup_ID.1.tbl PS_temp/species_contigs.fa/alignDir/DNGNGWU*.codon.updated.1.fasta > elites.bed 
    """
}

// makes pileups for the marker gene ranges
process elitePileups {
    input:
    file('elites.bed') from elitebed
    file('assembly.fa') from assembly4pileup
    each bam from bamfiles

    output:
    file('*.pileup') into pileup
    file '*.sort.bam' into q40bam

    """
    X=`basename ${bam} .bam`
    samtools faidx assembly.fa
    samtools index ${bam}
    cp ${bam} .
    samtools mpileup -l elites.bed -f assembly.fa ${bam} > \$X.pileup 
    """
}

// combine all the read pileups to a single input for the next step
allpileup = pileup.reduce([]){acc, it -> acc << it}

// call high quality variants from the counts
process callEliteVariants {
    input:
    file '*.pileup' from allpileup
    file 'assembly.fa' from assembly4elvar

    output:
    set file('dfreqssel_var.csv'),file('dfreqstran_df.csv') into desmanfreqs, desmanfreqs2
    
    """
    ${DESMANHOME}/scripts/pileups_to_freq_table.py assembly.fa *.pileup desmanfreqs.csv
    ${DESMANHOME}/desman/Variant_Filter.py desmanfreqs.csv -o dfreqs -p
    """
}

// estimate the true number of strains
process estimateStrainCountDesman {
    input:
    set file('outputsel_var.csv'),file('outputtran_df.csv') from desmanfreqs
    each g from 1,2,3,4,5,6,7,8,9,10
    each repid from 1,2,3,4,5
    
    output:
    file('fit_*.txt') into desman_dic
    
    when:
    params.straincount == 0

    script:
    """
    ${DESMANHOME}/bin/desman outputsel_var.csv -e outputtran_df.csv -o cluster_${g}_${repid} -r 1000 -i 100 -g $g -s $repid > cluster_${g}_${repid}.out
    cp */fit.txt fit_${g}_${repid}.txt
    """
}

// collate all results from strain estimation
alldics = desman_dic.reduce([]){acc, it -> acc << it}

// select the best strain count
process plotDev {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    file '*' from alldics

    output:
    file 'Dev.pdf' into devpdf

    when:
    params.straincount == 0

    script:
    """
    cat fit*.txt | cut -d"," -f2- > Dev_no_header.csv
    cat <(echo H,G,LP,Dev) Dev_no_header.csv > Dev.csv
    mkdir -p ~/.Rlibs
    ${DESMANHOME}/scripts/PlotDev.R -l Dev.csv -o Dev.pdf
    """
}

devpdf.subscribe { 
    println("\nStrain count evaluation complete. Please inspect the file ${params.output}/Dev.pdf and select the lowest number of strains at the knee in the curve.\nThen re-run desmanflow.nf with the following command, specifying the chosen number of strains with --straincount:\n")
    println(workflow.commandLine + " --straincount=<number> ") 
}

process desman {
    input:
    set file('dfreqssel_var.csv'),file('dfreqstran_df.csv') from desmanfreqs2
    
    output:
    set file('Gamma_star.csv'), file('Eta_star.csv') into desman
    
    when:
    params.straincount != 0

    script:
    """
    ${DESMANHOME}/bin/desman dfreqssel_var.csv -e dfreqstran_df.csv -o cluster -r 1000 -i 100 -g ${params.straincount}
    mv cluster/Gamma_star.csv cluster/Eta_star.csv .
    """
}


allq40bams = q40bam.reduce([]){acc, it -> acc << it}


process straintigs {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set file('Gamma_star.csv'), file('Eta_star.csv')  from desman    
    file('species.fa') from speciesfa
    file('assembly.fa') from assembly4straintigs
    file('*') from allq40bams
    file('elite.bed') from elitebed2

    output:
    file('straincluster*') into straintigs

    """
    ${DESMANHOME}/scripts/Lengths.py -i species.fa > species_contigs.bed
    perl -p -i -e "s/\t/\t1\t/g" species_contigs.bed

    samtools faidx assembly.fa
    for file in *.bam
    do
        bname=`basename \$file .bam`
        samtools index \$file
        samtools mpileup -l species_contigs.bed -f assembly.fa \$file > \$bname.pileup         
    done

    ${DESMANHOME}/scripts/pileups_to_freq_table.py assembly.fa *.pileup contigfreqs.csv
    rm *.pileup
    ${DESMANHOME}/desman/Variant_Filter.py contigfreqs.csv -m 0.0 -v 0.03
    ${DESMANHOME}/scripts/CalcGeneCov.py contigfreqs.csv > contig_cov.csv

    cut -f 1 elite.bed | sort | uniq > core_genes.txt
    ${DESMANHOME}/scripts/CalcDelta.py contig_cov.csv core_genes.txt cluster_core
    ${DESMANHOME}/bin/desman outputsel_var.csv -e outputtran_df.csv -o straincluster -r 1000 -i 100 -g ${params.straincount}
    ${DESMANHOME}/desman/GeneAssign.py cluster_coremean_sd_df.csv straincluster/Gamma_star.csv contig_cov.csv straincluster/Eta_star.csv -m 20 -v outputsel_var.csv -o straincluster --assign_tau > cluster.cout
    ${DESMANHOME}/scripts/write_strain_fasta.py species.fa straincluster_tau_star.csv strainclusteretaD_df.csv straincluster
    """
}

