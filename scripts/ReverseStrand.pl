#!/usr/bin/perl

use strict;

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

my $cogFile = $ARGV[0];

my %hashStrand = {};

open(COGFILE, $cogFile) or die "Can't open $cogFile\n";

while(my $line = <COGFILE>){
    chomp($line);
    
    my @tokens = split(",",$line);
    #print "$tokens[0] $tokens[6]\n";
    $hashStrand{$tokens[0]} = $tokens[5];
    
}

close(COGFILE);

my @alnFiles = <Select/*.gfa>;

for my $alnFile(@alnFiles){
    print "$alnFile\n";
    
    $alnFile=~/Select\/(.*).gfa/;
    my $cog = $1;
    #print "$alnFile $cog $hashStrand{$cog}\n";
    if($hashStrand{$cog} == -1){
        open(AFILE,$alnFile) or die "Can't open $alnFile\n";
        my @Seq = ();
        my @id       = ();

        my $count = 0;
        my $seq = "";

        while(my $line = <AFILE>){
            chomp($line);

            if($line =~ />(.*)/){

                $id[$count] = $1;

                if($seq ne ""){
                    $Seq[$count - 1] = $seq;

                    $seq = "";
                }

                $count++;
            }
            else{
                $seq .= $line;
            }
        }
        close(AFILE);

        $Seq[$count - 1] = $seq;
        my $seqTotal = $count;
    
        my $revFile = ">Select/${cog}_R.gfa";
    
        open(RFILE,$revFile) or die "Can't open $revFile\n";
        
        for(my $i = 0; $i < $seqTotal; $i++){
            my $rev_Seq = reverse_complement($Seq[$i]);
            print RFILE ">$id[$i]\n$rev_Seq\n";
        }
    
        close(RFILE);
    }
    
}
