#!/usr/bin/perl

use strict;
#use warnings;

my %hashContigCogs = ();

while(my $line = <STDIN>){
    chomp($line);
    
    my @tokens = split(/,/,$line);
        
    #$tokens[0]=~/(.*)_\d+/;
        
    my $contig = $tokens[4];
        
    my $cog = $tokens[0];
        
    $hashContigCogs{$cog} = $contig;
    
}

my $fastaFile = shift(@ARGV);
my %hashSeq = {};
my @Seq = ();
my @id  = ();

my $count = 0;
#Read in reference fasta file
open(FILE, $fastaFile) or die "Can't open $fastaFile\n";

my $seq = "";

while(my $line = <FILE>){ 
    chomp($line);
    
    if($line =~ />(.*?)\s.*/){
	    
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
close(FILE);

$Seq[$count - 1] = $seq;

for(my $i = 0; $i < $count; $i++){
    $hashSeq{$id[$i]} = $Seq[$i];
}

for my $cog (keys %hashContigCogs){
    my $contig = $hashContigCogs{$cog};
    my $seq = $hashSeq{$contig};
    
    open(FILE, ">SCGs/${cog}.fa");
    
    print FILE ">contig\n$seq\n";
    
    close(FILE);
     
}
