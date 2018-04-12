#!/usr/bin/perl

use strict;
#use warnings;

my %hashSCGs = ();

my $scgFile = shift(@ARGV);

open(FILE, $scgFile) or die "Can't open $scgFile\n";

while(my $line = <FILE>){
    chomp($line);
    
    $hashSCGs{$line} = 1;
}   

close(FILE);

my %hashContigCogs = ();

while(my $line = <STDIN>){
    chomp($line);
    
    my @tokens = split(/,/,$line);
        
    $tokens[0]=~/(.*)_\d+/;
    my @contig_items = split(/_/,$tokens[0]);
    pop @contig_items;
    my $contig = join("_", @contig_items);
        
    my $cog = $tokens[1];
    #k99_40717_1,COG2207,1,370,9.0,116.0,1        
    if($hashSCGs{$cog} == 1){
        if($hashContigCogs{$cog} eq undef){
            my @temp = ($contig,$tokens[2],$tokens[3],$tokens[0],$tokens[6]);
            my @temp2 = ();
            push(@temp2,\@temp);
            $hashContigCogs{$cog} = \@temp2; 
        }
        else{
            my @temp = ($contig,$tokens[2],$tokens[3],$tokens[0],$tokens[6]);
            push(@{$hashContigCogs{$cog}},\@temp);
        }
    }
}

foreach my $cog (sort keys %hashContigCogs){
    my @contigs = @{$hashContigCogs{$cog}};
    
    my $nM = scalar(@contigs);
    #print "$cog,$nM\n";
    if($nM == 1){
        my @hit = @{$contigs[0]};
        my $hstring = join(",",@hit);
        print "$cog,$hstring\n";
    }
    
}






