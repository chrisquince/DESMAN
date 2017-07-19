#!/usr/bin/perl

use strict;


my $cogFile = $ARGV[0];

my %hashContig = {};

open(COGFILE, $cogFile) or die "Can't open $cogFile\n";

while(my $line = <COGFILE>){
    chomp($line);
    
    my @tokens = split(",",$line);
    #COG0552,k141_866352,4207,5161,k141_866352_7,-1
    
    $hashContig{$tokens[0]} = $tokens[4];
    
}

close(COGFILE);

my $header = <STDIN>;

print "$header";

while(my $line = <STDIN>){
    chomp($line);

    my @tokens = split(/,/,$line);

    my $cog = shift(@tokens);
    

    my $contig = $hashContig{$cog};

    if ($contig ne undef){
        unshift(@tokens,$contig);

        my $tString = join(",",@tokens);

        print "$tString\n";
    }
}


