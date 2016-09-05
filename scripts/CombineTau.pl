#!/usr/bin/perl

use strict;

my @tau_files = <Select/*_tau.csv>;

my %hashCogVar = ();

my @store_ids = ();
my $bFirst = 1;

for my $tau_file(@tau_files){
    open(FILE, $tau_file) or die "Can't open $tau_file\n";
    
    my $line = <FILE>;
    chomp($line);
    
    my @tokens = split(/,/,$line);
    
    shift(@tokens); shift(@tokens);
    my @tids = ();
    my $tLen = scalar(@tokens);
    for(my $i = 0; $i < $tLen; $i+=4){
        push(@tids,$tokens[$i]);
    }
    my $nG = $tLen/4;
    
    my @ids = ();
    for(my $i = 0; $i < $nG; $i++){
        $tids[$i]=~/(.*)_.*?/;
        push(@ids,$1);
    }
    
    if($bFirst == 1){
        @store_ids = @ids;
        $bFirst = 0;
    }
    
    while($line = <FILE>){
        chomp($line);
    
        my @tokens = split(/,/,$line);
        my $cog = shift(@tokens);
        my $pos = shift(@tokens);
        my %tHash = ();
        for(my $i = 0; $i < $nG; $i++){
            my $nS = $i*4;
            my $nE = ($i + 1)*4 - 1;
            my @var = @tokens[$nS..$nE];
            #print "$cog,$pos,$ids[$i],@var\n";
            $tHash{$ids[$i]} = \@var;  
        }
        $hashCogVar{$cog}{$pos} = \%tHash; 
    }
}

my @sampleTags = ();
foreach my $sample(@store_ids){
    my @tags = ("${sample}-A","${sample}-C","${sample}-G","${sample}-T");
    push(@sampleTags,@tags);
}
    
my $gString = join(",",@sampleTags);
print "Contig,Position,$gString\n";

foreach my $cog (sort keys %hashCogVar){
    foreach my $pos(sort {$a <=> $b} keys %{$hashCogVar{$cog}}){
        my %vHash = %{$hashCogVar{$cog}{$pos}};
        my @allVar = ();
        foreach my $sample(@store_ids){
            #print "$sample\n";
            my @var = @{$vHash{$sample}};
            push(@allVar,@var);
            
        }
        my $vString = join(",",@allVar);
        print "$cog,$pos,$vString\n";
    }
}
