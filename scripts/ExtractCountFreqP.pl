#!/usr/bin/env perl

use strict;

my $posFile = $ARGV[0];
my $inputDir = $ARGV[1];

my $minLength = 0;

my %hashStart = ();
my %hashEnd = ();
my %hashCogs = ();
my %hashCogStart = ();
my %hashCogEnd = ();
open(FILE, $posFile) or die "Can't open $posFile\n";

while(my $line = <FILE>){
    chomp($line);
#COG0085,contig-50001822,0,2913,contig-50001822_1    
    my @tokens = split(/,/,$line);
    my $contig = $tokens[1];
    my $cog = $tokens[0];
    my $start = $tokens[2];
    my $end = $tokens[3];

    if($end - $start > $minLength){
    $hashCogStart{$cog} = $start;
    
    $hashCogEnd{$cog} = $end;	    

    if($hashStart{$contig} eq undef){
        my @start = ($start);
        my @end = ($end);
	         
        $hashStart{$contig} = \@start;
        $hashEnd{$contig} = \@end;
    }
    else{
        push(@{$hashStart{$contig}},$start);
        push(@{$hashEnd{$contig}},$end);
    }
    
    $hashCogs{$contig}{$start} = $cog;
    }
}

foreach my $contig(sort keys %hashStart){
        my @start = @{$hashStart{$contig}};
        my @end =  @{$hashEnd{$contig}};
        
        @start = (sort {$a <=> $b} @start);
        @end = (sort {$a <=> $b} @end);
        $hashStart{$contig} = \@start;
        $hashEnd{$contig} = \@end;
}

my @countFiles = <${inputDir}/*.cnt>; 
my @samples = ();

my %hashCogPosSample = ();
foreach my $countFile(@countFiles){
	open(FILE,$countFile) or die "Can't open $countFile\n";

	$countFile =~ /${inputDir}\/(.*).cnt/;
	
	my $sample = $1;
	push(@samples,$sample);
	while(my $line = <FILE>){
		chomp($line);

		my @tokens = split(/\t/,$line);

		my $contig=shift(@tokens);
		my $pos=shift(@tokens);

		shift(@tokens);
		shift(@tokens);
		shift(@tokens);
		my @var = ();
		for(my $i = 0; $i < 4;$i++){
			my $tok = shift(@tokens);
			my @tokens2 = split(/:/,$tok);
			push(@var,$tokens2[1]);
		}
	

  		if($hashStart{$contig} ne undef){
        		#print "$contig $pos\n";
        		my $idx = 0;
        		my @start = @{$hashStart{$contig}};
        		my @end =  @{$hashEnd{$contig}};
        		my $nLen = scalar(@start);
        		#print "@start\n";
        		#print "@end\n";
        		while($pos >= $end[$idx] && $idx < $nLen){
         			#   print "$idx $pos $end[$idx]\n";
            			$idx++;
        		}


        		if($pos >= $start[$idx] && $pos < $end[$idx]){
            			my $cog = $hashCogs{$contig}{$start[$idx]};
				
				$hashCogPosSample{$cog}{$pos}{$sample} = \@var;
			}
		}
		
	} 
}

my @sampleTags = ();
foreach my $sample(@samples){
    my @tags = ("${sample}-A","${sample}-C","${sample}-G","${sample}-T");
    push(@sampleTags,@tags);
}

my $sString = join(",",@sampleTags);
print "Cog,Position,${sString}\n";


foreach my $cog(sort keys %hashCogStart){
	my $start = $hashCogStart{$cog};
	my $end = $hashCogEnd{$cog};
		
	for(my $pos = $start; $pos < $end; $pos++){
		my @sFreqs = ();
		foreach my $sample(@samples){
			my @freqs = (0,0,0,0);
			if($hashCogPosSample{$cog}{$pos}{$sample} ne undef){
				@freqs = @{$hashCogPosSample{$cog}{$pos}{$sample}};
			}
			push(@sFreqs,@freqs);
		}

		my $sString = join(",",@sFreqs);
		my $lpos = $pos - $start;
		print "$cog,$lpos,$sString\n";
	} 
}

