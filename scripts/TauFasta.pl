#!/usr/bin/perl

use strict;

my @files = <Select/COG*.gfa>;

my %hashBase = (
        "a" => '0',
        "c" => '1',
        "g" => '2',
	    "t" => '3',
	    "-" => '4');

foreach my $fastaFile(@files){
    $fastaFile=~/(.*).gfa/;
    my $stub = $1;
    
    $stub=~/Select\/(.*)/;
    my $contig = $1;
    
    print "$stub\n";

    if( -e "${stub}_R.gfn" ){
    	print "Reverse\n";
	    $fastaFile = "${stub}_R.gfn";
    }

    open(FILE, $fastaFile) or die "Can't open $fastaFile\n";
    
    my @Seq = ();
    my @id       = ();

    my $count = 0;
    my $seq = "";

    while(my $line = <FILE>){
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
    close(FILE);

    $Seq[$count - 1] = $seq;
    my $refTotal = $count - 1;
    my $seqTotal = $count;
    
    my @allHits = ();
    for(my $i = 0; $i < $seqTotal; $i++){
        my @hit = split(//,$Seq[$i]);
        push(@allHits,\@hit);
    }
    

    my @compHits = ();
    my @pos = ();
    my $length = scalar(@{$allHits[0]});
    
    my $nHits = scalar(@allHits);
    my $nVar = 0;
    my $p = 0;
    my $l = 0;
    #print "@{$allHits[0]}\n";
    while($l < $length){
    
        if($allHits[0][$l] ne '-'){
            #print "$allHits[0][$l]\n";
            my @col = ();
            my @var = ();
            for(my $n = 0; $n < 4; $n++){
                $col[$n] = 0;
            }
            
            my $nA = 0;
            for(my $n = 1; $n < $seqTotal; $n++){
                my $d = $hashBase{$allHits[$n][$l]};
                if($d < 4){    
                    $col[$d]++;
            
                    push(@var,$d);
                    $nA++;
                }
            }  
            
            my $maxn = 0;
            for(my $n = 0; $n < 4; $n++){
                if($col[$n] > $maxn){
                    $maxn = $col[$n];
                }
            }
            #print "$p @var\n";
            if($maxn < $nA){
                #print "Variant $maxn $nA\n";
                push(@compHits,\@var);
                push(@pos,$p);
                $nVar ++;
            }
            
            $p++;
        }
        
        $l++;
    }
    
#    shift(@id);
    
    my @refs = @id;
    
    close(FILE);  
    
    open(TFILE, ">Select/${contig}_tau.csv") or die "Can't open >Select/${contig}_tau.csv\n";
    
    
    print TFILE "Contig,Pos,";
    
    my @sampleTags = ();
    foreach my $sample(@refs){
        my @tags = ("${sample}-A","${sample}-C","${sample}-G","${sample}-T");
        push(@sampleTags,@tags);
    }
    
    my $gString = join(",",@sampleTags);
    print TFILE "$gString\n";
     
    for(my $v = 0; $v < $nVar;$v++){
        
        print TFILE "$contig,$pos[$v],";
    
        my @vArray = ();
        
        for(my $n = 0; $n < $nHits; $n++){
            for(my $i = 0; $i < 4; $i++){
                if($i == $compHits[$v][$n]){
                    push(@vArray,1);
                }
                else{
                    push(@vArray,0);
                }
            } 
        }
        
        my $vString = join(",",@vArray);
        print TFILE "$vString\n";
    }   
        
    close(TFILE); 
}
