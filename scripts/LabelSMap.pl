#!/usr/bin/env perl

my $tFile = $ARGV[0];

my $zFile = $ARGV[1];

my $minZ = $ARGV[2];
my $minCount = $ARGV[3];

my @t = ();
my $maxt = 0;

my %hashCluster = {};

open(FILE, $tFile) or die "Can't open $xyFile";

while($line = <FILE>){
  
  chomp($line);
  
  my @tokens = split(/,/,$line);

  my $name = $tokens[0];
  my $cluster = $tokens[1];

  $hashCluster{$name} = $cluster;
  #print "$name $cluster\n";
  if($cluster > $maxt){
    $maxt = $cluster;
  }
}

close(FILE);

open(FILE, $zFile) or die "Can't open $zFile";

$line = <FILE>;
chomp($line);

my @tokens = split(/\t/,$line);

shift(@tokens);


my @genomes = ();
foreach $tok(@tokens){
  #print "$tok\n";
  if($tok =~ /unamb_read_count_(.*)/){
    #print "matched $tok $1\n"; 
    push(@genomes,$1);
  }
}
my $NGenomes = scalar(@genomes);
#print "$NGenomes @genomes\n";

my %hashC = {};
my $count = 0;
while($line = <FILE>){
  chomp($line);

  my @tokens = split(/\t/,$line);

  my $name = shift(@tokens);

  if($hashCluster{$name} ne undef){
   
    my $tcluster = $hashCluster{$name};
    my $MTotal = 0;
    my @unamb = ();

    for($i = 0; $i < $NGenomes; $i++){
      $unamb[$i] = $tokens[$i];
      $MTotal+=$unamb[$i];
    }

    #print "$name $tcluster $MTotal\n";

    my $label = "NA";
    my $max = 0;
    if($MTotal > $minCount){
      my $maxI = 0;
      for($i = 0; $i < $NGenomes; $i++){
   
    my $genus = $genomes[$i];
    #print "$i $genus $unamb[$i]\n";
    if($unamb[$i] > 0){
      my $f = $unamb[$i]/$MTotal;
      if($f > $max){
        $max  = $f;
        $maxI = $i;
      }
    }
      }

      if($max > $minZ){
    $label = $genomes[$maxI];
      }
    }
    print "$name,$label,$max,$MTotal\n";
  }
}

close(FILE);
