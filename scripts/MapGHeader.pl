#!/usr/bin/env perl

my %mapOTU = {};
my %mapGenome = {};
my %mapID = {};

my $mapFile = $ARGV[0];

open(FILE,$mapFile) or die "Can't open $mapFile\n";

while($line = <FILE>){
  chomp($line);

  my @tokens = split(/\t/,$line);

  my $genome = $tokens[1];

  my $id = $tokens[0];

  $mapGenome{$genome} = $id;
}

close(FILE);

$line = <STDIN>;

my @tokens = split(/\t/,$line);
my @rtokens = ();

push(@rtokens,shift(@tokens));

foreach my $tok(@tokens){
    my $rtok = "";

    if($tok =~/unamb_read_count_(.*)/){
        $rtok = "unamb_read_count_$mapGenome{$1}";
    }
    elsif($tok=~/amb_read_count_(.*)/){
        $rtok = "amb_read_count_$mapGenome{$1}";
    }
    else{
        print "Error!\n";
    }
    push(@rtokens,$rtok);
}
my $rstring = join("\t",@rtokens);
print "$rstring\n";

while($line = <STDIN>){
  chomp($line);
  print "$line\n";  
}
