#!/usr/bin/env perl

my $token = $ARGV[0];
my $min = 0.5;

while($line = <STDIN>){
    chomp($line);

    my @tokens = split(/,/,$line);

    $tok = $tokens[$token];

    if($tok =~/(.*)->(.*)->(.*)/){
        if ($2 > $min && $3 > $min){
            print "$tokens[0],$1\n";
        }
    }
}
