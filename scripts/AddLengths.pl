#!/usr/bin/env perl


while($line = <STDIN>){
    chomp($line);

    my @tokens = split(/\t/,$line);

    print "$tokens[0]\t1\t$tokens[1]\n";
}
