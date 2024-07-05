#!/usr/local/bin/perl
use strict;

my @genotype = ("homo","het");
my $g_flag = 0;
while(<>){
    my @x = split /\t/; 
    my @out;

    if($x[10] eq "SNP"){
        @out = ("s","control",$x[0], $x[1], $x[4], $x[9], $genotype[$g_flag]);
    } 
    elsif( $x[10] eq "INDEL" && length($x[4]) > length($x[9]) ){
        $x[1]++;
        my $len = length($x[4]) - 1;
        @out = ("d", "control", $x[0], $x[1], $len, $genotype[$g_flag]);
    } 
    elsif( $x[10] eq "INDEL" && length($x[4]) < length($x[9]) ){
        $x[1]++;
        my $insert = substr($x[9], 1, length($x[9]) - 1);
        @out = ("i", "control", $x[0], $x[1], $insert, $genotype[$g_flag]);
    }  

    print join("\t",@out)."\n";

    $g_flag = $g_flag == 0? 1 : 0;
}


