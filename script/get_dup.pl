#!/usr/bin/perl -w
open IN,$ARGV[0];
open O1,">$ARGV[1]";
open O2,">$ARGV[2]";
while (<IN>){
	chomp;
	$bam = "$ARGV[5]/$ARGV[6].sort.bam";
	$out = "$ARGV[5]/$ARGV[6].marked_duplicates.$_.bam";
	$line = "$ARGV[3] $ARGV[4] $bam $out $_\n";
	print O1 "$line";
	print O2 "$out\n";
}
close IN;
close O1;
close O2;
