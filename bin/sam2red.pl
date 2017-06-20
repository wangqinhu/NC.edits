#!/usr/bin/env perl

use strict;

my $sam_dir = "data/sam";
my $red_dir = "data/sam.red";

system("mkdir -p $red_dir");

opendir (DIR, "$sam_dir") or die "Cannot open $sam_dir: $!\n";
foreach my $file (sort readdir DIR) {
	next unless $file =~ /.txt$/;
	sam2red("$sam_dir/$file", "$red_dir/$file");
}
closedir DIR;


sub sam2red {
	my $sam = shift;
	my $red = shift;
	open (SAM, $sam) or die $!;
	open (RED, ">$red") or die $!;
	while (<SAM>) {
		chomp;
		next if /^\#/;
		next if /^\s*$/;
		my @w = split /\t/;
		print RED "-\t-\t-\t-\t0\t0\t0\t" . $w[0] . ":c." . $w[1] . "A>G\n";
	}
	close SAM;
	close RED;
}
