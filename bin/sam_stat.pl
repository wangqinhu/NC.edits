#!/usr/bin/env perl

use strict;
use warnings;

my $sam_dir = "dnds_vs_fnfs";
my %dnds_fnfs = ();

sam_stat($sam_dir);

sub sam_stat {
	my $sam_dir = shift;
	opendir (DIR, "$sam_dir") or die "Cannot open $sam_dir: $!\n";
	foreach my $file (sort readdir DIR) {
		next unless $file =~ /^sam.\d+/;
		load_dnds_fnfs("$sam_dir/$file");
	}
	closedir DIR;
	
	print_dnds();
	print_fnfs();
}

sub print_dnds {
	my @dnds = ();
	open (DNDS, ">$sam_dir/sam.dnds.txt") or die $!;
	for (my $i = 0; $i < 20; $i++) {
		foreach my $file (sort keys %dnds_fnfs) {
			push @dnds, $dnds_fnfs{$file}{$i}{'dnds'};
		}
		my $dnds_mean = mean(\@dnds);
		print DNDS $dnds_mean, "\n";
		@dnds = ();
	
	}
	close DNDS;
}

sub print_fnfs {
	open (FN, ">$sam_dir/sam.fn.txt") or die $!;
	open (FS, ">$sam_dir/sam.fs.txt") or die $!;
	foreach my $file (sort keys %dnds_fnfs) {
		for (my $i = 0; $i < 20; $i++) {
			print FN $dnds_fnfs{$file}{$i}{'fn'}, "\t";
			print FS $dnds_fnfs{$file}{$i}{'fs'}, "\t";
		}
		print FN "\n";
		print FS "\n";
	}
	close FN;
	close FS;
}

sub load_dnds_fnfs {
	my $file = shift;
	open (SAM, $file) or die $!;
	while (<SAM>) {
		chomp;
		next if /^\#/;
		next if /^\s*$/;
		my @w = split /\t/;
		$dnds_fnfs{$file}{$w[0]}{'dnds'} = $w[5];
		$dnds_fnfs{$file}{$w[0]}{'fn'} = $w[6];
		$dnds_fnfs{$file}{$w[0]}{'fs'} = $w[7];
	}
	close SAM;
}

sub mean {
	my $x = shift;
	my @num = @{$x};
	my $mean = 0;
	foreach (@num) {
		$mean += $_;
	}
	$mean /= @num;
	return $mean;
}
