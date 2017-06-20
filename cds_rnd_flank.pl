#!/usr/bin/env perl

use strict;
use warnings;

my $cds_rxd_dir = "cds_rxd_dir";
my $cds_file = $ARGV[0] || "data/NC.edit.cds.fa";
my $red_file = $ARGV[1] || "data/sam30/rand.1.txt";
my $flank_len = $ARGV[2] || 30;

system("mkdir -p $cds_rxd_dir");

get_flank();

sub get_flank {
	my $cds = load_fasta($cds_file);
	my $red = load_red($red_file);
	open (OUT, ">$cds_rxd_dir/NC.rand.cds.flank" . 2*$flank_len . ".fa") or die $!;
	foreach my $seq_id (sort keys $red) {
		foreach my $pos (sort keys $red->{$seq_id}) {
			my $offset = $pos - $flank_len - 1;
			next if ($offset < 0);
			next if ($pos+$flank_len > length($cds->{$seq_id}));
			my $flank_seq = substr($cds->{$seq_id}, $offset, 2*$flank_len+1);
			print OUT ">", $seq_id, ".", $pos, "\n", $flank_seq, "\n";
		}
	}
	close OUT;
}

sub load_red {
	my $file = shift;
	my %hash = ();
	open (IN, $file) or die $!;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		my @w = split /\t/;
		$hash{$w[0]}{$w[1]} = 1;
	}
	close IN;
	return \%hash;
}

sub load_fasta {
	my $file = shift;
	my %seq = ();
	my $seq_id = undef;
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		if (/^\>(\S+)/) {
			$seq_id = $1;
		} else {
			s/\s+//g;
			$seq{$seq_id} .= $_;
		}
	}
	close IN;
	return \%seq;
}
