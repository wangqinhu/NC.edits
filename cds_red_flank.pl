#!/usr/bin/env perl

use strict;
use warnings;

my $cds_rxd_dir = "cds_rxd_dir";
my $edit_file = "data/NC.edit.cds.v20161024.txt";
my $cds_seq_file = "data/NC.cds.seq.fas";

my $up_position = $ARGV[0] || 2;
my $motif_len = $ARGV[1] || 6;

system("mkdir -p $cds_rxd_dir");

my $edit = load_edit($edit_file);
my $seq = load_fasta("$cds_seq_file");

open (OUT, ">$cds_rxd_dir/NC.edit.cds.flank$motif_len.fa") or die $!;
foreach my $seq_id (keys $edit) {
	foreach my $site (keys $edit->{$seq_id}) {
		my $str = substr($seq->{$seq_id}, $site - $up_position - 1, $motif_len);
		print OUT ">$seq_id\.$site\n";
		print OUT $str, "\n";
	}
}
close OUT;

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

sub load_edit {
	my $file = shift;
	my %edit = ();
	open (IN, $file) or die "Cannot open file $file: $!\n";
	while (<IN>) {
		chomp;
		next if /^\#/;
		next if /^Chromosome/;
		next if /^\s*$/;
		my @w = split /\t/;
		my $gene_info = $w[7];
		my $gene_id = undef;
		my $edit_site = undef;
		$gene_info =~ s/\[//g;
		$gene_info =~ s/\]//g;
		if ($gene_info =~ /(\S+)\:\S\.(\d+)A\>G/) {
			$gene_id = $1;
			$edit_site = $2;
		} else {
			warn "invlaid edit site found in $_\n";
		}
		next unless (defined $gene_id);
		$w[4] = 0 if $w[4] =~ /^\D+$/;
		$w[6] = 0 if $w[6] =~ /^\D+$/;
		$edit{$gene_id}{$edit_site}{"level"} = $w[6];
		$edit{$gene_id}{$edit_site}{"coverage"} = $w[5];
		$edit{$gene_id}{$edit_site}{"count"} = $w[4];
	}
	close IN;
	return \%edit;
}
