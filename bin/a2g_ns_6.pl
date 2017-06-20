#!/usr/bin/env perl

use strict;
use warnings;

# Cal N/S

my $cds_file = $ARGV[0] || "data/NC.edit.cds.fa";
my $red_flank = $ARGV[1] || "data/NC.edit.cds.flank6.fa";
my $len_xmer = 6;
my %codon_table = codon_table();

cal_nonsyn_syn();

sub cal_nonsyn_syn {
	my ($N, $S) = (0, 0);
	my $cds = load_fasta($cds_file);
	my $mer = generate_x_mer($cds, $len_xmer);
	my $red = load_red_oligo($red_flank);
	
	foreach my $oligo (sort keys %{$red}) {
		for (my $phase = 0; $phase < 3; $phase++) {
			if (exists $mer->{$oligo}->{$phase}) {
				my $codon = substr($oligo, $phase, 3);
				my $aa = codon2aa($codon);
				my @base = split //, $codon;
				$base[2-$phase] = 'I';
				my $codon_edit = $base[0] . $base[1] . $base[2];
				my $aa_edit = codon2aa($codon_edit);
				if ($aa ne $aa_edit) {
					$N += $mer->{$oligo}->{$phase};
				} else {
					$S += $mer->{$oligo}->{$phase};
				}
			}
		}
	}

	print "$N\t$S";
}

sub generate_x_mer {
	my ($cds, $len) = @_;
	my %xmer = ();
	my $oligo = undef;
	my $phase = undef;
	foreach my $seq_id (sort keys %{$cds}) {
		my $seq = $cds->{$seq_id};
		my $seq_len = length($seq);
		for (my $i = 0; $i < $seq_len-$len; $i++) {
			$oligo = substr($seq, $i, $len);
			$phase = $i % 3;
			$xmer{$oligo}{$phase}++;
		}
	}
	return \%xmer;
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

sub load_red_oligo {
	my $file = shift;
	my %oligo = ();
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		s/\s+//g;
		next if (/^\>/);
		$oligo{$_}++;
	}
	close IN;
	return \%oligo;
}

sub codon_table {
	my %table = (
	'TCA' => 'S',
	'TCC' => 'S',
	'TCG' => 'S',
	'TCT' => 'S',
	'TTC' => 'F',
	'TTT' => 'F',
	'TTA' => 'L',
	'TTG' => 'L',
	'TAC' => 'Y',
	'TAT' => 'Y',
	'TAA' => '*',
	'TAG' => '*',
	'TGC' => 'C',
	'TGT' => 'C',
	'TGA' => '*',
	'TGG' => 'W',
	'CTA' => 'L',
	'CTC' => 'L',
	'CTG' => 'L',
	'CTT' => 'L',
	'CCA' => 'P',
	'CCC' => 'P',
	'CCG' => 'P',
	'CCT' => 'P',
	'CAC' => 'H',
	'CAT' => 'H',
	'CAA' => 'Q',
	'CAG' => 'Q',
	'CGA' => 'R',
	'CGC' => 'R',
	'CGG' => 'R',
	'CGT' => 'R',
	'ATA' => 'I',
	'ATC' => 'I',
	'ATT' => 'I',
	'ATG' => 'M',
	'ACA' => 'T',
	'ACC' => 'T',
	'ACG' => 'T',
	'ACT' => 'T',
	'AAC' => 'N',
	'AAT' => 'N',
	'AAA' => 'K',
	'AAG' => 'K',
	'AGC' => 'S',
	'AGT' => 'S',
	'AGA' => 'R',
	'AGG' => 'R',
	'GTA' => 'V',
	'GTC' => 'V',
	'GTG' => 'V',
	'GTT' => 'V',
	'GCA' => 'A',
	'GCC' => 'A',
	'GCG' => 'A',
	'GCT' => 'A',
	'GAC' => 'D',
	'GAT' => 'D',
	'GAA' => 'E',
	'GAG' => 'E',
	'GGA' => 'G',
	'GGC' => 'G',
	'GGG' => 'G',
	'GGT' => 'G',
	);
	return %table;
}

sub codon2aa {
	my $codon = shift;
	$codon = uc($codon);
	$codon =~ s/I/G/;
	my $aa = $codon_table{$codon};
	return $aa;
}
