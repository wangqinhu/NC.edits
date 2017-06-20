#!/usr/bin/env perl

use strict;
use warnings;

my $edit_file = "data/NC.edit.cds.v20161024.txt";
my $num_of_bin = $ARGV[0] || 20;
my $cds_seq_file = $ARGV[1] || "data/NC.cds.seq.fas";

my $expr_file = "data/fpkm.txt";
my $expr_bin = "expr_bin";
my $out_dir = "expr_vs_fnfs";

my $cds_seq = load_fasta($cds_seq_file);
my $codon_table = codon_table();

system("mkdir -p $out_dir");
expr_red();

sub expr_red {

	# load data
	my $expr = load_expr($expr_file);
	my $edit = load_edit($edit_file);

	# expr group ids
	my @grp_id = ();
	my $num_of_genes = keys $expr;
	my $sub_num = int($num_of_genes / $num_of_bin + 0.5) ;
	my ($i, $j, $k) = (0, 0, 0);
	for my $gene_id (sort {$expr->{$a} cmp $expr->{$b}} keys $expr) {
		$i = int ($k / $sub_num);
		next if $i >= $num_of_bin;
		$j = $k % $sub_num;
		$grp_id[$i] .= $gene_id . "\n";
		$k++;
	}

	# edited var
	my @genes = ();
	my @N = ();
	my @S = ();
	my @n = ();
	my @s = ();
	my @expr = ();

	open (TAB, ">$out_dir/expr_fnfs.txt") or die $!;
	# print header
	print TAB "grp\t\texpr\tN\tS\tn\ts\tfn\tfs\n";

	system("mkdir -p $expr_bin");
	
	for my $i (0..$#grp_id) {
		# calculate n and s
		my $id_all = $grp_id[$i];
		my @sub_ids = split /\n/, $id_all;
		foreach my $gene_id (@sub_ids) {
			next unless (exists $edit->{$gene_id});
			foreach my $site (sort keys $edit->{$gene_id}) {
				if (is_nonsyn($gene_id, $site)) {
					$n[$i]++;
				} else {
					$s[$i]++;
				}
			}
			$genes[$i] .= $gene_id . "\n";
			$expr[$i] .= $expr->{$gene_id} . "\n";
		}

		# calculate N and S
		open (OUT, ">$expr_bin/grp.$i.id.txt") or die $!;
		print OUT $genes[$i];
		close OUT;
		system("perl ./bin/id2seq.pl $expr_bin/grp.$i.id.txt $cds_seq_file $expr_bin/grp.$i.cds.fa");
		my $ns = `perl ./bin/a2g_ns_6.pl $expr_bin/grp.$i.cds.fa`;
		($N[$i], $S[$i]) = split /\t/, $ns;

		# calculate fn and fs
		my $fn = $n[$i]/$N[$i];
		my $fs = $s[$i]/$S[$i];
		
		# calculate mean fnfs
		my @expr_vals = split /\n/, $expr[$i];
		my $median_expr = median(@expr_vals);

		my $line = join("\t", $i,  $median_expr, $N[$i], $S[$i], $n[$i], $s[$i], $fn*1000, $fs*1000);
		$line =~ s/\t$//;
		print TAB $line, "\n";
	}

	system("rm -rf $expr_bin");
	close TAB;

}

sub median {
	my @vals = sort {$a <=> $b} @_;
	my $len = @vals;
	if($len % 2) {
		return $vals[int($len/2)];
	} else {
		return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
	}
}

sub load_expr {
	my $file = shift;
	my %hash = ();
	open (IN, $file) or die $!;
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		next if /\tNA/;
		my @w = split /\t/;
		$hash{$w[0]} = $w[1];
	}
	close IN;
	return \%hash;
}

sub is_nonsyn {
	my $gene_id = shift;
	my $site = shift;
	my $codon = substr($cds_seq->{$gene_id}, int(($site-1)/3)*3, 3);
	my $a_site = ($site-1) % 3;
	my @codon = split //, $codon;
	$codon[$a_site] = 'G';
	my $codon_edit = join('', @codon);
	if (codon2aa($codon) eq codon2aa($codon_edit)) {
		return 0;
	} else {
		return 1;
	}
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

sub codon2aa {
	my $codon = shift;
	$codon = uc($codon);
	$codon =~ s/I/G/;
	my $aa = $codon_table->{$codon};
	return $aa;
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
	return \%table;
}
