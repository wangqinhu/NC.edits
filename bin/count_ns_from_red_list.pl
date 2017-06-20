#!/usr/bin/env perl

my $red_file = $ARGV[0] || "sam.tsv";
my $cds_file = "data/NC.cds.seq.fas";
my $codon_table = codon_table();
my $cds_seq = load_fasta($cds_file);

my $red = load_red_tsv($red_file);

nonsyn_syn($red);

sub nonsyn_syn {
	my $red = shift;
	my ($n, $s) = (0, 0);
	foreach my $gene_id (sort keys $red) {
		foreach my $site (sort keys $red->{$gene_id}) {
			if (is_nonsyn($gene_id, $site)) {
				$n++;
			} else {
				$s++;
			}
		}
	}
	print $n/$s;
}

sub load_red_tsv {
	my $file = shift;
	my %hash = ();
	open (IN, $file) or die $!;
	while (<IN>) {
		chomp;
		next if /^\#/;
		next if /^\s*$/;
		my @w = split /\t/;
		$hash{$w[0]}{$w[1]} = 1;
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
		$edit{$gene_id}{$edit_site}{"count"} = $w[4];
	}
	close IN;
	return \%edit;
}
