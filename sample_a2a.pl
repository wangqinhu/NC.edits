#!/usr/bin/perl

use strict;
use warnings;

my $red_file = $ARGV[0] || "data/NC.edit.cds.v20161024.txt";
my $cds_file = $ARGV[1] || "data/NC.cds.seq.fas";
my $red_dir = $ARGV[2] || "data/red";
my $sam_dir = $ARGV[3] || "data/sam";
my $a2a_dir = $ARGV[4] || "sample_a2a";
my $upstream = 2;
my $dnstream = 3;
my $num_sams = 100;
my $len_rand = $upstream + $dnstream + 1;
my $codon_table = codon_table();
my %a2a = a2a_ids();
my %a2a_color = id2score_color();
my $dig = 5;

run();

sub run {
	print "Start at ", scalar localtime, "\n";
	system("mkdir -p $red_dir");
	system("mkdir -p $sam_dir");
	system("mkdir -p $a2a_dir");
	my $cds = load_fasta($cds_file);
	my $red = load_edit($red_file);
	my $motif = fetch_motif($cds, $red);
	motif2sam($cds, $motif, $num_sams);
	cal_a2a($cds, $a2a_dir);
	print "End at ", scalar localtime, "\n";
}

sub cal_a2a {
	print "Generating report files ...\n";
	my $cds = shift;
	my $a2a_dir = shift;
	my %freq_red = ();
	my %freq_sam = ();
	
	my $a2a_prefix = $a2a_dir . '/a2a';

	# red
	my $red_tsv = "$red_dir/red.txt";
	my $edit = load_tsv($red_tsv);
	foreach my $seq_id (sort keys $edit) {
		foreach my $site (sort keys $edit->{$seq_id}) {
			my $codon = substr($cds->{$seq_id}, int(($site-1)/3)*3, 3);
			my $aa = codon2aa($codon);
			my @base = split //, $codon;
			$base[($site - 1) % 3] = 'I';
			my $codon_edit = $base[0] . $base[1] . $base[2];
			my $aa_edit = codon2aa($codon_edit);
			my $a2a = $aa . '2' . $aa_edit;
			$freq_red{$a2a}++;
		}
	}
	# sam
	for (my $i = 1; $i <= $num_sams; $i++) {
		my $sam_tsv = "$sam_dir/sam.$i.txt";
		my $edit = load_tsv($sam_tsv);
		foreach my $seq_id (sort keys $edit) {
			foreach my $site (sort keys $edit->{$seq_id}) {
				my $codon = substr($cds->{$seq_id}, int(($site-1)/3)*3, 3);
				my $aa = codon2aa($codon);
				my @base = split //, $codon;
				$base[($site - 1) % 3] = 'I';
				my $codon_edit = $base[0] . $base[1] . $base[2];
				my $aa_edit = codon2aa($codon_edit);
				my $a2a = $aa . '2' . $aa_edit;
				$freq_sam{$a2a}{$i}++;					
			}
		}
	}
	my %buffer = ();
	foreach my $a2a (sort keys %freq_red) {
			$freq_red{$a2a} = 0 if (!exists $freq_red{$a2a});
			$buffer{$a2a}{'red'}{'mean'} = $freq_red{$a2a};
			my @sam_a2a = ();
			for (my $i = 1; $i <= $num_sams; $i++) {
				$freq_sam{$a2a}{$i} = 0 if (!exists $freq_sam{$a2a}{$i});
				push @sam_a2a, $freq_sam{$a2a}{$i};
			}
			$buffer{$a2a}{'sam'}{'mean'} = mean(\@sam_a2a);
			$buffer{$a2a}{'sam'}{'sd'} = sd(\@sam_a2a);
	}
	# output/plot
	# red
	open (OUT, ">$a2a_prefix.red.txt") or die $!;
	print OUT "\tFreq\tColor\n";
	foreach my $a2a (sort by_a2a keys %freq_red) {
		print OUT $a2a;
		$buffer{$a2a}{'red'}{'mean'} = 0 if (!exists $buffer{$a2a}{'red'}{'mean'});
		print OUT "\t", $buffer{$a2a}{'red'}{'mean'}, "\t", $a2a_color{$a2a}, "\n";
	}
	close OUT;
	# sam_mean
	open (OUT, ">$a2a_prefix.sam.mean.txt") or die $!;
	print OUT "\tFreq\n";
	foreach my $a2a (sort by_a2a keys %freq_red) {
		print OUT $a2a;
		$buffer{$a2a}{'sam'}{'mean'} = 0 if (!exists $buffer{$a2a}{'sam'}{'mean'});
		print OUT "\t", $buffer{$a2a}{'sam'}{'mean'}, "\n";
	}
	close OUT;
	# sam_sd
	open (OUT, ">$a2a_prefix.sam.sd.txt") or die $!;
	print OUT "\tsd\n";
	foreach my $a2a (sort by_a2a keys %freq_red) {
		print OUT $a2a;
		$buffer{$a2a}{'sam'}{'sd'} = 0 if (!exists $buffer{$a2a}{'sam'}{'sd'});
		print OUT "\t", $buffer{$a2a}{'sam'}{'sd'}, "\n";
	}
	close OUT;
	# output/stat
	# red
	open (OUT, ">$a2a_prefix.red.stat.txt") or die $!;
	foreach my $a2a (sort by_a2a keys %freq_red) {
		print OUT $a2a, "\t";
	}
	print OUT "\n";
	foreach my $a2a (sort by_a2a keys %freq_red) {
		$freq_red{$a2a} = 0 if (!exists $freq_red{$a2a});
		print OUT $freq_red{$a2a}, "\t";
	}
	print OUT "\n";
	close OUT;
	# sam
	open (OUT, ">$a2a_prefix.sam.stat.txt") or die $!;
	foreach my $a2a (sort by_a2a keys %freq_red) {
		print OUT $a2a, "\t";
	}
	print OUT "\n";
	for (my $i = 1; $i <= $num_sams; $i++) {
		foreach my $a2a (sort by_a2a keys %freq_red) {
			$freq_sam{$a2a}{$i} = 0 if (!exists $freq_sam{$a2a});
			print OUT $freq_sam{$a2a}{$i}, "\t";
		}
		print OUT "\n";
	}
	close OUT;
}

sub by_a2a {
	$a2a{substr($a,0,3)} <=> $a2a{substr($b,0,3)};
}

sub a2a_ids {
	my %a2a = ();
	my @a2a = qw/K2E S2G N2D R2G Y2C I2V T2A M2V I2M K2R Q2R E2G H2R N2S D2G *2W *2* L2L V2V Q2Q P2P S2S E2E T2T A2A R2R K2K G2G/;
    my $i = 0;
	foreach my $id (@a2a) {
	    $a2a{$id} = $i;
	    $i++;
	}
	return %a2a;
}

sub id2score_color {
	my @a2a = qw/K2E S2G N2D R2G Y2C I2V T2A M2V I2M K2R Q2R E2G H2R N2S D2G *2W *2* L2L V2V Q2Q P2P S2S E2E T2T A2A R2R K2K G2G/;	
	# Accroding to doi: 10.1038/srep11550 Tabel 1
	my @a2a_score = qw/2 2 3 1 1 3 2 3 3 3 2 1 1 2 1 0 4 4 4 4 4 4 4 4 4 4 4 4/;
	my %hash = ();
	for (my $i = 0; $i < @a2a; $i++) {
		$hash{$a2a[$i]} = $a2a_score[$i];
	}
	return %hash;
}

sub motif2sam {
	print "Sampling, please wait a few minutes ...\n";
	my $cds = shift;
	my $motif = shift;
	my $num_sams = shift;
	write_red($motif);
	my $frq = emit_from($motif);
	for (my $i = 1; $i <= $num_sams; $i++) {
		dicer($cds, $motif, $i, $frq);
	}
}

sub write_red {
	my $motif = shift;
	open (RED, ">$red_dir/red.txt") or die "Cannot open file $red_dir/red.txt: $!\n";
	foreach my $red_id (sort keys $motif) {
		print RED $red_id, "\n";
	}
	close RED;
}

sub dicer {
	my $cds = shift;
	my $motif = shift;
	my $num = shift;
	my $frq = shift;
	my $ids = fetch_ids($motif);
	open (SAM, ">$sam_dir/sam.$num.txt") or die "Cannot open file: $!\n";
	foreach my $seq (reverse sort {$frq->{$a} <=> $frq->{$b}} keys $frq) {
		next if (length($seq) != $len_rand);
		my $nor = $frq->{$seq};
		my $sam = sampling($cds, $ids, $seq, $nor);
		foreach my $seq_id (keys $sam) {
			foreach my $pos (keys $sam->{$seq_id}) {
				print SAM $seq_id, "\t", $pos+$upstream+1, "\n";
			}
		}
	}
	close SAM;
}

sub fetch_ids {
	my $motif = shift;
	my %ids = ();
	foreach my $red_id (sort keys $motif) {
		my ($gene_id, $site) = split /\t/, $red_id;
		$ids{$gene_id}++;
	}
	return id_filter(\%ids);
}

sub id_filter {
	my $ids = shift;
	my %filtered_ids = ();
	foreach my $id (keys $ids) {
		if ($id =~ /(NCU\d{5}T)(\d)/) {
			$filtered_ids{$1}{$2}++;
		}
	}
	my @filtred_ids = ();
	foreach my $id (sort keys %filtered_ids) {
		my @num = sort keys $filtered_ids{$id};
		push @filtred_ids, $id . $num[0];
	}
	return \@filtred_ids;
}

sub sampling {
	my ($cds, $ids, $seq, $nor) = @_;
	my %uniq_sam = ();
	for (my $i = 0; $i < $nor; $i++) {
		my $offsets = -1;
		my $seq_id = undef;
		while ($offsets == -1) {
			$seq_id = $ids->[int(rand() * @{$ids})];
			$offsets = offsets($cds, $seq_id, $seq);
		}
		my $offset = $offsets->[int(rand() * @{$offsets})];
		if ($offset <= 0) {
			$i--;
			next;
		}
		if (exists $uniq_sam{$seq_id}{$offset}) {
			$i--;
			next;
		} else {
			$uniq_sam{$seq_id}{$offset} = 1;
		}
	}
	return \%uniq_sam;
}

sub offsets {
	my $cds = shift;
	my $seq_id = shift;
	my $motif = shift;
	my $seq = $cds->{$seq_id};
	my @offsets = ();
	my $offset = 0;
	while ($offset != -1) {
		$offset = index($seq, $motif, $offset+1);
		unless ($offset == -1){
			push @offsets, $offset;
		}
	}
	return -1 if @offsets < 1;
	return \@offsets;
}

sub emit_from {
	my $motif = shift;
	my $num_seq = 0;
	my %motif_seq = ();
	foreach my $red_id (sort keys $motif) {
		$motif_seq{$red_id} = $motif->{$red_id};
		$num_seq++;
	}
	# calculate prob.
	my %prob = ();
	my %sum = ();
	foreach my $seq (values %motif_seq) {
		next if (length($seq) != $len_rand);
		for (my $i = 0; $i < $len_rand; $i++) {
			my $char = substr($seq, $i, 1);
			$prob{$i}{$char}++;
			$sum{$i}++;
		}
	}
	# freq. to frac.
	for (my $i = 0; $i < $len_rand; $i++) {
		foreach my $char (keys $prob{$i}) {
			$prob{$i}{$char} = $prob{$i}{$char} / $sum{$i};
		}
	}
	# scale to 1, accumulated prob.
	for (my $i = 0; $i < $len_rand; $i++) {
		# A C G T
		if ($prob{$i}{'A'} == 1) {
			$prob{$i}{'C'} = 1;
			$prob{$i}{'G'} = 1;
			$prob{$i}{'T'} = 1;
		} else {
			$prob{$i}{'C'} = $prob{$i}{'A'} + $prob{$i}{'C'};
			$prob{$i}{'G'} = $prob{$i}{'C'} + $prob{$i}{'G'};
			$prob{$i}{'T'} = $prob{$i}{'G'} + $prob{$i}{'T'};
		}		
	}
	# emit
	my %seq_freq = ();
	my @base = qw/A C G T/;
	for (my $i = 0; $i < $num_seq; $i++) {
		my @seq = ();
		for (my $j = 0; $j < $len_rand; $j++) {
			my $base_prob = rand();			
			if ($base_prob <= $prob{$j}{'A'}) {
				$seq[$j] = 'A';
			} elsif ($base_prob <= $prob{$j}{'C'}) {
				$seq[$j] = 'C';
			} elsif ($base_prob <= $prob{$j}{'G'}) {
				$seq[$j] = 'G';
			} elsif ($base_prob <= $prob{$j}{'T'}) {
				$seq[$j] = 'T';
			} else {
				warn "Unknown prob. of base found!\n";
			}
		}
		if ($seq[$upstream] ne 'A') {
			$i--;
			next;
		}
		my $seq = join('', @seq);
		$seq_freq{$seq}++;
	}
	return \%seq_freq;
}

sub fetch_motif {
	print "Fetching edited motif ...\n";
	my $cds = shift;
	my $red = shift;
	my %motif = ();
	foreach my $seq_id (sort keys $red) {
		foreach my $site (sort keys $red->{$seq_id}) {
			my $oligo = substr($cds->{$seq_id}, $site - $upstream - 1, $len_rand);
			my $red_id = $seq_id . "\t" . $site;
			$motif{$red_id} = $oligo;
		}
	}
	return \%motif;
}

sub load_edit {
	print "Loading RNA editing file ...\n";
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
		$edit{$gene_id}{$edit_site}{'level'} = $w[6];
		$edit{$gene_id}{$edit_site}{'count'} = $w[4];
	}
	close IN;
	return \%edit;
}

sub load_fasta {
	print "Loading CDS file ...\n";
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

sub load_tsv {
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

sub by_num {
	$a <=> $b;
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

sub sd {
	my $x = shift;
	my @num = @{$x};
	my $mean = mean(\@num);
	my $sqt = 0;
	foreach (@num) {
		$sqt += ($mean - $_) ** 2;
	}
	my $std = ($sqt / (@num-1)) ** 0.5;
	return $std;
}
