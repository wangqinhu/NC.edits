#!/usr/bin/env perl

use warnings;
use strict;

# i/o
my $cds_rxd_dir = "cds_rxd_dir";
my $src_file = $ARGV[0] || "data/NC.edit.cds.flank6.fa";
my $cds_file = $ARGV[1] || "data/NC.edit.cds.fa";
my $num_sams = 1;
my $len_rand = 6;
my $red_site = 3;

# main
system("mkdir -p $cds_rxd_dir");
sam_edit();

# subroutine
sub sam_edit {
	system("mkdir -p $cds_rxd_dir/sam30");
	for (my $i = 0; $i < $num_sams; $i++) {
		dicer($i);
	}
}

sub dicer {
	my $num = shift;
	$num = $num + 1;
	open (SAM, ">$cds_rxd_dir/sam30/rand.$num.txt") or die "Cannot open file: $!\n";
	my $cds = load_fasta($cds_file);
	my $ids = read_ids($cds_file);
	my $frq = emit_from($src_file);
	foreach my $seq (reverse sort {$frq->{$a} <=> $frq->{$b}} keys $frq) {
		next if (length($seq) != $len_rand);
		my $nor = $frq->{$seq};
		my $sam = sampling($cds, $ids, $seq, $nor);
		foreach my $seq_id (keys $sam) {
			foreach my $pos (keys $sam->{$seq_id}) {
				print SAM $seq_id, "\t", $pos+$red_site, "\n";
			}
		}
	}
	close SAM;
}

sub emit_from {
	my $src_file = shift;
	my $num_seq = `grep -c '>' $src_file`;
	my %prob = ();
	my %sum = ();
	# calculate prob.
	open (IN, $src_file) or die "Cannot open $src_file: $!\n";
	while (<IN>) {
		chomp;
		unless (/^\>/) {
			my $seq = $1 if (/(\S+)/);
			next if (length($seq) != $len_rand);
			for (my $i = 0; $i < $len_rand; $i++) {
				my $char = substr($seq, $i, 1);
				$prob{$i}{$char}++;
				$sum{$i}++;
			}
		}
	}
	close IN;
	# freq. to frac.
	for (my $i = 0; $i < $len_rand; $i++) {
		foreach my $char (keys $prob{$i}) {
			$prob{$i}{$char} = $prob{$i}{$char} / $sum{$i};
		}
	}
	# scale to 1
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
		if ($seq[$red_site-1] ne 'A') {
			$i--;
			next;
		}
		my $seq = join('', @seq);
		$seq_freq{$seq}++;
	}
	
	foreach my $seq (keys %seq_freq) {
		print $seq, "\t", $seq_freq{$seq}, "\n";
	}

	return \%seq_freq;
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

sub read_ids {
	my $file = shift;
	my $text = `grep '>' $file | cut -f2 -d '\>'`;
	my @ids = split /\n/, $text;
	return \@ids;
}

sub read_frq {
	my $file = shift;
	my %freq = ();
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		unless (/^\>/) {
			$freq{$1}++ if (/(\S+)/);
		}
	}
	close IN;
	return \%freq;
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
