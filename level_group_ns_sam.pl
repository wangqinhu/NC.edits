#!/usr/bin/env perl

use strict;
use warnings;

my $red_file = "data/NC.edit.cds.v20161024.txt";
my $cds_file = "data/NC.cds.seq.fas";

my $num_groups = $ARGV[0] || 5;
my $red_dir = $ARGV[1] || "level_group_ns/$num_groups/red/";
my $sam_dir = $ARGV[2] || "level_group_ns/$num_groups/sam/";

my $upstream = 2;
my $dnstream = 3;
my $num_sams = 10;
my $len_rand = $upstream + $dnstream + 1;

system("mkdir -p $red_dir");
system("mkdir -p $sam_dir");

sample_in_group();

sub sample_in_group {
	my $red = load_edit($red_file);
	my $cds = load_fasta($cds_file);
	my $grp = level2group($red, $num_groups);
	my $motif = fetch_motif($cds, $grp);
	motif2sam($cds, $motif, $num_sams);	
}

sub motif2sam {
	my $cds = shift;
	my $motif = shift;
	my $num_sams = shift;
	print "Sampling ...\n";
	foreach my $grp_id (sort by_num keys $motif) {
		print "Group $grp_id :";
		write_red($motif, $grp_id);
		my $frq = emit_from($motif, $grp_id);
		for (my $i = 1; $i <= $num_sams; $i++) {
			print " ", $i;
			dicer($cds, $motif, $grp_id, $i, $frq);
		}
		print "\n";
	}
	print "Done\n";
}

sub write_red {
	my $motif = shift;
	my $grp_id = shift;
	open (RED, ">$red_dir/red.$grp_id.txt") or die "Cannot open file $red_dir/red.$grp_id.txt: $!\n";
	foreach my $red_id (sort keys $motif->{$grp_id}) {
		print RED $red_id, "\n";
	}
	close RED;
}

sub dicer {
	my $cds = shift;
	my $motif = shift;
	my $grp_id = shift;
	my $num = shift;
	my $frq = shift;
	
	my $ids = fetch_ids($motif, $grp_id);

	open (SAM, ">$sam_dir/sam.$grp_id.$num.txt") or die "Cannot open file: $!\n";
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
	my $grp_id = shift;
	my %ids = ();
	foreach my $red_id (sort keys $motif->{$grp_id}) {
		my ($gene_id, $site) = split /\t/, $red_id;
		$ids{$gene_id}++;
	}
	my @ids = keys %ids;
	return \@ids;
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
	my $grp_id = shift;

	my $num_seq = 0;
	my %grp_motif = ();
	foreach my $red_id (sort keys $motif->{$grp_id}) {
		$grp_motif{$red_id} = $motif->{$grp_id}->{$red_id};
		$num_seq++;
	}

	# calculate prob.
	my %prob = ();
	my %sum = ();
	foreach my $seq (values %grp_motif) {
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
	my $cds = shift;
	my $grp = shift;
	my %motif = ();
	
	foreach my $grp_id (sort by_num keys %{$grp}) {
		foreach my $red_id (sort keys $grp->{$grp_id}) {
			my ($seq_id, $site) = split /\t/, $red_id;
			my $oligo = substr($cds->{$seq_id}, $site - $upstream - 1, $len_rand);
			$motif{$grp_id}{$red_id} = $oligo;
		}
	}
	return \%motif;
}

sub level2group {
	my $red = shift;
	my $num = shift;

	my %groups = ();
	my $grp_id = undef;
	my $red_id = undef;

	if ($num == 0) {
		die "Zero group specified!\n";
	} elsif ($num == 1) {
		foreach my $gene_id (sort keys $red) {
			foreach my $site (sort keys $red->{$gene_id}) {
				$grp_id = 1;
				$red_id = $gene_id . "\t" . $site;
				$groups{$grp_id}{$red_id} = $red->{$gene_id}->{$site}->{'level'};
			}
		}
		return \%groups;
	} else {
		foreach my $gene_id (sort keys $red) {
			foreach my $site (sort keys $red->{$gene_id}) {
				$grp_id = int($red->{$gene_id}->{$site}->{'level'} * $num - 0.0001) + 1;
				$red_id = $gene_id . "\t" . $site;
				$groups{$grp_id}{$red_id} = $red->{$gene_id}->{$site}->{'level'};
			}
		}
		return \%groups;
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
		$edit{$gene_id}{$edit_site}{'level'} = $w[6];
		$edit{$gene_id}{$edit_site}{'count'} = $w[4];
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

sub by_num {
	$a <=> $b;
}
