#!/usr/bin/env perl

use strict;
use warnings;

my $app = "bin/count_ns_from_red_list.pl";
my $level_group = "level_group_ns";
my $sam = 10;
my $dig = 2;

my %nsr = ();

# 1 group
my $red_1 = ns_ratio_red(1);
my $sam_1 = ns_ratio_sam(1);
# 5 group
my $red_5 = ns_ratio_red(5);
my $sam_5 = ns_ratio_sam(5);


$nsr{"red"}{"all"} = $red_1->[1];
for (my $i = 1; $i <= $sam; $i++) {
	$nsr{"sam" . $i}{"all"} = $sam_1->[1][$i];
}

for (my $i = 1; $i <= 5; $i++) {
	$nsr{"red"}{"itvl" . $i} = $red_5->[$i];
	for (my $j = 1; $j <= $sam; $j++) {
		$nsr{"sam" . $j}{"itvl" . $i} = $sam_5->[$i]->[$j];
	}	
}

open (STAT, ">$level_group/ns_ratio.txt") or die $!;
print STAT "total\titvl1\titvl2\titvl3\titvl4\titvl5\n";
foreach my $a (sort keys %nsr) {
	print STAT $a;
	foreach my $b (sort keys $nsr{$a}) {
		print STAT "\t", format_percentage($nsr{$a}{$b});
	}
	print STAT "\n";
}
close STAT;

sub ns_ratio_red {
	my $grp = shift;
	my @red = ();
	for (my $i = 1; $i <= $grp; $i++) {
		my $file = "$level_group/$grp/red/red.$i.txt";
		my $out = `perl $app $file`;
		$red[$i] = $out;
	}
	return \@red;
}

sub ns_ratio_sam {
	my $grp = shift;
	my $prefix = shift;
	my @sam = ();
	for (my $i = 1; $i <= $grp; $i++) {
		for (my $j = 1; $j <= $sam; $j++) {
			my $file = "$level_group/$grp/sam/sam.$i.$j.txt";
			my $out = `perl $app $file`;
			$sam[$i][$j] = $out;
		}
	}
	return \@sam;
}

sub format_percentage {
	my($num) = @_;
	my $percent = undef;
	if ( $num > 0 ) {
		$percent = 	int($num * 10**$dig + 0.5 ) / 10**$dig;
	} else {
		$percent = 	int($num * 10**$dig - 0.4 ) / 10**$dig;
	}
	return $percent;
}
