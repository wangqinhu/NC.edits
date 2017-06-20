#!/usr/bin/env perl

use strict;
use warnings;

extract_data();

sub extract_data {
	my $file = $ARGV[0];
	open (IN, $file) or die "Cannot open $file: $!\n";
	my %fgvs = ();
	my %type = ();
	my $id = undef;
	my $site = undef;
	my $buffer = '';
	my %print_once = ();
	while (<IN>) {
		next if /^\s*$/;
		if (/^#\s+\S+\s+(\S+)\s+(\S+)\s+/) {
			if (!exists $print_once{$buffer}) {
				print $buffer;
				$print_once{$buffer} = 1;
			}
			$id = $1;
			$site = $2;
			if (/syn/) {
				$type{$id}{$site} = "syn";
			} else {
				$type{$id}{$site} = "non";
			}
		} elsif (/^N/ or /^Sm/) {
			my @w = split /\s+/, $_, 3;
			my $str = '';
			if ($w[1] =~ /([A|G|C|T])/) {
				$str = $1;
			}
			$fgvs{$id}{$site} .= $str;
		} else {
			if (exists $fgvs{$id}{$site}) {
				if (length($fgvs{$id}{$site}) == 3 && $fgvs{$id}{$site} =~ /A\SA/) {
					$buffer = "$id\t$site\t$fgvs{$id}{$site}\t$type{$id}{$site}\n";
				}
			}
		}
	}
	if (!exists $print_once{$buffer}) {
		print $buffer;
	}
}
