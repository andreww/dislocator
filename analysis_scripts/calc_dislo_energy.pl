#!/usr/bin/perl 

#
#    Copyright 2010 Andrew Walker
#
#    This file is part of the Dislocator package.
#
#    The Dislocator package is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    The Dislocator package is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with the Dislocator package.  If not, see <http://www.gnu.org/licenses/>.
#

use warnings;
use strict;

my ($r, $bulk_file, $dislo_file) = @ARGV;

my $num = '[+-]?\d+\.\d+';
my $bulk_energy;
my $dislo_energy;


open (BULK, "< $bulk_file");
while (<BULK>) {
	if (/Polymer energy \(region 1\)  =\s+($num)\s+eV/) {
		$bulk_energy = $1;
	}
}
close BULK;

open (DISLO, "< $dislo_file");
while (<DISLO>) {
	if (/Polymer energy \(region 1\)  =\s+($num)\s+eV/) {
		$dislo_energy = $1;
	}
}
close DISLO;

print "$r $bulk_energy $dislo_energy " . ($dislo_energy - $bulk_energy) . "\n";
