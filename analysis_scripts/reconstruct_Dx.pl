#! /usr/bin/perl -w
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

use Math::Trig;

sub get_elastic {

#
# this subroutine is intended to extract the elastic constas
# from a gulp output file
#

    my $filename = $_[0];
    my $compl_flag = 0;
    my @s;

    open GULPOUT, "<$filename" 
        or die "could not open $filename : $!\n";

    while (<GULPOUT>) {
        
# Get the Elastic Compliance Matrix
        if (/Elastic Compliance Matrix:/) {
            $compl_flag = 3;
        }
        
        if ((/-{79}/) and ($compl_flag > 0)) {
            $compl_flag --;
            next;
        }
        
        if ($compl_flag == 1) {
            my @line = split;
            shift @line;
            push (@s, @line); # @s is the matrix 11,12,13..16..21,22,23..65,66
        }
        
    }
    close GULPOUT;
    
#
# calculate components of reduced elactic complience matrix 
#

    my $S_44 = $s[21] - (($s[15]*$s[15])/$s[14]); # S44 = s44-(s34s34/s33) I hope - ref steeds pg 5
    my $S_55 = $s[28] - (($s[16]*$s[16])/$s[14]); # S55 = s55-(s35s35/s33)
    
    print "\nReduced elastic compliences: S44 = $S_44, S55 = $S_55\n";
    
    if ($S_44 == $S_55) {
print "This reduces to an isotropic solution ... so use isotropic insted!\n";
    }
    
    if (($S_44 == 0) or ($S_55 == 0)) {
        die "Reduced complience matrix has zero in 44 or 55!\n";
    }

    return ($S_44, $S_55);
}

sub anisotropic_1 {
# based on eq 3.13 from Steeds book. puts a screw dislocation into
# material along z (NB must then rotate xyz so it is allong x!)
    my ($x, $y, $z, $S_44, $S_55, $burg) = @_;

    # calculate angle around z from x.
    # (such that theta runs from 0 to 2pi
    # not -pi to pi twice!)

    my $angle = 0;
    if ($x == 0.0) {
        $x = 0.00000000000000000000000000000000000000000000001;
    }
    if (($y >= 0) and ($x > 0)) {
        $angle = atan (($y/$x) *
                       sqrt ($S_44/$S_55)) ;
    } elsif (($y >= 0) and ($x < 0)) {
        $angle = atan (($y/$x) *
                       sqrt ($S_44/$S_55)) + pi;
    } elsif (($y < 0) and ($x < 0)) {
        $angle = atan (($y/$x) * 
                       sqrt ($S_44/$S_55)) + pi;
    } elsif (($y < 0) and ($x > 0)) {
        $angle = atan (($y/$x) *
                       sqrt ($S_44/$S_55)) + 2*pi;
    } else {
        die "error in angle calc (sub anisotropic_1)!\n";
    }

    # work out displacments
    #IN THIS SCRIPT ONLY WE WORK BACKWARDS... SO - rahter than plus.
    $z -= $burg * ($angle/(2*pi)); 

    return $z;
}


sub angle {
# Calculates angle around x axis in radians
# result is between 0 and 2pi.
# call by &angle($x, $y, $z).

    my $z = pop(@_);
    my $y = pop(@_);
    my $x = pop(@_);
    my $angle = 3*pi;
    if ($z == 0.0) {
        $z = 0.00000000000000000000000000000000000000000000001;
    }
    if (($y >= 0) and ($z > 0)) {
        $angle = atan ($y/$z);
    } elsif (($y >= 0) and ($z < 0)) {
        $angle = atan ($y/$z) + pi;
    } elsif (($y < 0) and ($z < 0)) {
        $angle = atan ($y/$z) + pi;
    } elsif (($y < 0) and ($z > 0)) {
        $angle = atan ($y/$z) + 2*pi;
    } else {
        die "Unexpected error in angle calc (sub angle)!\n";
    }

    if (($angle > 2*pi) or ($angle < 0)) {
	die "error in angle calc (angle is $angle)\n";
    }   
  
    return $angle;
}

# Start of code - what a pile of crap this is...
my ($elast_file) = @ARGV;
my ($s44, $s55) = &get_elastic ($elast_file);

{ # A block to deal with input for a polymeric system
    print "Parsing a polymeric cell with two regions\n";
    
    my $start_flag = 0;
    my $end_flag = 0;
    my $cell_flag = 0;
    my $pot_flag = 0;
    while (<STDIN>) {
	
# Get atomic co-ordinates    
	if (/Polymer cell parameter/) {
	    $cell_flag ++;
	  ##  print "incremented cell flag...\n"; # debugging
	}

	if (($cell_flag == 1) and (/a =\s+(\d+\.\d+)/)) {
	    $a_cell = $1;
	    print "cell paramiter = $a_cell\n"; 
	}
# Get sataring atomic co-ordinates
	
# Does not work unlss co-ordinates include regions - need to fix...
	if (/Mixed fractional\/Cartesian coordinates of polymer :/) {
	    $start_flag = 6;
	}
	
	if ((/-{79}/) and ($start_flag > 0)) {
	    $start_flag --;
	    next;
	}
	
	if (($start_flag == 3) or ($start_flag == 1)) {
	    my @split_line = split;            #
	    push (@atom_number, $split_line[0]);
	    
	    $atom_name{$split_line[0]} = $split_line[1];
	    $atom_type{$split_line[0]} = $split_line[2]; # core or shel
	    $atom_start_x{$split_line[0]} = ($split_line[3] * $a_cell);
	    $atom_start_y{$split_line[0]} = $split_line[4];
	    $atom_start_z{$split_line[0]} = $split_line[5];
	    if ($atom_start_y{$split_line[0]} eq "*") {
		$atom_start_y{$split_line[0]} = $split_line[5];
		$atom_start_z{$split_line[0]} = $split_line[7];
	    }

	    next;
	}
	
# Get final atomic co-ordinates
	
	if (/Final fractional\/Cartesian coordinates of atoms :/) {
	    $end_flag = 3;
	}
	
	if ((/-{79}/) and ($end_flag > 0)) {
	    $end_flag --;
	    next;
	}
	
	if ($end_flag == 1) {
	    my @split_line = split;            #
	    
	    $atom_end_x{$split_line[0]} = ($split_line[3] * $a_cell);
	    $atom_end_y{$split_line[0]} = $split_line[4];
	    $atom_end_z{$split_line[0]} = $split_line[5];
	    next;
	}

    }   
}


# Calculate displacments..

#print "12345 12345 12 123456789abc 123456789abc 123456789abc 123456789abc 123456789abc 123456789abc\n";
my %bulk_x;
print "No    name   cs bulk cell x (a) start x (a)  start y (a)  start z (a)     end x (a)   end y (a)    end z (a) total dx (a)\n";
foreach (@atom_number) {

    $bulk_x{$_} = &anisotropic_1 ($atom_start_y{$_}, $atom_start_z{$_}, $atom_start_x{$_},
                                  $s44, $s55, $a_cell);
# Print basic atomic info
    printf "%-5d %-5s %-2s %12g %12g %12g %12g %12g %12g %12g %12g \n",
    $_, $atom_name{$_}, $atom_type{$_}, $bulk_x{$_}, $atom_start_x{$_}, $atom_start_y{$_},
    $atom_start_z{$_}, $atom_end_x{$_}, $atom_end_y{$_}, $atom_end_z{$_}, ($atom_end_x{$_}-$bulk_x{$_}); 
}

my $r = 3.5;
my $r2 = $r**2;
foreach my $OA (@atom_number) {
	if (($atom_name{$OA} =~ /O/) and ($atom_type{$OA} =~ /c/)) 
	{
		foreach $IA (@atom_number) {
			next if ($IA <= $OA);
			if (($atom_name{$IA} =~ /O/) and ($atom_type{$IA} =~ /c/)) 
			{
				my $dy2 = ($atom_end_y{$OA} - $atom_end_y{$IA})**2;
				if ($dy2 < $r2) 
				{
					my $dz2 = ($atom_end_z{$OA} - $atom_end_z{$IA})**2;
					if ($dz2 < $r2)
					{
						foreach my $xmult (-1 .. 1)
						{
							my $dx2 = (($atom_end_x{$OA}+($xmult*$a_cell)) - $atom_end_x{$IA})**2;
							if ($dx2 < $r2)
							{
								if (($dx2+$dy2+$dz2) < $r2)
								{
									my $mid_y = $atom_end_y{$OA} - 0.5*($atom_end_y{$OA}-$atom_end_y{$IA});
									my $mid_z = $atom_end_z{$OA} - 0.5*($atom_end_z{$OA}-$atom_end_z{$IA});
									my $length = ($atom_end_x{$OA}-$bulk_x{$OA}) - ($atom_end_x{$IA}-$bulk_x{$IA});
									my $angle = &angle(0, ($atom_end_y{$OA}-$atom_end_y{$IA}), 
												($atom_end_z{$OA}-$atom_end_z{$IA}));
									if ($length > (0.5*$a_cell)) 
									{
										$length -= $a_cell;
									}
									elsif ((-1*$length) > (0.5*$a_cell))
									{
										$length += $a_cell;
									}
									my $len_y = $length * cos($angle);
									my $len_z = $length * sin($angle);
									printf "ARROW %12f %12f %12f %12f %12f %12f \n", 
										$mid_y, $mid_z, $length, $angle, $len_y, $len_z;
								}
							}
						}
					}
				}
			}
		}
	}
}
	

