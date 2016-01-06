#!/usr/bin/perl -w

use strict;
use autodie;

my $win = 200;
my $nshifts = 200;
my $nfiles = 4;

my $outprefix = "simulated.shape.data.";
my $outfile;

my @ref = ();

for (my $i = 0; $i < $win; $i++){
	my $n = int(rand(1000));
	push @ref, $n;
}

for (my $f = 0; $f < $nfiles; $f++){
	$outfile = $outprefix . $f . ".txt"; 
	open OUTFILE, '>', $outfile;
	print OUTFILE "ref\tNA\tNA";
	for (my $i = 0; $i < $win; $i++){
		print OUTFILE "\t";
		print OUTFILE $ref[$i];
	}
	print OUTFILE "\tNA\tNA\n";
	close OUTFILE;
}



for (my $n = 0; $n < $nshifts; $n++){

	my $neg = int(rand(2));
	my $shift = int(rand(26));
	my $flip = int(rand(2));
	$shift = -1*$shift if ($neg);

	my $site_name = $n . "_" . $shift . "_" . $flip;
	
	@ref = reverse(@ref) if ($flip);
	
	for (my $f = 0; $f < $nfiles; $f++){
		$outfile = $outprefix . $f . ".txt"; 
		
		open OUTFILE, '>>', $outfile;
		
	
		print OUTFILE $site_name, "\tNA\tNA\t";
	
		for (my $i = 0; $i < $win; $i++){
			print OUTFILE "\t";
			if ($i + $shift >= 0 && $i + $shift < $win){
				print OUTFILE $ref[$i+$shift];
			} else {
				print OUTFILE int(rand(1000));
			}
		}
		print OUTFILE "\tNA\tNA\n";
		
		close OUTFILE;
	}
	
	# flip back if flipped
	@ref = reverse(@ref) if ($flip);
}

exit;