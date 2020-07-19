#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 3)
{
	print "Parameters required: [file list] [target repository] [pdb repository]\n\n";
	print "Description: This script generates psi-blast output files.\n";
	exit();
}

# read path to files
my $flist = $ARGV[0];
my $patha = $ARGV[1] . "/fasta/";
my $pathb = $ARGV[2] . "/database/fasta.db";
my $pathc = $ARGV[1] . "/neighbor/blast/";

# load first directory
print "Loading file list: $flist\n";
open DIR, $flist or die $!;
my $local = trim(`pwd`);

# make
cmd ("mkdir -p $pathc");

# loop
my $id = 0;
while (<DIR>)
{
	# file
       my @values = split(" ", $_);
	my $fa = trim($values[0]);

       # path
       my $inppath = $patha . substr($fa, 0, 2) . "/";
       my $outpath = $pathc . substr($fa, 0, 2) . "/";
     
       # check
	unless (-e "$outpath$fa.out")
	{
		# wait
		while ((`qstat | grep "spectrum_blast" | wc -l`) >= 250)
		{
			`sleep 10`;
		}

              # make
              cmd ("mkdir -p $outpath");
        
		# print   
		print "Job $id submitted.\n";
		cmd("./psiblast_xsub.sh $local $inppath$fa $pathb $outpath$fa.out");

		# die
		#die();
	}

	# increase counter
	$id++;
}

# close file
close DIR;

# additional function
sub cmd($)
{
	my $string = shift;
	system($string);
}

# trim
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# additional function
sub suffix($)
{
	my $f = shift;
	return (split(/\./,$f))[-1];
}
