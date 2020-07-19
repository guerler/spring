#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 3)
{
	print "Parameters required: [list] [target repository] [database]\n\n";
	print "Description: This script performs HHSearch.\n";	
	exit();
}

# read path to files
my $flist = $ARGV[0];
my $path  = $ARGV[1];
my $db    = $ARGV[2];

# get local
my $local = trim(`pwd`);

# load first directory
print "Loading file list: $flist\n";
open DIR, $flist or die $!;

# loop
my $id = 0;
while (<DIR>)
{
	my $f = trim($_);
	my $fout = $path . "/hhr/" . substr($f, 0, 2) . "/" . $f . ".hhr";
	
	# check
	unless (-e $fout)
	{
		# wait
		while ((`qstat -uguerler | grep "spectrum_hhs" | wc -l`) >= 300)
		{
			`sleep 20`;
		}

		# print   
		print "Job $f submitted.\n";
		cmd("./xsub.sh $local ./hhsearch.sh $f $path $db");

		# increase counter
		$id++;
		#die();
	}
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
