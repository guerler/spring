#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [file list] [target repository directory]\n\n";
	print "Description: This script generates HHM files results.\n";	
	exit();
}

# read path to files
my $list = $ARGV[0];
my $path = $ARGV[1];
my $local = trim(`pwd`);
print $local;

# make directory
cmd ("mkdir -p $path");

# load first directory
print "Loading file list: $list\n";
open DIR, $list or die $!;

# loop
my $id = 0;
while (<DIR>)
{
       my @values = split(" ", $_);
	my $f = trim($values[0]);
	unless (-e $path . "hhm/" . substr($f, 0, 2) . "/" . $f . ".hhm")
	{
		# print   
		print "Job $f submitted.\n";
		cmd("echo $f >> log.txt");
		cmd("./xsub.sh $local ./hhmake.sh $f $path");
		
		# increase counter
		$id++;

		# wait
		while ((`qstat | grep "spectrum_hhs" | wc -l`) >= 200)
		{
		    `sleep 10`;
		}
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
