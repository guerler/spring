#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 3)
{
	print "Parameters required: [list of templates] [input path] [output path]\n\n";
	print "Description: This script aligns one PDB file with all contained in the template directory and parses the results to the output directory. Required file format is ppppx.pdb.\n";
	exit();
}

# read path to files
my $flist = $ARGV[0];
my $path  = $ARGV[1] . "/";
my $outp  = $ARGV[2] . "/";
my $local = trim(`pwd`);

# open file
open (data, $flist);
my @list = ();
my $n = 0;
while (my $record = <data>)
{
    # data
    my @values = split(" ", $record);

    # data
    @list[$n] = trim($values[0]);
    $n++;
}
close (data);

# loop
my $ida = "";
my $cx  = 0;
foreach $ida (@list)
{
	my $ofile = $outp . "/" . substr($ida, 0, 2) . "/" . $ida . ".tmalign";
	if (!(-e $ofile))
	{
		# wait
		while ((`qstat | grep "spectrum_tmalign" | wc -l`) >= 300)
		{
		    `sleep 60`;
		}

		# print
		print "Job $cx submitted.\n";
		cmd("./xsub.sh $local ./pdbalign_tmalign.sh $ida $flist $path $outp");
	}
	# increase counter
	$cx++;
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
