#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 4)
{
	print "Parameters required: [target list] [target repository] [pdb repository] [output path]\n\n";
	print "Description: This script runs SPRING on all given protein pairs.\n";
	exit();
}

# read path to files
my $inp   = $ARGV[0];
my $path  = $ARGV[1] . "/";
my $tmpl  = $ARGV[2] . "/";
my $pout  = $ARGV[3] . "/";

# size
my $nsize = 2;
my $local = trim(`pwd`);

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>)
{
	# data
	my @values = split(" ", $record);
        
	# data
	@list[$n] = trim($values[0]);
	$n++;

	# add second column
	if (@values > 1)
	{
		# data
		@list[$n] = trim($values[1]);
		$n++;
	}
}
close (data);

# names
my $feature = $pout . "/features/";
my $batch = $pout . "/batch/";

# make
cmd("mkdir -p $feature");
cmd("mkdir -p $batch");

# loop
for (my $i = 0; $i < $n; $i += $nsize)
{
    # update batch name if there is only one entry per file
    my $batchname = $i;
    if ($nsize == 2)
    {
        $batchname = $list[$i] . "_" . $list[$i+1];
    }

    # split
    my $bf = $batch . $batchname . ".txt";
            
    # check            
    if (!(-e $bf))
    {
        # generate file
        open (bfile, ">$bf");
        for (my $j = 0; $j < $nsize; $j++)
        {
            # done
            if ($j + $i >= $n)
            {
                last;
            }
            
            # add
            print bfile $list[$j + $i] . "\n";
        }
        close (bfile);

        # wait
        while ((`qstat | grep "spring_rx" | wc -l`) >= 100)
        {
            `sleep 20`;
        }

        # print
        print "Job $i submitted.\n";
        cmd("./xsub.sh $local ./spring_tasser_refinement $bf $path $tmpl $feature");
        #die();
    }
}

# close file
close DIR;

# additional function
sub cmd($)
{
	my $string = shift;
	print ($string . "\n");
	system($string);
}

# additional function
sub suffix($)
{
	my $f = shift;
	return (split(/\./,$f))[-1];
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
