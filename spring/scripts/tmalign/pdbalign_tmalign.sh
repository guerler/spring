#!/usr/bin/perl

# includes
use File::Basename;
use strict;

# username
my $localpath = "/tmp/" . trim(`whoami`) . "/";

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 4)
{
	print "Parameters required: [target chain] [template list] [input path] [pdb repository]\n\n";
	print "Description: This script aligns one PDB file with all contained in the template directory and parses the results to the output directory. Required file format is ppppx.pdb.\n";	
	exit();
}

# read parameters
my $ida    = $ARGV[0];
my $inp    = $ARGV[1];
my $ipath  = $ARGV[2];
my $opath  = $ARGV[3] . "/" . substr($ida, 0, 2) . "/";

# make output directory
cmd ("mkdir -p $opath");

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
}
close (data);

# loop
my $count = 0;

# open output file
my $localout = "$localpath/$ida.tmalign";
open(OUTFILE, ">$localout");

# copy tmalign
cmd ("cp tmalign $localpath");

# align target to templates
my $idb = "";
foreach $idb (@list)
{
   	# log
   	if (++$count % 100 == 0)
   	{
        print "Aligning chains of pdb pair no. " . $count . " " . $ida . ", " . $idb . "\n";
    }

    # skip
    if ($ida eq $idb)
    {
        next;
    }

    # template chain path
    my $fa = $ipath . substr($ida, 0, 2) . "/" . $ida;
    my $fb = $ipath . substr($idb, 0, 2) . "/" . $idb;

    # check if file exists
    if (!(-e $fa && -e $fb))
    {
         print OUTFILE $idb . " -\n"
    }

   	# execute command
    open(IN, "/$localpath/tmalign $fa $fb 2>&1 |");
   	while (<IN>)
    {
        my $l = $_;

        if (length($l) < 50)
        {
            next;
        }
        
        if (substr($l, 0, 7) eq "Aligned")
        {
            my $tmscore = substr($l, 43, 7);
            if ($tmscore > 0.35)
            {
                my $out = $idb . " " . $tmscore . "\n";
                print OUTFILE $out;
            }
            last;
        }
    }
    close(IN);
}

# done   
print OUTFILE "DONE";
        
# close output
close(OUTFILE);

# move
cmd ("mv $localout $opath");

# copy and untar file
sub copy($)
{
    # get pars
    my $sour = shift;
    my $dest = shift;
    my $tmpfile = $dest . "/content.tar";

    # check if database is available
    my $tag = "$dest/tag.txt";
    unless (-e $tag)
    {
        cmd("mkdir -p $dest");
        cmd("echo . > $tag");
        cmd("cp $sour $tmpfile");
        cmd("tar -C $dest -xf $tmpfile");
        cmd("rm $tmpfile");
        cmd("echo .. > $tag");
    } else {
        # wait
        while (-s $tag == 2)
        {
            print ".";
            `sleep 20`;
        }
    }
}

# additional function
sub cmd($)
{
	my $string = shift;
	system($string);
}

# additional function
sub suffix($)
{
	my $f = shift;
	return (split(/\./,$f))[-1];
}

sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}

# trim
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
