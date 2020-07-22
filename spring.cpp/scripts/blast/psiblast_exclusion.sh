#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# sequence identity cut off
my $cut   = 0.30;

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [list] [target repository]\n\n";
	print "Description: This script generates exclusion files.\n";
	exit();
}

# read path to files
my $names = $ARGV[0];
my $path  = $ARGV[1] . "/neighbor/blast/";
my $pout  = $ARGV[1] . "/neighbor/identity/";

# open file
open (data, $names);
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

# make output
cmd ("mkdir -p $pout");

# loop
my $id = 0;
for (my $i = 0; $i < @list; $i++)
{
    my $ida = prefix($list[$i]);
    my $tag = substr($ida, 0, 2) . "/";
 	my $inp = $path . $tag . $ida . ".out";
	my $out = $pout . $tag . $ida;

    # make output
    cmd ("mkdir -p $pout$tag");

	# verify
	if (!(-e $out))
	{
        # data
        my $name = "";
        my $lng = 0;

        # write result file
        open(OUTFILE, ">$out");
        open INFILE, $inp or die $inp;
        
        # loop
        my $stage = 0;
        while ($record = <INFILE>)
        {
            # new line
            if (trim ($record) eq "")
            {
                next;
            }
            
            # look for first length declaration
            if ($stage == 0 && substr($record, 0, 7) eq "Length=")
            {
                $lng = substr($record, 7);
                if ($lng == 0)
                {
                    print "Length invalid in $inp.\n";
                    die();
                }
                print $ida . " " . $lng;
                $stage = 1; 
            }
            
            # get target name
            if ($stage == 1 && substr($record, 0, 1) eq ">")
            {
                $name = trim(substr($record, 1));
                $stage = 2;
            }

            # add
            if ($stage == 2 && substr($record, 0, 13) eq " Identities =")
            {
                my @values = split (" ", $record);
                my @frac = trim(split ("/", $values[2]));
                
                # check
                my $identity = int(($frac[0] / $lng) * 100.0) / 100.0;
                if ($identity > $cut)
                {
                    print OUTFILE $name . " " . $identity . "\n";
                }
                $stage = 1;                
            }
            
            # finalize
            if (substr($record, 0, 6) eq "Window")
            {
                $stage = 3;
            }
        }
        close(INFILE);
        close(OUTFILE);
        
        # check
        if ($stage != 3)
        {
            print "File $inp incomplete.\n";
            `mv $inp ~/temp/`;
        }
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

sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}

# additional function
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
