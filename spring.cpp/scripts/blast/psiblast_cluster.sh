#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# sequence identity cut off
my $cut   = 0.90;

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [list] [target repository]\n\n";
	print "Description: This script generates exclusion files.\n";
	exit();
}

# read path to files
my $names = $ARGV[0];
my $path  = $ARGV[1] . "/exclusion/identity/";
my $fout  = $ARGV[1] . "/exclusion/cluster";

# open file
open (data, $names);
my @list = ();
my %mark = ();
my $n = 0;
while (my $record = <data>)
{
    # data
    my @values = split(" ", $record);

    # data
    my $id = trim($values[0]);
    $list[$n] = $id;
    $mark{$id} = 1;
    $n++;
}
close (data);

# open
open(OUTFILE, ">$fout");

# loop
my $id = 0;
for (my $i = 0; $i < @list; $i++)
{
    my $id  = $list[$i];
    my $tag = substr($id, 0, 2) . "/";
 	my $inp = $path . $tag . $id;

	# verify
	if ($mark{$id} == 1)
	{
        # print main id
        print OUTFILE $id . " ";
        $mark{$id} = 0;
                    
        # check
        unless (-e $inp)
        {
            next;
        }
                    
        # write result file
        open INFILE, $inp or die $inp;

        # loop
        while ($record = <INFILE>)
        {
            # values
            my @values = split (" ", $record);

            # check size
            if (scalar(@values) == 2)
            {
                # items
                my $ida = $values[0];
                my $sim = $values[1];
                
                # check
                if ($mark{$ida} == 1 && $sim > $cut)
                {
                    print OUTFILE $ida . " ";
                    $mark{$ida} = 0;
                }
            }
        }
        
        # next line
        print OUTFILE "\n";
                    
        # close input
        close(INFILE);
    }
}

# close
close(OUTFILE);

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
