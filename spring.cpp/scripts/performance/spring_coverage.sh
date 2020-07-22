#!/usr/bin/perl

# use
use POSIX;

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# verify
if ($#ARGV + 1 < 1)
{
    print "Parameters required: [input file]\n\n";
    print "Description: This script generates derives a histogram.\n";	
    exit();
}

# read parameters
my $inp = $ARGV[0];
my $out = $ARGV[1];

# resolution
my $scale = 100.0;

# histogram
my $h0 = ();
my $h1 = ();
my $h2 = ();
my $h3 = ();

# reset
my $m = int($scale) + 1;
for (my $i = 0; $i < $m; $i++)
{
    $h0[$i] = $h1[$i] = $h2[$i] = $h3[$i] = 0.0;
}

# open file
open (data, $inp);
my @list = ();
my $n = 0;
my $hits = 0;
while (my $record = <data>)
{    
    # data
    @values = split(" ", $record);

    # tm-score
    my $index = int($values[5] * 10.0);
    $h0[$index]++;    

    # fnat
    my $index = int($values[6] * 10.0);
    $h1[$index]++;    

    # rmsd
    if ($values[9] > 0.50 && $values[9] != 2.15E+09)
    {
        my $index = int($values[8]);
        $h2[$index]++;
    }

    if ($values[9] > 0.90)
    {
        # hits below 2.5
        if ($values[8] < 2.5)
        {
            $hits++;
        }
    }

    # rmsd
    if ($values[11] > 0.50)
    {
        my $index = int($values[10]);
        $h3[$index]++;    
    }
        
    # count
    $hall++;
}

# close
close (data);

# number of hits
print "hits $hits\n";

# print out
for (my $i = 0; $i < $m; $i++)
{
    # distributions
    print $i . " ";
    print $h0[$i] . " ";
    print $h1[$i] . " ";
    print $h2[$i] . " ";
    print $h3[$i] . " ";

    # next
    print "\n";
}

sub max($)
{
	my $a = shift;
	my $b = shift;
	if ($a > $b)
	{
		return $a;
	}
	return $b;
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
