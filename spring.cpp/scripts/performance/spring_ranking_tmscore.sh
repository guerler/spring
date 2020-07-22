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
my $scale = 1000.0;

# histogram
my $h0 = ();
my $h1 = ();
my $h2 = ();
my $hn = ();

# count all
my $hall = 0;

# reset
my $m = int($scale) + 1;
for (my $i = 0; $i < $m; $i++)
{
    $h0[$i] = $h1[$i] = $h2[$i] = $hn[$i] = 0.0;
}

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>)
{    
    # data
    @values = split(" ", $record);

    # index
    my $index = $values[4];

    # check
    if ($index < 0)
    {
        $index = 0.0;
    }
    
    # value
    my $q = $values[5];

    # count
    $hn[$index]++;    
    $hall++;

    # backup
    if ($q > 0.6)
    {
        $h0[$index]++;
    } else {
       if ($q > 0.5)
       {
           $h1[$index]++;
       } else {
           $h2[$index]++;
       }
    }
}

# close
close (data);

# print out
for (my $i = 0; $i < $m; $i++)
{
    # distributions
    print $i . " ";
    print $h0[$i] . " ";
    print $h1[$i] . " ";
    print $h2[$i] . " ";
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
