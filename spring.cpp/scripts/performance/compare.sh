#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# my
my $mn = 0.0;

# verify
if ($#ARGV + 1 < 2)
{
    print "Parameters required: [input file] [input file]\n\n";
    print "Description: This script generates compares two prediction files.\n";	
    exit();
}

# read parameters
my $inpa = $ARGV[0];
my $inpb = $ARGV[1];

# open file
open (data, $inpa);
my @namea = ();
my @refa = ();
my $na = 0;
while (my $record = <data>)
{    
    # data
    @values = split(" ", $record);
    $refa[$na]  = trim($record);
    $namea[$na] = $values[0] . "_" . $values[1];
    $na++;
}
close (data);

# open file
open (data, $inpb);
my @nameb = ();
my @refb = ();
my $nb = 0;
while (my $record = <data>)
{
    # data
    @values = split(" ", $record);
    $refb[$nb]  = trim($record);
    $nameb[$nb] = $values[0] . "_" . $values[1];
    $nb++;
}
close (data);

# open file
for (my $j = 0; $j < @namea; $j++)
{
	# search
	for (my $k = 0; $k < @nameb; $k++)        
	{
        if ($namea[$j] eq $nameb[$k])
	    {
            print $namea[$j] . " :1 " . $refa[$j] . " :2 " . $refb[$k] . "\n";
	        last;
	    }
	}
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
