#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# verify
if ($#ARGV + 1 < 1)
{
    print "Parameters required: [input fasta datbase file]\n\n";
    print "Description: This script copies protein information from fasta database files into an sql.\n";
    exit();
}

# read parameters
my $inp = $ARGV[0];
my $out = $inp . ".sql";

# check
if ($inp eq $out)
{
    print "Input and Output have to be different.\n";
    die();
}

# counter
my $n = 0;

# open file
open (INPUT, $inp);
open (OUTPUT, ">$out");
while (my $record = <INPUT>)
{    
    # check
    if (substr($record, 0, 1) ne ">")
    {
        next;
    }

    # current format
    #>b0002 aspartokinase I, homoserine dehydrogenase I (thrA) [2.7.2.4] {Escherichia coli K12-MG1655}

    # get entries
    my $code        = substr($record, 1, 5);
    my $info        = trim(substr($record, 6));

    # make sql
    print OUTPUT "INSERT INTO proteins VALUES('$code', '$info');\n";
    $n++;
}
close (INPUT);
close (OUTPUT);

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
