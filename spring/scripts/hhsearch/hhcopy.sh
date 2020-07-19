#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# verify
if ($#ARGV + 1 < 3)
{
    print "Parameters required: [dimer list] [input directory] [output directory]\n\n";
    print "Description: This script renames the proteins in hhm files.\n";	
    exit();
}

# read parameters
my $inp = $ARGV[0];
my $path = $ARGV[1];
my $pout = $ARGV[2];

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>)
{    
    # data
    @list[$n] = trim($record);
    $n++;
}
close (data);

# open file
for (my $i = 0; $i < @list; $i ++)
{
    unless (-e $pout . $list[$i] . ".hhm")
    {
        if (-e $path . $list[$i] . ".hhm")
        {
            cmd("cp " . $path . $list[$i] . ".hhm " . $pout . $list[$i] . ".hhm");
            #die();
        }
    }
}

# additional function
sub cmd($)
{
	my $string = shift;
	print $string . "\n";
	system($string);
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
