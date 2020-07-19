#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [fasta sequence database file] [output list file]\n\n";
	print "Description: This script retrieves all pdb identifiers from a resource file.\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $out = $ARGV[1];

# output file
open (OUT, ">$out");

# open file
open (data, $inp);
while (my $record = <data>)
{    
    # data
    my @col = split(" ", $record);
    my $id = trim(substr($col[0], 1));
    if (substr($record, 0, 1) eq ">")
    {
        # split name
        my @namevec = split ("\\|", $id);
        my $namevecsize = $#namevec + 1;

		# fix name
        if ($namevecsize > 1)
        {
			$id = $namevec[1];
		}        

		# print name
		print OUT $id . "\n";
    }
}
close (data);
close (OUT);

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

# additional function
sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}
# additional function
sub getname($)
{
	my $f = shift;
	return substr($f, rindex($f,"/") + 1, length($f) - rindex($f,"/") - 1);
}
# additional function
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
