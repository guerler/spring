#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# verify
if ($#ARGV + 1 < 2)
{
    print "Parameters required: [input directory] [output file]\n\n";
    print "Description: This script reads and concats files.\n";	
    exit();
}

# read parameters
my $inp = $ARGV[0];
my $out = $ARGV[1];

# open directory
opendir(D, $inp) || die "Can not open the directory : $!\n";

# loop
my @list = ();
my $n = 0;
while (my $f = readdir(D))
{
    if ($f ne "." && $f ne ".." && suffix($f) ne "pdb")
    {
        # data
        @list[$n] = trim($f);
        $n++;
    }
}

@list = sort (@list);

# open file
open (OUT, ">$out");

# open file
for (my $i = 0; $i < @list; $i++)
{
    # name
    my $name = $inp . $list[$i];

    # open file
    open (DATA, $name);
    my $key = -1;
    my $content = "";
    while (my $record = <DATA>)
    {
        if (substr($record, 0, 1) eq "#")
        {
            next;
        }

        if (substr($record, 0, 4) eq "DONE")
        {
            $key = 1;
        } else {
            if (trim($record) ne "")
	        {
                 $content .= $record;
            }
        }
    }
    close (DATA);

    # key
    if ($key == -1)
    {
        print "File $name is incomplete.\n";
        #system("mv $name ~/package/results/tmp/");
    }

    # content
    if ($content ne "")
    {
        print OUT $content;
    }
}
close (OUT);

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# additional function
sub suffix($)
{
	my $f = shift;
	return (split(/\./,$f))[-1];
}
