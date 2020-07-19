#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);
strict;

# verify
if ($#ARGV + 1 < 1)
{
    print "Parameters required: [input SPRING results directory]\n\n";
    print "Description: This script converts a SPRING output into a valid SQL statement for our database.\n";
    exit();
}

# read parameters
my $indir = $ARGV[0];
my $out = "interactions.sql";

# open directory
opendir(D, $indir) || die "Can not open the directory : $!\n";

# loop
my @list = ();
my $n = 0;
while (my $f = readdir(D))
{
    if ($f ne "." && $f ne ".." && suffix($f) eq "results")
    {
        # data
        @list[$n] = trim($f);
        $n++;
    }
}

# found
print "Found $n result files.\n";

# loop
open (OUTPUT, ">$out");
for (my $i = 0; $i < $n; $i++)
{
    # open file
    my $inp = $list[$i];

    # split to get name
    @values = split("_", prefix($inp));
    my $ida = $values[0];
    my $idb = $values[1];

    # log
    # print $inp . " " . $ida . " " . $idb . "\n";

    # check
    if ($ida eq "" || $idb eq "")
    {
        next;
    }

    # check for file
    if (!(-e $indir . "/" . $ida . "_" . $idb . "_1.pdb"))
    {
        print "PDB file for $inp does not exist.";
        next;
    }

    # open file
    open (INPUT, $indir . "/" . $inp);
    while (my $record = <INPUT>)
    {
        # split columns
        @values = split(" ", $record);

        # check
        if ($values[0] ne "1")
        {
            next;
        }

        # get entries
        my $templatea   = $values[1];
        my $templateb   = $values[2];
        my $score       = $values[3];
        my $energy      = $values[4];
        my $tmscore     = $values[5];
        my $aligned     = $values[6];
        my $identity    = substr($values[7], 0, -1) / 100.0;

        # check
        if ($score < 5)
        {
            print "Score below expected threshold: $ida $idb $score.\n";
            next;
        }

        # make sql
        print OUTPUT "INSERT INTO interactions VALUES('$n', '$ida','$idb','$templatea','$templateb','$score','$energy','$tmscore','$aligned','$identity');\n";
        $n++;
    }
    close (INPUT);
}

close (OUTPUT);

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

# additional function
sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}
