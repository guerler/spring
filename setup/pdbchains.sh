#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script separates multiple chains contained in a single PDB into multiple PDB files using the TER argument\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $path = $ARGV[1] . "files/";
my $outpath = $ARGV[1] . "chains/";

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>) {
    # data
    my @col = split(" ", $record);
    my $id = trim($col[0]);
    if ($id ne "") {
        @list[$n] = $id;
        $n++;
    }
}
close (data);

# loop
for (my $i = 0; $i < $n; $i++) {
    my $entry = $list[$i];
    print "Loading " . $entry . ".\n";

    # length of identifier
    my $lng = length($entry);

    # split
    my $id = lc(substr($entry, 0, $lng - 1));
    my $ch = substr($entry, $lng - 1);

    # prepare name
    my $od = $outpath . substr($entry, 0, 2) . "/";
    my $of = $od . $entry;

    # check
    if (-e $of) {
        print "Skipping.\n";
        next;
    }

	# make directory
	cmd ("mkdir -p $od");

    #
    # load chain to array
    #
    my $fc = "";
    my $targetfile = $path . $id . ".pdb";
    print "$targetfile\n";
    open (data, $targetfile);
    while ($record = <data>) {
        # get key argument
        my $key = substr($record, 0, 3);
        my $kch = substr($record, 21, 1);
        my $ktp = trim(substr($record, 12, 4));
        my $kfl = substr($record, 16, 1);
        
        # open and close new files
        if (($kfl eq " " || $kfl eq "A") && $key eq "ATO" && $ktp eq "CA" && ($kch eq $ch || ($kch eq " " && $ch eq "A"))) {
            $fc .= substr($record, 0, 55) . "\n";
        }
        
        # end
        if($key eq "END") {
            last;
        }
    }
    
    # write file
    if ($fc ne "") {
        open (out, ">>" . $of);
        print out ">" . $entry . "\n";
        print out $fc . "TER";
        close (out);
    }
}

# additional function
sub cmd($) {
	my $string = shift;
	print ($string . "\n");
	system($string);
}

# additional function
sub suffix($) {
	my $f = shift;
	return (split(/\./,$f))[-1];
}

# additional function
sub prefix($) {
	my $f = shift;
	return (split(/\./,$f))[0];
}

# additional function
sub getname($) {
	my $f = shift;
	return substr($f, rindex($f,"/") + 1, length($f) - rindex($f,"/") - 1);
}

# additional function
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
