#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
    print "Description: This script separates multiple chains contained in a single PDB into multiple PDB files using the TER argument\n";	
	print "Additional : Find out how many complexes are available by using: find . -name '*_0_*.pdb' | wc -l\n";
	exit();
}

# maximum number of chains
my $MAXCHAINS = 50;

# read path to files
my $inp = $ARGV[0];
my $path = $ARGV[1] . "complex/";
my $outpath = $ARGV[1] . "splits/";

# make directory
cmd ("mkdir -p $outpath");

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
    # length of identifier
    my $lng = length($list[$i]);   

    # split
    my $id = substr($list[$i], 0, $lng - 1);
    my $ch = substr($list[$i], $lng - 1);

	# sub path
	my $pathsub = substr($id, 0, 2) . "/" . $id;

    # prepare name
    my $od = $outpath . $pathsub . $ch . "/";

    # check
    if (-e $od) {
        next;
    }

    # get files
    my @files = glob($path . $pathsub . "/*.*");
    for my $file (@files) {
        #
        # load file to array
        #
        my $nchains = 0;
        my @fc = ();
        open (data, $file);
        while ($record = <data>) {
            # get key argument
            my $key = substr($record, 0, 3);

            # get atom coordinates
            if ($key eq "ATO") {
                push (@fc, $record);                
            }

            # count chains
            if ($key eq "TER") {
                push (@fc, "TER");
                $nchains++;
                
                # set maximum number of chains
                if ($nchains >= $MAXCHAINS) {
					last;
				}
            }
        }

        # verify size
        my $fcn = @fc;
        if ($fcn == 0) {
            next;
        }

        # make sure file end with a TER
        if (substr($fc[$fcn - 1], 0, 3) ne "TER") {
            push (@fc, "TER");
            $nchains++;
        }
        close ($file);

        #
        # check for core chain
        #
        my $foundchain = 0;
        for (my $i = 0; $i < @fc; $i++) {
            # get key argument
            my $key = substr($fc[$i], 0, 3);

            # find chains
            if ($key eq "ATO") {
                my $cc = substr($fc[$i], 21, 1);
                if($ch eq $cc) {
                    $foundchain = 1;
                    last;
                }
            }
        }

        # verify that its at least a dimer
        if ($foundchain == 0 || $nchains < 2) {
            next;
        }

        # log
        print $file . "\n";

        # id
        my $biomol = prefix(getname($file));

        # make directory
        cmd ("mkdir -p $od");

        # separate chains
        my $nc = 0;
        my $cn = "";
        my $cc = "";

        # open pdb files
        for (my $i = 0; $i < @fc; $i++) {
            # get key argument
            my $key = substr($fc[$i], 0, 3);

            # open and close new files
            if ($key ne "TER" && $key ne "ATO") {
                next;
            }

            # verify
            if ($key eq "ATO" && ($cc eq "" || $cc eq substr($fc[$i], 21, 1))) {
                # get current chain
                $cc = substr($fc[$i], 21, 1);
                # write data
                $cn .= $fc[$i];                
            } else {
                # check
                my $on = "_1_";
                if ($cc eq $ch) {
                    $on = "_0_";
                }

                # write file
                open (out, ">" . $od . $biomol . $on . $nc . ".pdb");
                print out $cn . "TER";
                close (out);

                # chain counter
                $nc++;

                # reset content
                $cn = "";
                $cc = "";
            }
        }
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
