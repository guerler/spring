#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script generates fasta sequence files from original pdb entries\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $path = $ARGV[1] . "files/";
my $outpath = $ARGV[1] . "fasta/";

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

#loop
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
	if (!(-e $od)) {
		cmd ("mkdir -p $od");
	}

    #
    # load chain to array
    #
    my $nc = 0;
    my $cx = 0;
    my $fc = "";
    my $targetfile = $path . $id . ".pdb";
    print "$targetfile\n";
    open (data, $targetfile);
    while ($record = <data>) {
        # get key argument
        my $key = substr($record, 0, 3);
        my $kch = substr($record, 21, 1);
        my $kat = trim(substr($record, 13, 3));

        # open and close new files
        if ($key eq "ATO" && $kch eq $ch && $kat eq "CA") {
        	# content
        	my $c = getCode(substr($record, 17, 3));

			# add amino acid        	
        	if ($c ne "") {
				# add
	            $fc .= $c;        
	            $nc++;

	            # add new line
    	        if ($cx != 0 && $cx % 50 == 0) {
					$fc .= "\n";
					$cx = 0;
				} else {
					$cx++;
				}
			}
        }

        # end
        if($key eq "END") {
            last;
        }
    }

    # write file if exceeds min length
    if ($nc >= 20) {
		$fc = ">" . $entry . "\n" . $fc;
	    open (out, ">" . $of);
    	print out $fc;
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
sub getCode ($) {
	my $seq = shift;
	my $c   = "";
	if($seq eq "GLY"){$c="G";}
	if($seq eq "ALA"){$c="A";}
	if($seq eq "VAL"){$c="V";}
	if($seq eq "LEU"){$c="L";}
	if($seq eq "ILE"){$c="I";}
	if($seq eq "MET"){$c="M";}
	if($seq eq "PHE"){$c="F";}
    if($seq eq "PRO"){$c="P";}
    if($seq eq "TYR"){$c="Y";}
    if($seq eq "TRP"){$c="W";}
    if($seq eq "LYS"){$c="K";}
    if($seq eq "SER"){$c="S";}
    if($seq eq "CYS"){$c="C";}
    if($seq eq "ASN"){$c="N";}
    if($seq eq "GLN"){$c="Q";}
    if($seq eq "HIS"){$c="H";}
    if($seq eq "THR"){$c="T";}
    if($seq eq "GLU"){$c="E";}
    if($seq eq "ASP"){$c="D";}
    if($seq eq "ARG"){$c="R";}
    if($c   eq ""   ){$c="X";}
    return $c;
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
