#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script runs download PDB files and extract chains.\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $pathin = $ARGV[1] . "files/";
my $pathout = $ARGV[1] . "complex/";

# make
cmd ("mkdir -p $pathout");

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

# alphabet
my @alpha = ("A".."Z");

# report
open (report, ">log.txt");

# loop
for (my $i = 0; $i < $n; $i++) {
    # split
    my $id = substr($list[$i], 0 , 4);

    # verify
    my $pd = $pathin . lc($id) . ".pdb";
    my $fo = $pathout . $id . ".pdb";

    # print name
    print "$pd\n";  

    # check
    if (-e $fo) {
        print "Skipping.\n";
        next;
    }

    # check
    if (!(-e $pd)) {
        print "PDB file not found.\n";
        next;
    }

    # get all chains
    my $protein = "";  

    # remark 350
    my $biodef = "";
    my $biofound = 0;

    # open file
    open (data, "$pd");
    while (my $record = <data>) {
        # check 
        my $key = trim(substr($record, 0, 6));
        if ($key eq "ENDMDL" || $key eq "END") {
            last;
        }

        # collect all chains
        if ($key eq "ATOM") {
            $protein .= substr($record, 0, 55) . "\n";
        }

        # collect all chains
        if ($key eq "TER") {
            $protein .= "TER\n";
        }

        # check for remark
        if (trim(substr($record, 0, 23)) eq "REMARK 350 BIOMOLECULE:") {
            if ($biofound) {
                last;
            } else {
                $biofound = 1;
                $biodef .= "#BIOMOL 1\n";
            }
        }

        # check for remark
        if (trim(substr($record, 0, 41)) eq "REMARK 350 APPLY THE FOLLOWING TO CHAINS:") {
            $biodef .= "#CHAINS " . substr($record, 42);
        }

        # get rotation            
        if (trim(substr($record, 0, 18)) eq "REMARK 350   BIOMT") {
            $biodef .= "#ROTMAT " . substr($record, 22);
        }                
    }

    # end tag
    if ($biodef ne "") {
        $biodef .= "#END\n";
    }

    # close
    close (data);

    # check
    if ($biofound == 0) {
        print "No definiton found.\n";
        next;
    }

    # debug
    print $biodef;

    # final molecule
    my $biomol = "";

    # biodefinition
    my @lb = split ("\n", $biodef);

    # chains
    my @rotchains = ();
    my @rotmat = ();
    my $rotmatid = 0;

    # generate complex files
    for (my $j = 0; $j < @lb; $j++) {
        # key
        my $key = substr($lb[$j], 0, 7);

        # select chains
        if ($key eq "#CHAINS") {
            my $str = substr($lb[$j], 7);
            $str =~ s/,//g;
            @rotchains = split (" ", $str);
        }

        # select chains
        if ($key eq "#ROTMAT") {
            my @col = split (" ", $lb[$j]);
            $rotmat[$rotmatid][0] = $col[2];
            $rotmat[$rotmatid][1] = $col[3];                
            $rotmat[$rotmatid][2] = $col[4];                
            $rotmat[$rotmatid][3] = $col[5];                                                
            $rotmatid++;
        }

        # is complete?
        if ($rotmatid < 3) {
            next;
        }

        # rotate            
        for (my $c = 0; $c < @rotchains; $c++) {
            my $chain = $rotchains[$c];
            my @lines = split ("\n", $protein);
            for (my $k = 0; $k < @lines; $k++) {
                # lines
                my $line = $lines[$k];

                # is atom
                if (substr($line, 0, 4) eq "ATOM") {
                    # select chain
                    if (substr($line, 21, 1) ne $chain) {
                        next;
                    }

                    # get coordinates
                    my $x = trim(substr($line, 30, 8));
                    my $y = trim(substr($line, 38, 8));
                    my $z = trim(substr($line, 46, 8));

                    #debug
                    #print $x . " " . $y . " " . $z . "\n";

                    # rotate
                    my $newx = $x * $rotmat[0][0] + $y * $rotmat[0][1] + $z * $rotmat[0][2] + $rotmat[0][3];
                    my $newy = $x * $rotmat[1][0] + $y * $rotmat[1][1] + $z * $rotmat[1][2] + $rotmat[1][3];
                    my $newz = $x * $rotmat[2][0] + $y * $rotmat[2][1] + $z * $rotmat[2][2] + $rotmat[2][3];

                    #debug
                    #print $newx . " " . $newy . " " . $newz . "\n";

                    # write
                    $biomol .= substr($line, 0, 30);
                    $biomol .= sprintf ("%8.3f%8.3f%8.3f\n", $newx, $newy, $newz);
                }
            }

            # add termination
            $biomol .= "TER\n";

            # reset
            @rotmat = ();
            $rotmatid = 0;
        }

        # write
        open (fout, ">$fo");
        print fout $biomol;
        close (fout);
    }
}

# close report
close(report);

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

sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
