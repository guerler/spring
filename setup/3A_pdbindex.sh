#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# strict
use strict;

# verify
if ($#ARGV + 1 < 3) {
	print "Parameters required: [list] [path to pdb repository] [database file]\n\n";
	print "Description: This script maps the extended monomer libraries complexation partners to the template database\n";	
	exit();
}

# read path to files
my $lst = $ARGV[0];
my $inp = $ARGV[1] . "splits";
my $sdb = $ARGV[1] . "database/fasta.db";
my $idx = $ARGV[2];

# temporary
my $record = "";

# open file
open (data, $lst);
my @lid = ();
my $n = 0;
while (my $record = <data>) {
    # data
    my @col = split(" ", $record);
    my $id = trim($col[0]);
    if ($id ne "") {
        @lid[$n] = $id;
        $n++;
    }
}
close (data);
               
# loop
open (idx, ">$idx");
for (my $i = 0; $i < $n; $i++) {
    # print
    print "Loading " . $lid[$i] . ".\n";

    # get id
    my $id = $lid[$i];

    # sub directory
    my $sub = $inp . "/" . substr($id, 0, 2) . "/" . $id;
    opendir(data, $sub);
    my @files = readdir(data);
    closedir(data);

    # sort files
    my @lis=sort @files;

    # print directory content
    print $sub . "\n";

    # split
    for (my $j = 0; $j < @lis; $j++) {
        # get sub name
        my $is = $lis[$j];

        # log
        print $is . "\n";

        # verify
        if ($is eq "." || $is eq "..") {
            print "Skipping.\n";
            next;
        }

        # get sequence
        my $seq = "";
        open (data, $sub . "/" . $is);
        while ($record = <data>) {
            # get key argument
            my $key = trim(substr($record, 0, 3));

            # read sequence
            if ($key eq "ATO") {
                my $kat = trim(substr($record, 12, 4));
                if ($kat eq "CA") {
                    my $kre = substr($record, 17, 3);
                    $seq .= getcharcode($kre);
                }
            }

            # end
            if($key eq "END") {
                last;
            }
        }
        close (data);

        # selected
        if ($seq ne "") {
            # make name
            my $tmpkey = "/tmp/" . $id . "." . $is;
            my $tmpseq = $tmpkey . ".sequence"; 
            my $tmpres = $tmpkey . ".result";

            # write file
            open (data, ">$tmpseq");	    
            print data ">$id\n$seq";
            close (data);

            # align
            cmd ("psiblast -query $tmpseq -db $sdb -out $tmpres");

            # open file
            open (data, "$tmpres");
            my $stage = 0;
            my $lng = 0;
            while ($record = <data>) {
                # get end
                if (length($record) > 1) {
                    # look for first length declaration
                    if ($stage == 0 && substr($record, 0, 7) eq "Length=") {
                        $lng = substr($record, 7);
                        if ($lng == 0) {
                            last;
                        }
                        $stage = 1; 
                    }
                    if ($stage == 1 && substr($record, 0, 1) eq ">") {
                        my $info = $id . " " . prefix($is) . " " . trim(substr($record, 1));
                        print idx "$info\n";
                        $stage = 2;
                    }
                    if ($stage == 2 && substr($record, 0, 13) eq " Identities =") {
                        my @values = split (" ", $record);
                        my @frac = trim(split ("/", $values[2]));                        
                        my $identity = int(($frac[0] / $lng) * 100.0) / 100.0;
                        last;
                    }
                }
            }
            close (data);

            # remove files
            cmd ("rm -f $tmpkey*");
        }
    }
}
print idx "DONE";
close (idx);

# residue codes
sub getcharcode ($) {
    my $s = shift;
	if ($s eq "ALA") { return "A"; }
	if ($s eq "CYS") { return "C"; }
	if ($s eq "ASP") { return "D"; }
	if ($s eq "GLU") { return "E"; }
	if ($s eq "PHE") { return "F"; }
	if ($s eq "GLY") { return "G"; }
	if ($s eq "HIS") { return "H"; }
	if ($s eq "ILE") { return "I"; }
	if ($s eq "LYS") { return "K"; }
	if ($s eq "LEU") { return "L"; }
	if ($s eq "MET") { return "M"; }
	if ($s eq "ASN") { return "N"; }
	if ($s eq "PRO") { return "P"; }
	if ($s eq "GLN") { return "Q"; }
	if ($s eq "ARG") { return "R"; }
	if ($s eq "SER") { return "S"; }
	if ($s eq "THR") { return "T"; }
	if ($s eq "VAL") { return "V"; }
	if ($s eq "TRP") { return "W"; }
	if ($s eq "TYR") { return "Y"; }
	return "X";
}

# additional function
sub cmd($) {
	my $string = shift;
	print ($string . "\n");
	system($string);
}

# additional function
sub prefix($) {
	my $f = shift;
	return (split(/\./,$f))[0];
}

# additional function
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
