#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# path
my $source_ftp = "https://files.rcsb.org/download/";
my $source_http = "https://files.rcsb.org/download/";

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script runs download PDB files and extract chains.\n";	
	exit();
}

# read path to files
my $inp   = $ARGV[0];
my $out  = $ARGV[1] . "files/";

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>) {    
    # data
    my @col = split(" ", $record);
    my $id = trim(substr($col[0], 0, 4));
    if ($id ne "") {
        @list[$n] = lc($id);
        $n++;
    }
}
close (data);

# make directory
cmd ("mkdir -p $out");
            
# loop
for (my $i = 0; $i < $n; $i++) {
   	# print 
	print "Loading " . $list[$i] . ".\n";

	# make names
   	my $id = $list[$i];
	my $pd = $id . ".pdb";
	my $zp = $pd . ".gz";

   	# add
	if (!(-e $out . $pd)) {
		# try to get it from ftp source
		cmd("wget $source_ftp$zp");
		cmd("mv " . $zp . ".* " . $zp);

		# verify
		if (-e $zp) {
			cmd("gunzip $zp");
		} else {
			cmd("wget -O $pd $source_http$pd");
		}
		if (-e $pd) {
			my $filesize = -s $pd;
			if ($filesize > 0) {
				cmd ("mv $pd $out");
			} else {
				cmd ("rm $pd");
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
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
