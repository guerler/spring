#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script copies the sequences in the fasta directory to the database file\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $path = $ARGV[1] . "fasta/";
my $outd = $ARGV[1] . "database";

# make
cmd ("mkdir -p $outd");

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
open (data, ">$outd/fasta.db");
for (my $i = 0; $i < $n; $i++) {
    # print 
    print "Loading " . $list[$i] . ".\n";

    # prepare name
    my $f = $path . substr($list[$i], 0, 2) . "/" . $list[$i];

	# file exists
	if (-e $f) {
		# open/read
		my $content = "";
		open (sequence, $f);
		while (my $record = <sequence>) {
			$content .= $record;
		}
		# print
		print data $content . "\n";
	} else {
		print "Skipping.\n";
	}
}
close (data);

# former formatdb
cmd ("makeblastdb -in $outd/fasta.db  -dbtype prot");

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
