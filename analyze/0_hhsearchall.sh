#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;
use File::Basename;

# verify
if ($#ARGV + 1 < 3) {
	print "Parameters required: [directory with files '.fasta'] [HH-search binary path] [HH-search database path]\n\n";
	print "Description: This script performs HH-search on all matched files in the directory.\n";
	exit();
}

# read path to files
my $path   = $ARGV[0];
my $binary = $ARGV[1];
my $db     = $ARGV[2];

# get files
my @files = glob($path . "*.fasta");
for my $file(@files) {
	my $name = fileparse($file, '\..*');
	my $output = $path . $name . ".hhr";
	print("\nProcessing: " . $name . "\n");
	if (-e $output) {
        print("Already available.\n");
        next;
    }
	my $ready = false;
	my $command = "$binary -i $file -d $db -ohhm $path$name.hhm -oa3m $path$name.a3m &";
	print ("Executing: $command\n");
	print "Press ENTER to continue or CTRL-C to exit.";
	<STDIN>;
	#system($command);
}