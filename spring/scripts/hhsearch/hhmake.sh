#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# NEEDS TO BE MODIFIED
my $tmpd    = "/tmp/guerler/";
my $db_path = "/library/yzhang/springdb/";

# use items
use strict;

# verify
if ($#ARGV + 1 < 2) {
	print "Parameters required: [id] [target repository directory]\n\n";
	print "Description: This script generates HHM files.\n";	
	exit();
}

# read path to files
my $f = $ARGV[0];
my $path_fasta = $ARGV[1] . "fasta/" . substr($f, 0, 2) . "/";
my $path_out = $ARGV[1] . "hhm/" . substr($f, 0, 2) . "/";
my $local = trim(`pwd`);
print $local . "\n";

# print   
print "$f...\n";
my $t = "$tmpd/$f";
my $s = "$t/$f.fasta";
my $o = "$t/$f.a3m";

#
# set par
#
my $db_name = "uniprot20_02Sep11";
my $db_tail = ".hhblits.tar.gz";
my $dr = "$local/hhsuite";
my $ds = "$dr/lib/hh/scripts";

# set env
$ENV{'HHLIB'} = "$dr/lib/hh";
$ENV{'PATH'} = "$ENV{'PATH'}:$dr/bin:$ds";

# make hhm
cmd("mkdir -p $tmpd");
cmd("mkdir -p $path_out");

# copy database
print "Preparing database ($db_path/$db_name$db_tail)...\n";
utar ("$db_path/$db_name$db_tail", "$tmpd/db");

# run
print "\nRunning query...\n";
cmd("mkdir -p $t");
cmd("cp $path_fasta/$f $s");
cmd("./hhsuite/bin/hhblits -i $s -d $tmpd/db/$db_name -oa3m $o -n 2");
chdir($ds) or die "Cant chdir to $ds $!";
cmd("./addss.pl $o") or "$!";
chdir($dr) or die "Cant chdir to $dr $!";
cmd("./bin/hhmake -i $o");
cmd("mv $t/$f.hhm $path_out");
cmd("rm -rf $t");

# copy and untar file
sub utar($)
{
    # get pars
    my $sour = shift;
    my $dest = shift;

    # check if database is available
    my $tag = "$dest/tag.txt";
    unless (-e $tag)
    {
        cmd("mkdir -p $dest");
        cmd("echo . > $tag");
        cmd("tar -C $dest -xzf $sour");
        cmd("echo .. > $tag");
    } else {
        # wait
        while (-s $tag == 2)
        {
            print ".";
            `sleep 20`;
        }
    }
}

# additional function
sub cmd($)
{
	my $string = shift;
	print $string . "\n";
	system($string);
}

# trim
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


