#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# use items
use strict;

# verify
if ($#ARGV + 1 < 4)
{
	print "Parameters required: [id] [target repository] [database]\n\n";
	print "Description: This script performs HHSearch.\n";	
	exit();
}

# read path to files
my $id   = $ARGV[0];
my $pinp = $ARGV[1] . "/hhm/" . substr($id, 0, 2) . "/";;
my $pout = $ARGV[1] . "/hhr/" . substr($id, 0, 2) . "/";
my $db   = $ARGV[2];

# print   
print "$id...\n";

# default
my $local_db = $db;

# copy database to local node?
#my $local_db_path = "/tmp/hhm/";
#my $local_db      = $local_db_path . getname($db);
#copy($db, $local_db_path);

# make output
cmd ("mkdir -p $pout");

# copy hhm to local
my $plocal = "/tmp/$id/";
cmd ("mkdir $plocal");
cmd ("cp $pinp$id.hhm $plocal");

# execute script
cmd("./hhsuite/bin/hhsearch -i $plocal$id.hhm -d $local_db -B 10000 -E 0.01");
cmd("mv -f $plocal$id.hhr $pout");

# clean local drive
cmd("rm -rf $plocal");

# copy and untar file
sub copy($)
{
    # get pars
    my $sour = shift;
    my $dest = shift;
    my $tmpfile = $dest . "/content.tar";

    # check if database is available
    my $tag = "$dest/tag.txt";
    unless (-e $tag)
    {
        cmd("mkdir -p $dest");
        cmd("echo . > $tag");
        cmd("cp $sour $tmpfile");
        cmd("tar -C $dest -xf $tmpfile");
        cmd("rm $tmpfile");
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
sub getname($)
{
    my $f = shift;
    my $start = rindex($f,"/");
    my $end   = rindex($f,".");
    return substr($f, $start, $end - $start);
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


