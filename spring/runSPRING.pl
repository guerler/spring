#!/usr/bin/env perl
my $docstring=<<EOF
SPRING USAGE:
=====================
./runSPRING.pl -datadir data_dir -seqname sequence_name 
=====================
Arguments:
=====================
    -datadir: This is the directory where your input sequence "seq.fasta" is
              located. "seq.fasta" is FASTA format file with two sequences.
              When you run multiple jobs, different targets must be put
              under different folders.
    -seqname: This is the name for each target. When running multiple jobs at
              the same time, make sure to use differnt name for each target.
EOF
;

use strict;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

#### check system ####
my $system=$^O;
unless (lc($system) eq "linux")
{
    printf "Your operating system $system is unsupported\n";
    printf "Currently only 64bit linux is supported\n";
    exit();
}

#### parse commandline argument ####
my $datadir='';
my $seqname='';
my $help='';

GetOptions("datadir:s" => \$datadir, "seqname:s" => \$seqname, "h"=>\$help);
$seqname=~s/[-\s]/_/g;
die "$docstring\n" if ($help || ! $datadir || ! $seqname );

$datadir=abs_path($datadir);
printf "\ndatadir=$datadir\nseqname=$seqname\n\n";

#### parse input files ####
die "ERROR! No such folder: $datadir\n" if (!-d "$datadir");
die "ERROR! No such file: $datadir/seq.fasta\n" if (!-s "$datadir/seq.fasta");

my @seqList;
my $seqnum=-1;
open(fp,"<$datadir/seq.fasta");
while(my $line=<fp>){
    if ($line=~/^>/){
        $seqnum++;
        $seqList[$seqnum]='';
    }else{
        chomp($line);
        $seqList[$seqnum].=$line;
    }
}
close(fp);
$seqnum++;
die "ERROR! Sequence number: $seqnum\t(expected 2)\n" if ($seqnum!=2);

my $seqnameA="${seqname}A";
my $seqnameB="${seqname}B";
my $tag="${seqnameA}-${seqnameB}";
my $tagdir="$datadir/$tag";
if (!-d "$tagdir"){
    system("mkdir -p $tagdir");
    system("unlink $datadir/SPRING");
    system("ln -s $tagdir/SPRING $datadir/SPRING");
}

open(fp,">$tagdir/$tag.txt");
printf fp "$tag\n";
close(fp);

open(fp,">$tagdir/$seqnameA.seq");
printf fp ">$seqnameA\n$seqList[0]\n";
close(fp);

open(fp,">$tagdir/$seqnameB.seq");
printf fp ">$seqnameB\n$seqList[1]\n";
close(fp);

#### run SPRING ####
my $sourceDirectory=abs_path(dirname(__FILE__));
my $spring_executable="$sourceDirectory/SPRING/spring.py";
my $cmd="$spring_executable -q $tagdir/$tag.txt -iDir $datadir -hhsearch -dockMono";
printf "$cmd\n";
system("$cmd");
