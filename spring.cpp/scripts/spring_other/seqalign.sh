#!/usr/bin/env perl

#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   nw.pl : Needlman and Wunsch
#   dynamic programming (special case, gop=0).
#   each sequence comes in a differrent file

# read two sequences from two fasta format file:
open F0, $ARGV[0];
while ( <F0>){$al0.=$_;}
close F0;

open F1, $ARGV[1];
while ( <F1>){$al1.=$_;}
close F1;

# extract the names and the sequences
@name_list0=($al0=~/>(.*)[^>]*/g);
@seq_list0 =($al0=~/>.*([^>]*)/g);

@name_list1=($al1=~/>(.*)[^>]*/g);
@seq_list1 =($al1=~/>.*([^>]*)/g);

# get rid of the newlines, spaces and numbers
foreach $seq (@seq_list0)
{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;	
}

foreach $seq (@seq_list1)
{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;	
}

# split the sequences
for ($i=0; $i<=$#name_list0; $i++)
{
	$res0[$i]=[$seq_list0[$i]=~/([a-zA-Z-]{1})/g];
}

for ($i=0; $i<=$#name_list1; $i++)
{
	$res1[$i]=[$seq_list1[$i]=~/([a-zA-Z-]{1})/g];
}

# evaluate substitutions
$len0 = $#{$res0[0]}+1;
$len1 = $#{$res1[0]}+1;

# initialize
for ($i=0; $i<=$len0; $i++){$smat[$i][0]=0.0; $tb[$i][0 ]= 1;}
for ($j=0; $j<=$len1; $j++){$smat[0][$j]=0.0; $tb[0 ][$j]=-1;}

# loop
for ($i=1; $i<=$len0; $i++)
{
	for ($j=1; $j<=$len1; $j++)
	{
		# calculate score
		if ($res0[0][$i-1] eq $res1[0][$j-1])
			{$s=1.0;}
		else
			{$s=0.0;}
		
		# calculate options
		$sub=$smat[$i-1][$j-1] + $s;
		$del=$smat[$i  ][$j-1];
		$ins=$smat[$i-1][$j  ];
		
		# select
		if ($sub>$del && $sub>$ins)
			{$smat[$i][$j]=$sub;$tb[$i][$j]=0;}
		elsif($del>$ins)
			{$smat[$i][$j]=$del;$tb[$i][$j]=-1;}
		else
			{$smat[$i][$j]=$ins;$tb[$i][$j]=1;}
	}
}

# write score
my $score = -1.0;
my $lng = $len0;
if ($lng < $len1)
{
	$lng = $len1;
}
if ($lng > 0.0)
{
	$score = $smat[$len0][$len1] / $lng;
}
print $score . "\n";
exit(1);

# retrieve alignment
$i=$len0;
$j=$len1;
$aln_len=0;
while (!($i==0 && $j==0))
{
	if ($tb[$i][$j]==0)
	{
		$aln0[$aln_len]=$res0[0][--$i];
		$aln1[$aln_len]=$res1[0][--$j];
	}
	elsif ($tb[$i][$j]==-1)
	{
		$aln0[$aln_len]='-';
		$aln1[$aln_len]=$res1[0][--$j];
	}
	elsif ($tb[$i][$j]==1)
	{
		$aln0[$aln_len]=$res0[0][--$i];
		$aln1[$aln_len]='-';
	}
	$aln_len++;	
}

# output fasta
print ">$name_list0[0]\n";
for ($i=$aln_len-1; $i>=0; $i--){print $aln0[$i];}
print "\n";
print ">$name_list1[0]\n";
for ($j=$aln_len-1; $j>=0; $j--){print $aln1[$j];}
print "\n";

