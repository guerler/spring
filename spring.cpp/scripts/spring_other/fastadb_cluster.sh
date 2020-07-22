#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# cutoff
my $cut = 0.90;

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [list of pdb files] [path to pdb repository]\n\n";
	print "Description: This script clusters redundant pdb chains (>$cut) within the same complex\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $inpath = $ARGV[1] . "fasta/";

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>)
{    
    # data
    my @col = split(" ", $record);
    my $id = trim($col[0]);
    if ($id ne "")
    {
        @list[$n] = $id;
        $n++;
    }
}
close (data);

# loop
my $cp = ();
my $cx = 0;
for (my $i = 0; $i < $n; $i++)
{
	# check if same pdb-file
	if ($cx != 0)
	{
		if (substr($cp[0], 0, 4) ne substr($list[$i], 0, 4))
		{
			# check subset
			for (my $j = 0; $j < $cx; $j++)			
			{
				for (my $k = $j + 1; $k < $cx; $k++)
				{
					if ($cp[$j] ne "" && $cp[$k] ne "")
					{
						# get sequence identity
						my $fa = $inpath . substr($cp[$j], 0, 2) . "/" . $cp[$j];
						my $fb = $inpath . substr($cp[$k], 0, 2) . "/" . $cp[$k];
						my $id = `./seqalign.sh $fa $fb`;

						# check cut off
						if ($id > $cut)
						{
							$cp[$k] = "";
						}
					}
				}
			}
			
			# provide output
			for (my $j = 0; $j < $cx; $j++)			
			{
				if ($cp[$j] ne "")
				{
					print $cp[$j] . "\n";
				}
			}		
			
			# reset counter
			$cx = 0;
		}
	}
	
	# check if file exists
	if (-e $inpath . substr($list[$i], 0, 2) . "/" . $list[$i])
	{	
		# add pdb-code to list
		$cp[$cx] = $list[$i];	
		$cx++;
	}
}

# additional function
sub cmd($)
{
	my $string = shift;
	print ($string . "\n");
	system($string);
}
	
# additional function
sub suffix($)
{
	my $f = shift;
	return (split(/\./,$f))[-1];
}

# additional function
sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}
# additional function
sub getname($)
{
	my $f = shift;
	return substr($f, rindex($f,"/") + 1, length($f) - rindex($f,"/") - 1);
}
# additional function
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
