#!/usr/bin/perl

# flush
select((select(STDOUT), $| = 1)[0]);

# verify
if ($#ARGV + 1 < 2)
{
	print "Parameters required: [fasta sequence database] [path to repository]\n\n";
	print "Description: This script separates multiple sequences contained in a single file\n";	
	exit();
}

# read path to files
my $inp = $ARGV[0];
my $outpath = $ARGV[1] . "fasta/";

# make directory
cmd ("mkdir -p $outpath");

# open file
open (data, $inp);
my @list = ();
my $n = 0;
while (my $record = <data>)
{    
    if (trim($record) ne "")
    {
        @list[$n] = $record;
        $n++;
    }
}
close (data);

# loop
my $i = 0;
while ($i < $n)
{
    # look for sequence
    if (substr($list[$i], 0, 1) eq ">")
    {
        # content
        my @vec     = split (" ", $list[$i]);
        
        # split name        
        my $name    = substr($vec[0], 1);
        
        # split name
        my @namevec = split ("\\|", $name);
        my $namevecsize = $#namevec + 1;

		# fix name
        if ($namevecsize > 1)
        {
			$name = $namevec[1];
		}        

        # load sequence    
        my $content = ">$name\n";

        # counter
        my $chars = 0;
        
        # add sequence information
        while ($i < $n)
        {
			# next line			
			$i++;

            # is new sequence?
            if (substr($list[$i], 0, 1) eq ">")
            {
                last;
            } else {
                my $line = trim($list[$i]);
                for (my $j = 0; $j < length($line); $j++)
                {
                    $content .= substr($line, $j, 1);
                    if (++$chars >= 50)
                    {
                        $chars = 0;
                        $content .= "\n";
                    }
                }
            }
        }
         
		# make directory
		my $outpath_final = $outpath . "/" . substr($name, 0, 2) . "/";
		if (!(-e $outpath_final))
		{
    		cmd ("mkdir -p $outpath_final");
        }
          
        # write file
        open (out, ">$outpath_final" . $name);
        print out $content;
        close (out);
        #die();
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
