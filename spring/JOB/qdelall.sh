#!/usr/bin/perl -w

use strict;
my $username = `whoami`;
chomp $username;

sub prefix($)
{
	my $f = shift;
	return (split(/\./,$f))[0];
}

#my $joblst = `qstat -u $username`;
my $joblst = `qstat | grep "bgovi"`;
my @jobs = split('\n', $joblst);

foreach my $job (@jobs) {
	my @fields = split(' ',$job);
	next if( scalar @fields < 4 );
	my $jobid  = prefix($fields[0]);
	my $name   = $fields[3];
	print "removing job $jobid $name\n";
	`qdel $jobid`;
}
