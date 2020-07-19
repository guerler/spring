#!/usr/bin/perl -w

use strict;
my $username = `whoami`;
chomp $username;

`tentakel /bin/rm -rf /tmp/$username/`;
`tentakel /bin/mkdir -p /tmp/$username`;