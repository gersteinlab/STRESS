#!/usr/bin/perl

##
## Author: Simon Mitternacht 2010 
##

# Takes normal-modes from stdin and outputs the 20 first
# modes to stdout. An optional command line argument changes the
# number of modes to be output.


if (scalar @ARGV) {
    $limit = shift @ARGV;
} else {
    $limit = 20;
}
$count = 0;
while (<>) {
    ++$count if (/BEGIN/);
    last if ($count > $limit);
    print $_;
}
