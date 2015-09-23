#!/usr/bin/perl
use strict;
use warnings;

my $input = shift @ARGV;
my $output = shift @ARGV;

open(my $in, "<", $input);
open(my $out, ">", $output);

my @lines;

print "Gathering IDs\n";
while (my $line = <$in>) {
	chomp $line;
	#print "line: $line\n";
	#<STDIN>;
	#my $IDunParsed = (split("_", $line))[1];
	#print "IDunParsed: $IDunParsed\n";
	#<STDIN>;
	my $IDparsed = (split(/\./, $line))[1];
	#print "IDparsed: $IDparsed\n\n";
	#<STDIN>; 
	my $sequence = <$in>;
	chomp $sequence;
	
	if (defined $lines[$IDparsed]) {
		$lines[$IDparsed] = $lines[$IDparsed] .  $line . "\n" . $sequence . "\n";
	} else {
		$lines[$IDparsed] = $line . "\n" . $sequence . "\n";
	}
}

print "printing to file\n";
for my $aLine (@lines) {
	if (defined $aLine) {
		print $out $aLine;
	}
}

close $in;
close $out;
