#!/usr/bin/perl

#runs SPAdes assembler in batch
#Written by Tom de Man

use warnings;
use strict;
use File::Basename;

use lib dirname(__FILE__);  # Import modules from file's directory
use TomDeManGenomicsWorkflow qw(find_read_pairs);

my @files=@ARGV;

# Configurable options for SPAdes run
my $spades_opts = "-t 12 --careful";

# Get read file pairs
my %paired_files = find_read_pairs(@files);

if (!%paired_files) {
	warn "No files given";
}

# Run SPAdes for each pair of reads
foreach my $name (sort keys %paired_files){
	unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
		warn "Couldn't find matching paired end files for file starting with: $name";
		next;
	}
	print "assembling your data....\n";
	print "----------------------\n";
	print "$paired_files{$name}[0]"." <--> "."$paired_files{$name}[1]"."\n";

	my $cmd="spades.py $spades_opts --only-assembler -1 $paired_files{$name}[0] -2 $paired_files{$name}[1] -o $name"."_assembly";
	print $cmd,"\n";
	die if system($cmd);
}
