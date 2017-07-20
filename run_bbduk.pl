#!/usr/bin/perl

#Needs bbmap in path
#Runs bbduk in batch
#Written by Tom de Man

use warnings;
use strict;
use File::Basename;

use lib dirname(__FILE__);  # Import modules from file's directory
use TomDeManGenomicsWorkflow qw(find_read_pairs);

# First argument is path to phiX174.fasta, rest are sequence files
my ($DB, @files) = @ARGV;

# Common arguments to each call to bbduk
my $bbduk_args = "-Xmx20g threads=12 k=31 hdist=1 ref=$DB";

# Get read file pairs
my $path = dirname($files[0]);
my %paired_files = find_read_pairs(@files);

if (!%paired_files) {
	warn "No files given";
}

# Run bbduk for each pair of reads
foreach my $name (sort keys %paired_files){
	unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
		warn "Couldn't find matching paired end files for file starting with: $name";
		next;
    }
	print "remove phiX174 from your data....\n";
	print "----------------------\n";
	print "$paired_files{$name}[0]"." <--> "."$paired_files{$name}[1]"."\n";

	system("bbduk.sh $bbduk_args in=$paired_files{$name}[0] in2=$paired_files{$name}[1] out=$path/$name"."_R1_noPhiX out2=$path/$name"."_R2_noPhiX");
}
