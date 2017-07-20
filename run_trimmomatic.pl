#!/usr/bin/perl

#Runs trimmomatic-0.33 in batch
#Written by Tom de Man

use warnings;
use strict;
use File::Basename;

my $path;

# First argument is path to adapters, rest are read files
my ($adaptors, @files) = @ARGV;

# Base command for trimmomatic
# my $trimmomatic_cmd = "java -jar trimmomatic-0.33.jar";
my $trimmomatic_cmd = "trimmomatic";

# Common options to each run of trimmomatic - configure here
my $options1 = "-phred33 -threads 12";
my $options2 = "ILLUMINACLIP:$adaptors:2:20:10:8:TRUE SLIDINGWINDOW:20:30 LEADING:20 TRAILING:20 MINLEN:50";

my %paired_files;
foreach my $file (@files){
	my ($file_name,$dir)=fileparse($file);
	$path = $dir;
	if($file_name =~ /(.+)_R([1|2])_/){
		$paired_files{$1}[$2-1]=$file;
	#attempt different naming scheme
	}elsif($file_name =~ /(.+)_([1|2])/){
		$paired_files{$1}[$2-1]=$file;
	}else{
		warn "Input file does not contain '_R1_' or '_R2_' in name: $file";
	}
}

foreach my $name (sort keys %paired_files){
	unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
		warn "Couldn't find matching paired end files for file starting with: $name";
		next;
	}
	print "trimming your data....\n";
	print "----------------------\n";
	print "$paired_files{$name}[0]"." <--> "."$paired_files{$name}[1]"."\n";

	my $cmd="$trimmomatic_cmd PE $options1 $paired_files{$name}[0] $paired_files{$name}[1] $path$name"."_R1_paired_trimmed.fastq $path$name"."_R1_single_trimmed.fastq $path$name"."_R2_paired_trimmed.fastq $path$name"."_R2_single_trimmed.fastq $options2";
	print $cmd,"\n";
	die if system($cmd);
}
