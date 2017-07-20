# Common code for perl scripts

package TomDeManGenomicsWorkflow;

use strict;
use warnings;

use File::Basename;
use Exporter qw(import);

our @EXPORT_OK = qw(find_read_pairs);


# Take in list of read file names, return associative array of file name pairs
# keyed by the common part of the file names
sub find_read_pairs {
	my @files = @_;
	my %paired_files;

	foreach my $file (@_){
		my ($file_name,$dir)=fileparse($file);

		if($file_name =~ /(.+)_R([1|2])_/){
			$paired_files{$1}[$2-1]=$file;

		#attempt different naming scheme
		}elsif($file_name =~ /(.+)_([1|2])/){
			$paired_files{$1}[$2-1]=$file;

		}else{
			warn "Input file does not contain '_R1_' or '_R2_' in name: $file";
		}

	}

	return %paired_files;
}


# Perl is weird
1;
