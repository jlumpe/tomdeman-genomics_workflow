#!/bin/bash

# Stop on first error
set -e

# Get directory of this script - don't necessarily want files to be in same directory as scripts
# Should support symlinks
# (readlink isn't POSIX - so may not be entirely portable. See
# https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself)
DIR=`dirname $(readlink -f $0)`

# These are still hard-coded for now - change them
BUSCO_LINEAGE=/path/to/BUSCO/lineage_data/species_odb9

# Inform BUSCO of the config file location using environment variable
export BUSCO_CONFIG_FILE=$DIR/busco_config.ini

# Create output directories
mkdir PhiX_free_reads
mkdir trimmed_PhiX_free_reads
mkdir assemblies
mkdir assemblies_500
mkdir Plasmid_assemblies
mkdir Prokka_stats
mkdir BUSCO_stats

# Run bbduk
perl $DIR/run_bbduk.pl $DIR/phiX174.fasta raw_reads/*
mv raw_reads/*_noPhiX PhiX_free_reads

# Run trimmomatic
perl $DIR/run_trimmomatic.pl $DIR/adapters.fasta PhiX_free_reads/*
mv PhiX_free_reads/*paired_trimmed* trimmed_PhiX_free_reads

# Run SPAdes
perl $DIR/run_SPAdes.pl trimmed_PhiX_free_reads/*
perl $DIR/run_plasmidSPAdes.pl trimmed_PhiX_free_reads/*

#move assemblies to different folder
for file in *_assembly/scaffolds.fasta; do new="$(echo "$file" | cut -d '_' -f 1)".scaffolds.fasta; cp "$file" "assemblies/$new"; done
for file in *_PLASMID_contigs/scaffolds.fasta; do new="$(echo "$file" | cut -d '_' -f 1)".scaffolds.fasta; cp "$file" "Plasmid_assemblies/$new"; done

#filter contigs on length
for i in assemblies/*.fasta; do perl contig_size_select.pl -low 500 $i > $i.500.fna; done
mv assemblies/*.fna assemblies_500

#Check if genomes are complete gene content wise. Make sure to add the right lineage db to -l
for i in assemblies_500/*.fna

do
	location="$(echo "$i" | cut -d '.' -f 1)"
	name="$(echo "$location" | cut -d '/' -f 2)"
	run_BUSCO.py -i $i -o BUSCO_$name -l $BUSCO_LINEAGE -m geno -c 12
done

#take all BUSCO summaries and plot data
for file in run_BUSCO_*/*.txt; do cp $file BUSCO_stats; done
# This is a stupid script name but it's the default for BUSCO installation
generate_plot.py -wd BUSCO_stats

#annotate HQ contigs
for i in assemblies_500/*.fna

do
	location="$(echo "$i" | cut -d '.' -f 1)"
	name="$(echo "$location" | cut -d '/' -f 2)"
	prokka --kingdom Bacteria --outdir prokkaDIR_$name --locustag $name $i
done

for file in prokkaDIR_*/*.txt; do cp $file Prokka_stats; done
