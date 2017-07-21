import os
from textwrap import dedent

from snakemake_utils import config_key_or_fail


#################### GLOBAL VARS ####################


BASE_DIR = config.get('base_dir', workflow.basedir)
READ_FILE_DIR = config.get('read_dir', 'raw_reads')


# Get input files
IDS = None
RAW_READS = None

if 'reads' in config:
	# Explicitly specified in config, dictionary of file name pairs keyed by ID
	RAW_READS = {id_: (read1, read2) for id_, (read1, read2) in config.input_files.items()}


elif 'ids' in config:
	# IDs given, infer file names automatically
	IDS = list(config.ids)
	RAW_READS = {
		id_: tuple('raw_reads/{id}_{r}.fastq'.format(id=id_, r=i + 1) for i in range(2))
		for id in IDS
	}

else:
	# Find read pairs in directory automatically
	from snakemake_utils import find_read_file_pairs

	_paired_files, _unpaired_files = find_read_file_pairs(os.listdir(READ_FILE_DIR))

	if not _paired_files:
		raise RuntimeError('No paired reads found in {}'.format(READ_FILE_DIR))

	if _unpaired_files:
		print('Warning! The following files in {} were not identified as read files: {}'.format(READ_FILE_DIR, ', '.join(_unpaired_files)))

	RAW_READS = _paired_files

if IDS is None:
	if 'ids' in config:
		IDS = list(config.ids)
		RAW_READS = {id_: RAW_READS[id_] for id_ in IDS}
	else:	
		IDS = list(RAW_READS)


# Paths to data files
PHIX_FASTA = os.path.join(BASE_DIR, 'phiX174.fasta')
ADAPTERS_FASTA = os.path.join(BASE_DIR, 'adapters.fasta')

# BUSCO configuration
BUSCO_CONFIG_FILE = config.get('busco_config_file', os.path.join(BASE_DIR, 'busco_config.ini'))


#################### RULE DEFINITIONS ####################


rule first:
	run:
		print('Use rule "all" to run all steps.')


rule show_vars:
	run:
		print('\n##### Values of script variables: #####\n')
		varnames = ['BASE_DIR', 'READ_FILE_DIR', 'IDS', 'RAW_READS', 'READ_FILE_DIR']
		for name in varnames:
			print('{}={!r}'.format(name, globals()[name]))
		print()


rule remove_PhiX:
	input:
		r1=expand('raw_reads/{id}_1.fastq', id=IDS),
		r2=expand('raw_reads/{id}_2.fastq', id=IDS),
	output:
		r1=expand('PhiX_free_reads/{id}_1_noPhiX.fastq', id=IDS),
		r2=expand('PhiX_free_reads/{id}_2_noPhiX.fastq', id=IDS),
	threads: 12
	shell:
		'bbduk.sh -Xmx20g threads={threads} k=31 hdist=1 ref={PHIX_FASTA} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2}' 


rule trim_reads:
	input: rules.remove_PhiX.output
	output:
		r1_paired=expand('trimmed_PhiX_free_reads/{id}_1_paired_trimmed.fastq', id=IDS),
		r1_single=expand('trimmed_PhiX_free_reads/{id}_1_single_trimmed.fastq', id=IDS),
		r2_paired=expand('trimmed_PhiX_free_reads/{id}_2_paired_trimmed.fastq', id=IDS),
		r2_single=expand('trimmed_PhiX_free_reads/{id}_2_single_trimmed.fastq', id=IDS),
	threads: 12
	shell:
		'trimmomatic PE -phred33 -threads {threads} {input} {output.r1_paired} {output.r1_single} {output.r2_paired} {output.r2_single} ILLUMINACLIP:{ADAPTERS_FASTA}:2:20:10:8:TRUE SLIDINGWINDOW:20:30 LEADING:20 TRAILING:20 MINLEN:50'


rule assemble_spades:
	input:
		r1=rules.trim_reads.output.r1_paired,
		r2=rules.trim_reads.output.r2_paired,
	params:
		out_dir=expand('assemblies/{id}/genomic', id=IDS)
	output:
		scaffolds=expand('assemblies/{id}/genomic/scaffolds.fasta', id=IDS)
	threads: 12
	shell:
		'spades.py -t {threads} --careful --only-assembler -1 {input.r1} -2 {input.r2} -o {params.out_dir}'


rule assemble_plasmid_spades:
	input:
		r1=rules.trim_reads.output.r1_paired,
		r2=rules.trim_reads.output.r2_paired,
	params:
		out_dir=expand('assemblies/{id}/plasmid', id=IDS)
	output:
		scaffolds=expand('assemblies/{id}/plasmid/scaffolds.fasta', id=IDS)
	threads: 12
	shell:
		'spades.py -t {threads} --plasmid --only-assembler -1 {input.r1} -2 {input.r2} -o {params.out_dir}'


rule filter_contig_length:
	input: rules.assemble_spades.output.scaffolds
	params:
		min_length=500
	output: expand('assemblies/{id}/genomic_filtered_min_length.fasta', id=IDS)
	run:
		from snakemake_utils import parse_fasta_basic, write_fasta_seq
		with open(input) as fin, open(output, 'w') as fout:
			for id_, seq in parse_fasta_basic(fin):
				if len(seq) > params.min_length:
					write_fasta_seq(fout, id_, seq)
#	shell:
#		'{BASE_DIR}/contig_size_select.pl -low {params.min_length} {input} > {output}'


# Since we are cd-ing into a different directory, need to make paths absolute
rule run_BUSCO:
	input: rules.filter_contig_length.output,
	params:
		name=expand('{id}', id=IDS),
		input_abs=[os.path.abspath(o) for o in rules.filter_contig_length.output],
		lineage=os.path.abspath(config_key_or_fail(config, 'busco_lineage')),
		tmpdir=expand('tmp_{id}', id=IDS),
	output: expand('BUSCO_output/{id}', id=IDS)
	threads: 12
	shell:
		'\n'.join([
			'cd BUSCO_output',
			'export BUSCO_CONFIG_FILE={BUSCO_CONFIG_FILE}',
			'run_BUSCO.py -i {params.input_abs} -o {params.name} -l {params.lineage} -m geno -c {threads} -t {params.tmpdir}',
			'mv run_{params.name} {params.name}',
			'rm -r {params.tmpdir}',
		])


rule plot_BUSCO:
	input: rules.run_BUSCO.output
	output:
		png=expand('BUSCO_plots/{id}.png', id=IDS),
		rscript=expand('BUSCO_plots/{id}.R', id=IDS),
	shell:
		'\n'.join([
			'export BUSCO_CONFIG_FILE={BUSCO_CONFIG_FILE}',
			'generate_plot.py -wd {input}',
			'mv {input}/busco_figure.png {output.png}',
			'mv {input}/busco_figure.R {output.rscript}',
		])


rule annotate_prokka:
	input: rules.filter_contig_length.output
	output: expand('prokka/{id}', id=IDS)
	

rule all:
	input:
		busco_stats=rules.run_BUSCO.output,
		busco_plots=rules.plot_BUSCO.output,
		plasmid_assembly=rules.assemble_plasmid_spades.output,
