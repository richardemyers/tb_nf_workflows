# tb_nf_workflows

M Y C O B A C T E R I A L  P I P E L I N E

Cleans and QCs reads with fastp and FastQC, removes human/viral reads using hostile/tb_or_not_tb,
does species classification with Mykrobe and resistance profiling with tb-profiler

Takes as input one directory containing pairs of fastq(.gz) or bam files.
Produces as output one directory per sample, containing the relevant reports & a pair of cleaned fastqs.

Mandatory and conditional parameters:
------------------------------------------------------------------------
--input_dir           Directory containing fastq OR bam files. Workflow will process one or the other, so don't mix
--filetype            File type in input_dir. One of either "fastq" or "bam". fastq files can be gzipped and do not
                      have to literally take the form "*.fastq"; see --pattern
--pattern             Regex to match files in input_dir, e.g. "*_R{1,2}.fq.gz". Only mandatory if --filetype is "fastq"
--output_dir          Output directory, in which will be created subdirectories matching base name of fastq/bam files
-resistance_profiler Tool to profile resistance with. At the moment options are "tb-profiler", tbt-amr or "none"



Examples:
------------------------------------------------------------------------
nextflow run tb-typing.SMED.nf --filetype fastq --input_dir fq_dir --pattern "*_{1,2}.fastq.gz" --output_dir output_dir


To note: current code includes clockwork as the assembler - this is due to be replaced in the live version with tb-profiler
