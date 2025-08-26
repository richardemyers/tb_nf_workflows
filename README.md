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

Containers in use:
------------------------------------------------------------------------
Current in use containers:
         - tb_qaat_container = "$container_path/tb_typing/tb_typing_qaat.sinimg" 
         - hostile_container = "$container_path/tb_typing/tb_typing_hostile.sinimg" - Due to be replaced with tb or not tb
         - tb_or_not_tb_container = "$container_path/tb_typing/tb_typing_tb_or_not_tb_1_1.sinimg"
         - mykrobe_container = "$container_path/tb_typing/tb_typing_mykrobe_0.4.sinimg" - placeholder during testing
         - afanc_container = "$container_path/tb_typing/tb_typing_afanc_0_3.sinimg" - placeholder during testing
         - clockwork_container = "$container_path/tb_typing/tb_typing_clockwork_v0.12.5.sinimg" - due to be removed
         - tbprofiler_container = "$container_path/tb_typing/tb_profiler_default_genome_v6.6.5.sinimg"
         - fastlin_container = "$container_path/tb_typing/tb_typing_fastlin_0_1.sinimg"
         - snippy_container = "$container_path/tb_typing/tb_typing_snippy.sif" - due to be replaced

Date to support container exection - will be made available over sftp
  
