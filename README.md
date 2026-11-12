# MyPi

## Contents
- [Introduction](#Introduction)
- [Quick Start](#quick-start)
- [Executors](#executors)
- [Parameters](#parameters)
- [Checkpoints](#checkpoints)
- [Appendix](#appendix)

## Introduction
MyPi is an analysis pipeline for genomic data suspected to be from Mycobacterial species. It takes a directory containing paired-end reads as input.

The pipeline conducts quality control, read cleaning, host read removal, speciation, variant calling and consensus generation. Only samples that are assigned as a member of the *Mycobacterium tuberculosis* complex (MTBC) proceed to variant calling and consensus generation. [tb-or-not-tb](https://gitlab.phe.gov.uk/bioinf/NGSS/tb-typing/tb_or_not_tb) and [d1avalo](https://gitlab.phe.gov.uk/bioinf/NGSS/tb-typing/d1avolo) are both tools developed within UKHSA, and are used for host read removal and speciation respectively. [tb-profiler](https://github.com/jodyphelan/TBProfiler) is used for variant calling.

**N.B** this pipeline is under active development, so occasionally there may be changes to the code that are not yet reflected in the README.

## Quick Start
MyPi is a containerised pipeline, written in Nextflow. This quick start assumes that the pipeline is being run on a UKHSA HPC i.e. has Nextflow and Singularity installed, and the containers have been built and stored in the required location. The container definition files are available in the [singularity_defs](https://gitlab.phe.gov.uk/bioinf/NGSS/tb-typing/singularity_defs) repository.

To run the pipeline with default settings, some environment variables need to be set. This can be done with:
```bash
module purge
module use /hpspor01/tenant1/ngsservice/modulefiles
module load slurm
module load phe/phengs-environment
module load singularity/singularity-3.4.1
module add nextflow/25.04.6.5954

module load phe/workflow_version_modules/tb-typing

export NGSSERVICE_CONTAINER_PATH=/hpspor01/tenant1/ngsservice/containers_production
```
Replacing `/hpspor01` with `/hpscol02` according to the HPC the pipeline is being run on.

**N.B** the nextflow version may also need to be changed.

A 'basic' Nextflow command is:
```bash
nextflow -log $NF_LOG_DIR/$NF_LOG_FILENAME run tb-typing.SMED.nf -work-dir $NF_WORK_DIR --indir $INPUT_DIR --readpat "*.tb-typing.ngsservice.R[12].fastq.gz" --output_dir $OUTPUT_DIR -with-trace -with-report --workflow "tb-typing.ngsservice"
```
Replacing the variables with paths accordingly.

**N.B** it is recommended that the above commands are submitted to the cluster using a batch script. 

## Executors
Currently, the pipeline has only been tested on the UKHSA SLURM clusters. To run elsewhere, modifications will likely have to be made to the `nextflow.config`. For example, amending the executor definition which is currently set to:

```nextflow
executor {
    name = 'slurm'
    exitReadTimeout = '7200 sec'
    //queueSize = "$SLURM_LARGE_QUEUE_SIZE"
    queueSize = 50
}
```
## Parameters
The following parameters have default settings in the `nextflow.config`, which can be amended.

| Option | Help | Notes |
|-|-|-|
| `--workflow` | which workflow to run | default='tb-typing.ngsservice'|
| `--module_version` | version of the containers being used | |
| `--indir` | directory containing fastqs/bam files | required |
| `--filetype` | input file format (can only be fastq) | default='fastq' |
| `--readpat` | glob pattern for input files | default='*_{1,2}.fastq.gz' |
| `--output_dir` | directory to store results | required |
| `--min_reads` | minimum read threshold for qc | default=100000|
| `--resistance_profiler` | method to use for resistotyping | default='tb-profiler' |
| `--update_tbprofiler` | if the tb-profiler reference genome should be updated | default='no' |
| `--collate` | if tb-profiler collate should be run | default='yes' |
| `--config_dir` | path to nextflow config file | default="${baseDir}/nextflow.config" |
| `--pherefs_tb` | path to tb references | default="${PHENGS_BASE_PATH}/tenant1/ngsservice/pherefs/tb" |
| `--resource_dir` | MyPi resources directory| default="${pherefs_tb}/ukhsa_lodestone/v0.9.9-gamma/resources" |
| `--tb_ref_fasta` | path to reference fasta | default="${resource_dir}/tuberculosis.fasta" |
| `--tb_ref_dir` | path to other index files | default="${resource_dir}/tuberculosis" |
| `--species_of_interest` | species reported if tb positive by afanc/mykrobe | default='Mycobacterium tuberculosis' |
| `--human_read_removal` | method for host read removal | default='tb_or_not_tb' |
| `--gem_index` | the reference used by tb_or_not_tb | default="${pherefs_tb}/ukhsa_lodestone/gem3_indexes/indexes/HVMC.v3.175.ds.gem" |
| `--speciation_method` | method to use for speciation | default='d1' (d1avalo) |
| `--d1_ref_fasta` | d1avalo-specific - the directory containing required references | default="${pherefs_tb}/ukhsa_lodestone/d1avolo_speciation/d1avolo_v2_index_final/" |
| `--fastlin_barcodes` | barcodes used for lineage assignment | default="${pherefs_tb}/ukhsa_lodestone/fastlin_barcodes/MTBC_barcodes_v3.tsv" |
| `--variant_caller` | method to use for variant calling | default='tb-profiler' |
| `--tbdb_default_reference_genome` | path to reference genome used by tbprofiler | default="${pherefs_tb}/tbprofiler_tbdb/genome.fasta" |
| `--lowdepth_threshold` | depth threshold below which a site is masked | default=5 |
| `--mapq_mask` | mapping quality threshold below which a site is masked | default=20 |

**WARNING:** While the `--human_read_removal` parameter is optional, it is highly recommended (and run by default) to reduce the risk of Patient Identifying Information (PII) being retained. This choice is the responsibility of the script-user.

## Checkpoints
There are a number of checkpoints to identify potential sample fails:
1. process tb_qaat:checkFqValidity
   - checks if the input fastq is valid with `fqtools validate`
   - if not, a `invalid_fastq_file` file will be produced with `sample did not pass fqtools validation check. ${fqtools_validate_stderr}`
2. process tb_qaat:countReads
   - checks if the sample has > 100k raw reads
   - if not, a `yield_raw_data_inadequate` file will be produced with `sample did not have not have > 100k pairs of raw reads it only contained $read_count`
3. process tb_qaat:fastp
   - checks if the sample has > 100k reads after cleaning
   - if not, a `yield_processed_inadequate` file will be produced with `After fastp sample did not have > 100k pairs of raw reads it only contained $read_count`
4. process tb_qaat:tb_or_not_tb
    - checks if the sample has > 100k reads after host read removal
    - if not, a `yield_processed_inadequate` file will be produced with `After tb_or_not_tb processing sample did not have > 100k pairs of raw reads, it only contained $num_reads`
5. process tb_d1avolo:d1_avolo
   - following the process, the workflow checks if the sample contains enough MTBC to proceed (i.e. the `${sample_name}__d1avolo_contigs_tuber.fasta` file is >= 4MB)
   - if not, the sample does not proceed further
