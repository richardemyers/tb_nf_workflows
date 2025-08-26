#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

include { tb_qaat } from './workflows/qc.nf'
include { tb_mykrobe } from './workflows/speciation_mykrobe.nf'
include { tb_afanc } from './workflows/speciation_afanc.nf'
include { tb_clockwork } from './workflows/clockwork.nf'
include { vcfpredict_clockwork } from './workflows/vcfpredict.nf'
include { vcfpredict_tbprofiler } from './workflows/vcfpredict.nf'


/*
 ANSI escape codes to allow colour-coded output messages
 This code is from https://github.com/angelovangel
 */

ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
log.info """
========================================================================
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
========================================================================
"""
.stripIndent()
}


resistance_profiler = params.resistance_profiler
human_read_removal = params.human_read_removal


// confirm that mandatory parameters have been set and that the conditional parameter, --pattern, has been used appropriately
if ( params.input_dir == "" ) {
    exit 1, "error: --input_dir is mandatory (run with --help to see parameters)"
}
if ( params.filetype == "" ) {
    exit 1, "error: --filetype is mandatory (run with --help to see parameters)"
}
if ( ( params.filetype == "fastq" ) && ( params.pattern == "" ) ) {
    exit 1, "error: --pattern is mandatory if you are providing fastq input; describes files in --input_dir (e.g. \"*_R{1,2}.fastq.gz\") (run with --help to see parameters)"
}
if ( ( params.filetype == "bam" ) && ( params.pattern != "" ) ) {
    exit 1, "error: --pattern should only be set if you are providing fastq input (run with --help to see parameters)"
}
if ( params.output_dir == "" ) {
    exit 1, "error: --output_dir is mandatory (run with --help to see parameters)"
}
if ( ( params.filetype != "fastq" ) && ( params.filetype != "bam" ) ) {
    exit 1, "error: --filetype is mandatory and must be either \"fastq\" or \"bam\""
}


log.info """
========================================================================
M Y C O B A C T E R I A L  P I P E L I N E

Parameters used:
------------------------------------------------------------------------
--input_dir             ${params.input_dir}
--filetype              ${params.filetype}
--pattern               ${params.pattern}
--output_dir            ${params.output_dir}
--human_read_removal    ${params.human_read_removal}
--speciation_method     ${params.speciation_method}
--resistance_profiler   ${params.resistance_profiler}
--variant_caller        ${params.variant_caller}


Runtime data:
------------------------------------------------------------------------
Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
"""
.stripIndent()

// main workflow
workflow {

    // add a trailing slash if it was not originally provided to --input_dir
    inputdir_amended = "${params.input_dir}".replaceFirst(/$/, "/")

    indir = inputdir_amended
    numfiles = 0

    if ( params.filetype == "bam" ) {
        reads = indir + "*.bam"
        numfiles = file(reads) // count the number of files

        Channel.fromPath(reads)
               .set{ input_files }
    }

    if ( params.filetype == "fastq" ) {
        pattern = params.pattern
        reads = indir + pattern
        numfiles = file(reads) // count the number of files

        Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
	       .ifEmpty { error "cannot find any reads matching ${pattern} in ${indir}" }
	       .set{ input_files }
    }

    // main workflow
    main:

    tb_qaat(input_files)

    if ( params.speciation_method == "mykrobe" && params.variant_caller == "clockwork") {
        //pass the tb_qaat output to the next step (speciation by mykrobe)
        tb_mykrobe(tb_qaat.out.qaat_output).branch{
            positive: it[5] =~ /TB\_SAMPLE\_POSITIVE/
            negative: it[5] =~ /TB\_SAMPLE\_NEGATIVE/

        }.set{ mykrobe_tb_ch }

        mykrobe_tb_ch.positive | tb_clockwork 
    }

    if ( params.speciation_method == "afanc" && params.variant_caller == "clockwork") {
        //pass the tb_qaat output to the next step (speciation by afanc)
        tb_afanc(tb_qaat.out.qaat_output).branch{
            positive: it[5] =~ /TB\_SAMPLE\_POSITIVE/
            negative: it[5] =~ /TB\_SAMPLE\_NEGATIVE/

        }.set{ afanc_tb_ch }

        afanc_tb_ch.positive | tb_clockwork 
    }

    // Calling tb-profiler after calling clockwork
    if ( params.variant_caller == "clockwork") {
    // CLOCKWORK AND VCFPREDICT SUB-WORKFLOW
    profiler_input_vcf = tb_clockwork.out.profiler_input_vcf 
    tbp_parser_input_bam = tb_clockwork.out.tbp_parser_input_bam
    bam_path = tbp_parser_input_bam.map{it[2]}
    bai_path = tb_clockwork.out.tbp_parser_input_bai

    vcfpredict_clockwork(profiler_input_vcf, bam_path, bai_path)
    }

    // Calling tb-profiler directly by-passing clockwork
    if ( params.speciation_method == "mykrobe" && params.variant_caller == "tb-profiler") {
        //pass the tb_qaat output to the next step (speciation by mykrobe)
        tb_mykrobe(tb_qaat.out.qaat_output).branch{
            positive: it[5] =~ /TB\_SAMPLE\_POSITIVE/
            negative: it[5] =~ /TB\_SAMPLE\_NEGATIVE/

        }.set{ mykrobe_tb_ch }

        mykrobe_tb_ch.positive | vcfpredict_tbprofiler 
    }

    if ( params.speciation_method == "afanc" && params.variant_caller == "tb-profiler") {
        //pass the tb_qaat output to the next step (speciation by afanc)
        tb_afanc(tb_qaat.out.qaat_output).branch{
            positive: it[5] =~ /TB\_SAMPLE\_POSITIVE/
            negative: it[5] =~ /TB\_SAMPLE\_NEGATIVE/

        }.set{ afanc_tb_ch }

        afanc_tb_ch.positive | vcfpredict_tbprofiler 
    }
    
}

workflow.onComplete {
    if ( workflow.success ) {
        log.info """
        ===========================================
        ${ANSI_GREEN}Workflow completed successfully
        """
        .stripIndent()
    }
    else {
        log.info """
        ===========================================
        ${ANSI_RED}Finished with errors${ANSI_RESET}
        """
        .stripIndent()
    }
}
