#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

include { tb_qaat }              from './workflows/qc.nf'
include { tb_d1avolo }           from './workflows/speciation_d1.nf'
include { tb_drug_resistance }from './workflows/vcfpredict.nf'
include { vcf2fasta_flow }       from './workflows/vcf2fasta.nf'
include { tb_lineage } from './workflows/fastlin.nf'
include { gen_summary_files as gen_summary_files_pos } from './workflows/summary_files.nf'
include { gen_summary_files as gen_summary_files_neg } from './workflows/summary_files.nf'

ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

if (params.help) { helpMessage(); exit(0) }

def helpMessage() {
  log.info """
========================================================================
MYPI - M Y C O B A C T E R I A L  P I P E L I N E
... (help text unchanged for brevity) ...
""".stripIndent()
}

resistance_profiler = params.resistance_profiler
human_read_removal  = params.human_read_removal

if ( params.indir == "" )           exit 1, "error: --indir is mandatory"
if ( params.filetype == "" )        exit 1, "error: --filetype is mandatory"
if ( params.filetype == "fastq" && params.readpat == "" )
  exit 1, "error: --pattern is mandatory for fastq input"
if ( params.output_dir == "" )      exit 1, "error: --output_dir is mandatory"
if ( !["fastq"].contains(params.filetype) )
  exit 1, "error: --filetype must be 'fastq'"

log.info """
========================================================================
MYPI - M Y C O B A C T E R I A L  P I P E L I N E

Parameters used:
------------------------------------------------------------------------
--indir                 ${params.indir}
--filetype              ${params.filetype}
--readpat               ${params.readpat}
--output_dir            ${params.output_dir}
--human_read_removal    ${params.human_read_removal}
--speciation_method     ${params.speciation_method}
--resistance_profiler   ${params.resistance_profiler}
--variant_caller        ${params.variant_caller}
--workflow              ${params.workflow}

Runtime data:
------------------------------------------------------------------------
Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
""".stripIndent()


process wf_complete_positive {

  input:
    val(sample_ids_list)

  output:
    val(sample_ids_list), emit: pos_sample_ids

  script:
  // Use id_list.join(' ') to convert the list ['ID1', 'ID2', 'ID3'] 
  // into the string "ID1 ID2 ID3" with space. This is what Bash needs for the 'for' loop.
    def sample_ids = sample_ids_list.join(' ')
    """
    # Access the full list of sample_ids
    echo "All sample_ids: $sample_ids"

    for sample_id in $sample_ids; do
      echo "Checking sample: \$sample_id"

      files=("${params.output_dir}/\$sample_id/output_reads_tb_or_not_tb/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/output_reads/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/megahit_assembly_reports/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/d1avolo_speciation_reports/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/antibiogram/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/tbp_parser/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/consensus/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/fastlin/ComponentComplete.txt"
      )

      all_exist=true

      for file in "\${files[@]}"; do
          if [ ! -f "\$file" ]; then
              echo "\$(date '+%Y-%m-%d %H:%M:%S') - File does not exist: \$file" >> ${params.output_dir}/\$sample_id/wf_tb_incomplete.txt
              all_exist=false
          fi
      done

      if [ "\$all_exist" = true ]; then
          echo "All files exist for \$sample_id"
          touch ${params.output_dir}/\$sample_id/WorkflowComplete.txt
          touch ${params.output_dir}/\$sample_id/wf_tb_positive_complete.txt
      else
          echo "Some files are missing (see above)" >> ${params.output_dir}/\$sample_id/wf_tb_incomplete.txt
      fi
    done 
    """
}


process wf_complete_negative {

  input:
    val(sample_ids_list)

  output:
    val(sample_ids_list), emit: neg_sample_ids

  script:
  // Use id_list.join(' ') to convert the list ['ID1', 'ID2', 'ID3'] 
  // into the string "ID1 ID2 ID3" with space. This is what Bash needs for the 'for' loop.
    def sample_ids = sample_ids_list.join(' ')
    """
    # Access the full list of sample_ids
    echo "All sample_ids: $sample_ids"

    for sample_id in $sample_ids; do
      echo "Checking sample: \$sample_id"

      files=("${params.output_dir}/\$sample_id/output_reads_tb_or_not_tb/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/output_reads/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/megahit_assembly_reports/ComponentComplete.txt"
            "${params.output_dir}/\$sample_id/d1avolo_speciation_reports/ComponentComplete.txt"
      )

      all_exist=true

      for file in "\${files[@]}"; do
          if [ ! -f "\$file" ]; then
              echo "\$(date '+%Y-%m-%d %H:%M:%S') - File does not exist: \$file" >> ${params.output_dir}/\$sample_id/wf_tb_incomplete.txt
              all_exist=false
          fi
      done

      if [ "\$all_exist" = true ]; then
          echo "All files exist for \$sample_id"
          touch ${params.output_dir}/\$sample_id/WorkflowComplete.txt
          touch ${params.output_dir}/\$sample_id/wf_tb_negative_complete.txt
      else
          echo "Some files are missing (see above)" >> ${params.output_dir}/\$sample_id/wf_tb_incomplete.txt
      fi
    done 
    """
}


workflow {
  inputdir_amended = "${params.indir}".replaceFirst(/$/, "/")
  indir = inputdir_amended

  if ( params.filetype == "fastq" ) {
    pattern = params.readpat
    reads = indir + pattern
    Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
           .ifEmpty { error "cannot find any reads matching ${pattern} in ${indir}" }
           .set{ input_files }
  }

  main:

  // qc steps
  tb_qaat(input_files)

  // lineage detection using fastlin
  tb_lineage(tb_qaat.out.qaat_output) 

  // default path with d1 as specition method and 'tb-profiler' for variant calling
  if ( params.speciation_method == "d1" && params.variant_caller == "tb-profiler") {
    tb_d1avolo(tb_lineage.out.fastlin_fq).branch{
      positive: it[5] =~ /TB\_SAMPLE\_POSITIVE/
      negative: it[5] =~ /TB\_SAMPLE\_NEGATIVE/
    }.set{ d1_tb_ch }


    //call tbprofiler for only positive samples
    tb_drug_resistance(d1_tb_ch.positive)

    // call the fasta file generation
    vcf2fasta_flow(tb_drug_resistance.out.targets_vcf_for_consensus)

    //collect all the sample_ids from negative branch and pass it as list to wf_complete_negative
    //ch_list_of_neg_ids = d1_tb_ch.negative.collect(flat : false).map { it.collect { element -> element[0] } }

    // collect all negative ids from branching channel
    ch_neg_ids = d1_tb_ch.negative.map { element -> element[0] }
    // get another channel with actual files in it
    ch_files = tb_d1avolo.out.speciation_fqs
    // filter the channel for just the elements from negative samples
    ch_filtered_negatives = ch_neg_ids.join(ch_files)
    // get only the ids and forget about the rest of the elements and pass that to wf_complete_negative()
    ch_list_of_neg_ids = ch_filtered_negatives.map { it[0] }.collect()
    wf_complete_negative(ch_list_of_neg_ids)
    gen_summary_files_neg(wf_complete_negative.out.neg_sample_ids)

    //   Create a channel that only emits the ID list from vcf2fasta_flow
    ch_list_of_pos_ids = vcf2fasta_flow.out.fasta_sample_id.collect()
    wf_complete_positive(ch_list_of_pos_ids)
    gen_summary_files_pos(wf_complete_positive.out.pos_sample_ids)

  }
}

workflow.onComplete {
  if ( workflow.success ) {
    log.info """
    ===========================================
    ${ANSI_GREEN}Workflow completed successfully
    """.stripIndent()
  } else {
    log.info """
    ===========================================
    ${ANSI_RED}Finished with errors${ANSI_RESET}
    """.stripIndent()
  }
}
