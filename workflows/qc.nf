// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../modules/tb_qaat.nf' params(params)
include {countReads} from '../modules/tb_qaat.nf' params(params)
include {fastp} from '../modules/tb_qaat.nf' params(params)
include {fastQC} from '../modules/tb_qaat.nf' params(params)
include {hostile} from '../modules/tb_hostile.nf' params(params)
include {tb_or_not_tb} from '../modules/tb_or_not_tb.nf' params(params)

// define workflow component
workflow tb_qaat {

    take:
      input_files


    main:

      if ( params.filetype == "bam" ) {

          checkBamValidity(input_files)

          bam2fastq(checkBamValidity.out.checkValidity_bam)

          countReads(bam2fastq.out.bam2fastq_fqs)

      }

      if ( params.filetype == "fastq" ) {

          checkFqValidity(input_files)

          countReads(checkFqValidity.out.checkValidity_fqs)

      }

      fastp(countReads.out.countReads_fqs)

      fastQC(fastp.out.fastp_fqs)

      if ( params.human_read_removal == "hostile" ) {
      // add hostile
      hostile(fastp.out.fastp_fqs)
      qaat_output = hostile.out.hostile_fqs
      }

      if ( params.human_read_removal == "tb_or_not_tb" ) {
      // add tb_or_not_tb, replacing hostile step
      tb_or_not_tb(fastp.out.fastp_fqs, params.gem_index)
      qaat_output = tb_or_not_tb.out.tbntb_fqs
      }

    emit:
      qaat_output
}
