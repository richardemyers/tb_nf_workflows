// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../modules/tb_qaat.nf' params(params)
include {countReads} from '../modules/tb_qaat.nf' params(params)
include {fastp} from '../modules/tb_qaat.nf' params(params)
include {fastQC} from '../modules/tb_qaat.nf' params(params)
include {tb_or_not_tb} from '../modules/tb_or_not_tb.nf' params(params)

// define workflow component
workflow tb_qaat {

    take:
      input_files

    main:

      // for FASTQ input, validate and count reads
      if ( params.filetype == "fastq" ) {
          checkFqValidity(input_files)
          countReads(checkFqValidity.out.checkValidity_fqs)
      }

      // run fastp on the countReads output; expects (sid, fq1, fq2, run_fastp)
      fastp(countReads.out.countReads_fqs)

      // generate FastQC reports on cleaned reads
      fastQC(fastp.out.fastp_fqs)

      // human read removal: use tb_or_not_tb only
      if ( params.human_read_removal == "tb_or_not_tb" ) {
        tb_or_not_tb(fastp.out.fastp_fqs, params.gem_index)
        qaat_output = tb_or_not_tb.out.tbntb_fqs
      } else {
        // no human read removal, forward cleaned reads
        qaat_output = fastp.out.fastp_fqs
      }

    emit:
      qaat_output
}