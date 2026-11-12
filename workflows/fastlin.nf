// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { fastlin }               from '../modules/tb_fastlin.nf' params(params)


workflow tb_lineage {
  take:
    // fqs after qaat and tb_or_not_tb
    qc_input_files

  main:

    // Mixed lineage typing
    fastlin(qc_input_files)
    fastlin_csv = fastlin.out.fastlin_csv
    fastlin_fq = fastlin.out.fastlin_fq
    fastlin_sample_id = fastlin_csv.map { it[0] }

    fastlin_process_complete = fastlin.out.fastlin_component_complete


  emit:
    fastlin_sample_id
    fastlin_csv
    fastlin_fq
    fastlin_process_complete
}

