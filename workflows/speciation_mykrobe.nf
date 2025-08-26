// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {mykrobe} from '../modules/tb_mykrobe.nf' params(params)

// define workflow component
workflow tb_mykrobe {

    take:
      input_files


    main:
      mykrobe(input_files)

    emit:
      speciation_fqs = mykrobe.out.speciation_fqs
}