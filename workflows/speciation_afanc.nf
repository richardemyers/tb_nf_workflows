// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {afanc} from '../modules/tb_afanc_kraken.nf' params(params)
include {kraken} from '../modules/tb_afanc_kraken.nf' params(params)

// define workflow component
workflow tb_afanc {

    take:
      input_files


    main:
      afanc(input_files)
      kraken(input_files)

    emit:
      speciation_fqs = afanc.out.speciation_fqs
}