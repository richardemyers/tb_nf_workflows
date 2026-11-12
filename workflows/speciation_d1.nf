// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {megahit} from '../modules/tb_megahit_d1avolo.nf' params(params)
include {d1_avolo} from '../modules/tb_megahit_d1avolo.nf' params(params)

// define workflow component
workflow tb_d1avolo {

    take:
      input_files


    main:
      megahit(input_files)
      d1_avolo(megahit.out.megahi_assembly_fasta)

    emit:
      speciation_fqs = d1_avolo.out.speciation_fqs
}