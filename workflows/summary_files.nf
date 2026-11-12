// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { tb_xml_generator } from '../modules/tb_xml_generator.nf' 
include { tb_combine_xml }  from '../modules/tb_xml_generator.nf'
include { tb_generate_json } from '../modules/tb_xml_generator.nf'

workflow gen_summary_files {

    take:
     sample_ids 

    main:

      tb_xml_generator(sample_ids)
      tb_combine_xml(tb_xml_generator.out.xml_sample_ids)
      tb_generate_json(tb_combine_xml.out.combine_xml_sample_ids)

}
