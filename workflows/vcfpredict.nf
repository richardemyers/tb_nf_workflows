// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {tbprofiler} from '../modules/tb_vcfpredict.nf' params(params)
include {tbprofiler_fastq} from '../modules/tb_vcfpredict.nf' params(params)
include {tbp_parser} from '../modules/tb_vcfpredict.nf' params(params)
include {tbprofiler_update_db} from '../modules/tb_vcfpredict.nf' params(params)
include {tbprofiler_collate} from '../modules/tb_vcfpredict.nf' params(params)
include {fastlin} from '../modules/tb_vcfpredict.nf' params(params)
include {snippy} from '../modules/tb_vcfpredict.nf' params(params)


// define workflow component
workflow vcfpredict_clockwork {

    take:
    profiler_input_vcf
    profiler_input_bam
    profiler_input_bai

    main:
      sample_name =  profiler_input_vcf.map{it[0]}
      if ( params.resistance_profiler == "tb-profiler"){

        //if we are local and want to match our references, run this
        if (params.update_tbprofiler == "yes"){
          tbprofiler_update_db(reference_fasta)
        }
        
        //run tb-profiler
        tbprofiler(profiler_input_vcf)

        tbprofiler_json = tbprofiler.out.tbprofiler_json

        tbp_parser(tbprofiler_json, profiler_input_bam, profiler_input_bai)

        if(params.collate == "yes"){
          collated_jsons = tbprofiler.out.collate_json.collect()
          tbprofiler_collate(collated_jsons)
        }
      } 

      //run fastlin - mixed lineage typing tool
      fastlin(sample_name)
      

}


workflow vcfpredict_tbprofiler {

    take:
        speciation_fqs //tuple of sample_name, fq1, fq2, speciation_json, mykrobe/afanc_report, do_we_proceed


    main:
      // //get just the json
      // json = speciation_fqs.map{it[3]}
      // is_tb = speciation_fqs.map{it[5]}
        sample_name = speciation_fqs.map{it[0]}


      if ( params.resistance_profiler == "tb-profiler"){

        //if we are local and want to match our references, run this
        if (params.update_tbprofiler == "yes"){
          tbprofiler_update_db(params.tbdb_default_reference_genome)
        }
        
        //run tb-profiler
        tbprofiler_fastq(speciation_fqs)

        tbprofiler_json = tbprofiler_fastq.out.tbprofiler_json
        profiler_input_bam = tbprofiler_json.map{it[3]}
        profiler_input_bai = tbprofiler_json.map{it[4]}
        profiler_input_vcf = tbprofiler_json.map{it[5]}

        tbp_parser(tbprofiler_json)
        if(params.collate == "yes"){
          collated_jsons = tbprofiler_fastq.out.collate_json.collect()
          tbprofiler_collate(collated_jsons)
        }
      } 

      //run fastlin - mixed lineage typing tool
      fastlin(sample_name, profiler_input_bam, profiler_input_bai)

      //run snippy - Rapid haploid variant calling and core genome alignment
      //snippy(sample_name, profiler_input_bam, profiler_input_bai)
      

}
