// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { tbprofiler_fastq }      from '../modules/tb_vcfpredict.nf' params(params)
include { tbp_parser }            from '../modules/tb_vcfpredict.nf' params(params)
include { tbprofiler_update_db }  from '../modules/tb_vcfpredict.nf' params(params)
include { tbprofiler_collate }    from '../modules/tb_vcfpredict.nf' params(params)


/* =========================
   TB-PROFILER PATH (direct from FASTQs)
   Input:
     - speciation_fqs : (sid, fq1, fq2, speciation_json, report, proceed_flag)
   Emits:
     - targets_vcf_for_consensus : (sid, vcf, bam, bai) for vcf2fasta_flow
   ========================= */
workflow tb_drug_resistance {
  take:
    speciation_fqs

  main:

    if ( params.resistance_profiler == "tb-profiler" ){
      if (params.update_tbprofiler == "yes"){
        tbprofiler_update_db(params.tbdb_default_reference_genome)
      }

      tbprofiler_fastq(speciation_fqs)

      // tbprofiler_fastq.out.tbprofiler_json yields:
      // (sid, tbprofiler_json, report_json, bam, bai, vcf)
      tbprofiler_json   = tbprofiler_fastq.out.tbprofiler_json

      // Build the 5‑tuple for tbp_parser directly from the tbprofiler_fastq output
      tbp_parser_input = tbprofiler_json.map { sid, tbjson, bam, bai, vcf ->
        tuple(sid, tbjson, bam, bai, vcf)
      }
      tbp_parser(tbp_parser_input)
      
      // // Build (sid, vcf, bam, bai) for vcf2fasta (tb-profiler path only)
      // targets_vcf_for_consensus = tbprofiler_json.map { sid, j, bam, bai, vcf ->
      //   tuple(sid, vcf, bam, bai)
      // }


      if (params.collate == "yes"){
        collated_jsons = tbprofiler_fastq.out.collate_json.collect()
        tbprofiler_collate(collated_jsons)
      }

      targets_vcf_for_consensus = tbp_parser.out.tbp_vcf_bam_bai
  
  }

  emit:
    targets_vcf_for_consensus

}

