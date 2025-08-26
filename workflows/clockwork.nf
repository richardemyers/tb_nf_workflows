// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {alignToRef} from '../modules/tb_clockwork.nf' params(params)
include {callVarsMpileup} from '../modules/tb_clockwork.nf' params(params)
include {callVarsCortex} from '../modules/tb_clockwork.nf' params(params)
include {minos} from '../modules/tb_clockwork.nf' params(params)
include {gvcf} from '../modules/tb_clockwork.nf' params(params)
include {getRefFromJSON} from '../modules/tb_clockwork.nf' params(params)
include {getRefCortex} from '../modules/tb_clockwork.nf' params(params)
         
// define workflow component
workflow tb_clockwork {

    take:
        speciation_fqs //tuple of sample_name, fq1, fq2, speciation_json, mykrobe/afanc_report, do_we_proceed


    main:
      //get just the json
      json = speciation_fqs.map{it[3]}
      is_tb = speciation_fqs.map{it[5]}
      sample_name = speciation_fqs.map{it[0]}

      getRefFromJSON(json, is_tb, sample_name)
      alignToRef(speciation_fqs, getRefFromJSON.out.ref_path)
      

      callVarsMpileup(alignToRef.out.alignToRef_bam)

      getRefCortex(alignToRef.out.alignToRef_bam)
      callVarsCortex(alignToRef.out.alignToRef_bam, getRefCortex.out)

      minos(alignToRef.out.alignToRef_bam
            .join(callVarsCortex.out.cortex_vcf, by: 0)
            .join(callVarsMpileup.out.mpileup_vcf, by: 0))

      gvcf(alignToRef.out.alignToRef_bam.join(minos.out.minos_vcf, by: 0))

      report_for_ntm = gvcf.out.gvcf_report_resistance
      sample_and_fqs = speciation_fqs.map{it[0,1,2]}
      profiler_input_fq = sample_and_fqs.join(report_for_ntm, by:0)

    emit:
      profiler_input_vcf = gvcf.out.tbprofiler
      profiler_input_fq = profiler_input_fq
      tbp_parser_input_bam = alignToRef.out.alignToRef_bam
      tbp_parser_input_bai = alignToRef.out.alignToRef_bai
}
