// workflows/vcf2fasta.nf
nextflow.enable.dsl = 2

include { vcf2fasta as run_vcf2fasta } from '../modules/tb_vcf2fasta.nf'

// Ensure a VCF is indexed (.csi preferred, fall back to .tbi)
process index_vcf {
  label 'tbprofiler'
  tag { sid }

  input:
  // expect (sid, vcf, bam, bai)
  tuple val(sid), path(vcf), path(bam), path(bai)

  output:
  // keep bam/bai in the tuple
  tuple val(sid), path(vcf), path(index), path(bam), path(bai)

  script:
  """
  set -euo pipefail
  bcftools index -f --csi "${vcf}" || tabix -f -p vcf "${vcf}"
  if [ -f "${vcf}.csi" ]; then
    ln -sf "${vcf}.csi" index
  else
    ln -sf "${vcf}.tbi" index
  fi
  """
}

workflow vcf2fasta_flow {
  take:
    // tuples from caller: (sid, vcf, bam, bai)
    targets_vcf_ch

  main:
    def ref_path = params.tbdb_default_reference_genome

    // 1) index VCFs, preserving bam/bai
    def indexed = targets_vcf_ch | index_vcf

    // 2) pad with ref → (sid, vcf, vcf_index, ref, bam, bai)
    def padded = indexed.map { sid, vcf, idx, bam, bai ->
      tuple(sid, vcf, idx, file(ref_path), bam, bai)
    }

    // 3) run the consensus module
    run_vcf2fasta(padded)

    fasta_sample_id = run_vcf2fasta.out.fasta_sample_id
    //fasta_process_complete = run_vcf2fasta.out.fasta_component_complete
    tbparser_bam_bai = run_vcf2fasta.out.tbparser_bam_bai
    consensus_fasta_files = run_vcf2fasta.out.v2f_fasta_files
    
  
  emit:
    fasta_sample_id
    tbparser_bam_bai
    consensus_fasta_files

}


