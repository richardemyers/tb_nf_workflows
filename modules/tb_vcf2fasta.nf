nextflow.enable.dsl = 2

process vcf2fasta {
  label 'tbprofiler'
  tag { sid }

  // publish into: <output_dir>/<sid>/consensus/
  publishDir "${params.output_dir}/${sid}/consensus", mode: 'copy', overwrite: true, pattern: '*{.fasta,*.bed,.txt}' 
  publishDir "${params.output_dir}/${sid}/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'


  input:
    tuple val(sid), path(vcf), path(vcf_index), path(ref_fa), path(bam), path(bai)

  output:
    val(sid), emit: fasta_sample_id
    tuple path("${sid}.consensus.fasta"), path("${sid}.consensus.lowdepth_masked.fasta"), path("${sid}.consensus.polished.fasta"), path("${sid}.consensus.complex_indel_fix.fasta"), emit: v2f_fasta_files
    tuple path("${sid}.consensus.lowdepth_lt${params.lowdepth_threshold ?: 3}.bed"), path("${sid}.consensus.complex_indel_fix.bed"), path("${sid}.consensus.final.bed"), emit: v2f_bed_files
    tuple val(sid), path(bam), path(bai), emit: tbparser_bam_bai
    path "ComponentComplete.txt", emit: fasta_component_complete, optional: true
    path "${sid}_vcf2fasta_command.err", emit: vcf2fasta_err optional true
    path "${sid}_vcf2fasta_command.log", emit: vcf2fasta_log optional true

  script:
  def lowDepth = params.lowdepth_threshold ?: 3
  def mapqMask = params.mapq_mask ?: 20

  """
  set -euo pipefail

  # 0) Localize reference & index
  cp -L "${ref_fa}" local_ref.fasta
  samtools faidx local_ref.fasta

  # 0.5) Localize BAM/BAI
  ln -sf "${bam}" in.bam
  if [ -f "${bai}" ]; then ln -sf "${bai}" in.bam.bai || true; fi

  # 1) Ensure gzipped VCF
  if [[ "${vcf}" == *.vcf.gz ]]; then
    ln -sf "${vcf}" targets.vcf.gz
  else
    bcftools view -Oz -o targets.vcf.gz "${vcf}"
  fi

  # 2) Use provided index if present, otherwise (re)index
  if [ -s "${vcf_index}" ]; then
    ln -sf "${vcf_index}" index || true
  fi
  bcftools index -f --csi targets.vcf.gz || tabix -f -p vcf targets.vcf.gz

  # 3) Normalize & split multi-allelic
  bcftools norm -m-any -f local_ref.fasta -Oz -o normalized.vcf.gz targets.vcf.gz
  bcftools index -f normalized.vcf.gz

  # 4) Drop insertions (keep SNPs + deletions only)
  bcftools view -i 'strlen(ALT)<=strlen(REF)' normalized.vcf.gz -Oz -o no_ins.vcf.gz
  bcftools index -f no_ins.vcf.gz

  # 4.5) Complex indel+mismatch BED (0-based, half-open)
  complex_bed="${sid}.consensus.complex_indel_fix.bed"
  : > "\${complex_bed}"

  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CIGAR\\n' normalized.vcf.gz \
    | awk 'BEGIN{OFS="\\t"}
           \$5 != "." && ((\$5 ~ /I/ || \$5 ~ /D/) && \$5 ~ /X/) {
             ref_len=length(\$3);
             alt_len=length(\$4);
             max_len=(ref_len>alt_len)?ref_len:alt_len;
             start=\$2-1;          # 0-based
             end=start+max_len;    # half-open
             print \$1, start, end
           }' \
    >> "\${complex_bed}"

  # 5) Low-depth BED via samtools depth (0-based, half-open)
  low_depth_bed="${sid}.consensus.lowdepth_lt${lowDepth}.bed"
  samtools depth -a -Q ${mapqMask} in.bam \
    | awk -v t=${lowDepth} 'BEGIN{OFS="\\t"} \$3 < t { print \$1, (\$2-1), \$2 }' > "\${low_depth_bed}"

  # 5.5) Publish a combined "final" BED (union of complex + low-depth) - sorted earlier
  final_bed="${sid}.consensus.final.bed"
  cat "\${complex_bed}" "\${low_depth_bed}" \
    | sort -k1,1 -k2,2n -k3,3n -u > "\${final_bed}"

  # 6) Consensus with variants applied ONCE (no masking)
  bcftools consensus -f local_ref.fasta --mark-del N --missing N \
    -o ${sid}.consensus_raw.fa no_ins.vcf.gz
  awk '/^>/ {print; next} {gsub(/[^ATCGNatcgn]/, "N"); print}' \
    ${sid}.consensus_raw.fa > ${sid}.consensus.fasta

  # 7) Mask low-depth regions using low_depth_bed
  cp ${sid}.consensus.fasta tmp_lowdepth.fa
  awk '
    BEGIN { FS="\\t" }
    # read BED into per-chromosome arrays (preserves order from sorted BED)
    NR==FNR {
      chrom = \$1; start = \$2; end = \$3;
      cnt[chrom]++
      starts[chrom "|" cnt[chrom]] = start
      ends[chrom "|" cnt[chrom]]   = end
      next
    }
    # FASTA header: use only the first token after ">" as seqname (avoids description mismatches)
    /^>/ {
      print
      header = substr(\$0,2)
      split(header, hparts, " ")
      seqname = hparts[1]
      pos = 0
      idx = 1
      next
    }
    {
      line = \$0
      len = length(line)
      for (i=1; i<=len; i++) {
        pos++
        # advance interval pointer while current pos is beyond current interval end
        while ( idx <= cnt[seqname] && (pos-1) >= ends[seqname "|" idx] ) idx++
        if ( idx <= cnt[seqname] ) {
          if ( (pos-1) >= starts[seqname "|" idx] && (pos-1) < ends[seqname "|" idx] ) {
            line = substr(line,1,i-1) "N" substr(line,i+1)
          }
        }
      }
      print line
    }
  ' "\${low_depth_bed}" tmp_lowdepth.fa > ${sid}.consensus.lowdepth_masked.fasta

  # 8) Mask complex indel + low-depth regions using final_bed
  cp ${sid}.consensus.fasta tmp_complex.fa
  awk '
    BEGIN { FS="\\t" }
    NR==FNR {
      chrom = \$1; start = \$2; end = \$3;
      cnt[chrom]++
      starts[chrom "|" cnt[chrom]] = start
      ends[chrom "|" cnt[chrom]]   = end
      next
    }
    /^>/ {
      print
      header = substr(\$0,2)
      split(header, hparts, " ")
      seqname = hparts[1]
      pos = 0
      idx = 1
      next
    }
    {
      line = \$0
      len = length(line)
      for (i=1; i<=len; i++) {
        pos++
        while ( idx <= cnt[seqname] && (pos-1) >= ends[seqname "|" idx] ) idx++
        if ( idx <= cnt[seqname] ) {
          if ( (pos-1) >= starts[seqname "|" idx] && (pos-1) < ends[seqname "|" idx] ) {
            line = substr(line,1,i-1) "N" substr(line,i+1)
          }
        }
      }
      print line
    }
  ' "\${final_bed}" tmp_complex.fa > ${sid}.consensus.complex_indel_fix.fasta

  # sanitize sequences (non-ATCGN -> N)
  awk '/^>/ {print; next} {gsub(/[^ATCGNatcgn]/,"N"); print}' ${sid}.consensus.lowdepth_masked.fasta > tmp && mv tmp ${sid}.consensus.lowdepth_masked.fasta
  awk '/^>/ {print; next} {gsub(/[^ATCGNatcgn]/,"N"); print}' ${sid}.consensus.complex_indel_fix.fasta > tmp && mv tmp ${sid}.consensus.complex_indel_fix.fasta

  # 9) Pileup-based restore (unchanged)
  #bam_from_antibiogram="${params.output_dir}/${sid}/antibiogram/${sid}.bam"
  bam_from_antibiogram=${bam}
  bcftools mpileup \
    -f ${sid}.consensus.lowdepth_masked.fasta \
    -R ${sid}.consensus.complex_indel_fix.bed \
    -a AD,DP,SP,INFO/AD -d 1000000 -B \
    \${bam_from_antibiogram} -Ou | \
  bcftools call -m -Ou | \
  bcftools view -e 'TYPE=\"indel\"' -Ou | \
  bcftools norm --check-ref s -f ${sid}.consensus.lowdepth_masked.fasta -Ou | \
  bcftools view | \
  awk 'BEGIN{OFS=\"\\t\"} \$5==\".\" { \$5=\$4 } {print}' > ${sid}.restoration.vcf

  bgzip -f ${sid}.restoration.vcf
  tabix -p vcf ${sid}.restoration.vcf.gz

  # 10) Apply restoration
  bcftools consensus -f ${sid}.consensus.lowdepth_masked.fasta \
    -o ${sid}.consensus.polished.fasta ${sid}.restoration.vcf.gz

  awk '/^>/ {print; next} {gsub(/[^ATCGNatcgn]/,\"N\"); print}' \
    ${sid}.consensus.polished.fasta > ${sid}.consensus.polished.tmp && \
  mv -f ${sid}.consensus.polished.tmp ${sid}.consensus.polished.fasta

  # mark done
  if [[ -f ${sid}.consensus.fasta && -f ${sid}.consensus.lowdepth_masked.fasta ]]; then
    touch ComponentComplete.txt
  fi

  cp .command.err "${sid}_vcf2fasta_command.err"
  cp .command.log "${sid}_vcf2fasta_command.log"
  """
}

