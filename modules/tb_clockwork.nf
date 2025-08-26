// modules for the clockwork workflow

process getRefFromJSON {
    tag { sample_name }
    label 'clockwork'
    
    input:
    path(species_json)
    val(doWeAlign)
    val(sample_name)

    when:
    doWeAlign =~ /TB\_SAMPLE\_POSITIVE\_${sample_name}/
    
    output:
    stdout emit: ref_path
    
    script:
    """
    ref_string=\$(jq -r '.top_hit.file_paths.ref_fa' ${species_json})
    printf "\$ref_string"
    """
    
    stub:
    """
    printf ${baseDir}/resources/tuberculosis.fasta
    """
    
}

process alignToRef {
    /**
    * @QCcheckpoint fail if insufficient number and/or quality of read alignments to the reference genome
    */

    tag { sample_name }
    label 'clockwork'

    publishDir "${params.output_dir}/$sample_name/output_bam", mode: 'copy', overwrite: 'true', pattern: '*{.bam,.bam.bai,_alignmentStats.json}'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(report_json), path(top_hit_report_json), val(doWeAlign), path(species_data)
    path(reference_path)

    output:
    tuple val(sample_name), path("${sample_name}_report.json"), path("${sample_name}.bam"), path(reference_path), stdout, emit: alignToRef_bam
    path("${sample_name}.bam.bai", emit: alignToRef_bai)
    path("${sample_name}_alignmentStats.json", emit: alignToRef_json)
    path "${sample_name}_err.json", emit: alignToRef_log optional true
    tuple val(sample_name), path("${sample_name}_report.json"), emit: alignToRef_report

    script:
    bam = "${sample_name}.bam"
    bai = "${sample_name}.bam.bai"
    stats = "${sample_name}.stats"
    stats_json = "${sample_name}_alignmentStats.json"
    report_json = "${sample_name}_report.json"
    error_log = "${sample_name}_err.json"

    """
    echo $reference_path

    minimap2 -ax sr $reference_path -t ${task.cpus} $fq1 $fq2 | samtools fixmate -m - - | samtools sort -T tmp - | samtools markdup --reference $reference_path - minimap.bam

    java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups INPUT=minimap.bam OUTPUT=${bam} RGID=${sample_name} RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample

    samtools index ${bam} ${bai}
    samtools stats ${bam} > ${stats}

    continue='yes'
    touch ${stats_json}

    if [ \$continue == 'yes' ]; then printf "NOW_VARCALL_${sample_name}"; elif [ \$continue == 'no' ]; then echo '{"error":"insufficient number and/or quality of read alignments to the reference genome"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${error_log} ${sample_name}_report_previous.json > ${report_json}; fi
    """

    stub:
    bam = "${sample_name}.bam"
    bai = "${sample_name}.bam.bai"
    stats = "${sample_name}.stats"
    stats_json = "${sample_name}_alignmentStats.json"
    out_json = "${sample_name}_report.json"
    error_log = "${sample_name}_err.json"

    """
    touch ${sample_name}.fa
    touch ${bam}
    touch ${bai}
    touch ${stats}
    touch ${stats_json}
    touch ${out_json}
    touch ${error_log}
    printf ${params.alignToRef_doWeVarCall}
    """
}

process callVarsMpileup {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'clockwork'

    publishDir "${params.output_dir}/$sample_name/output_vcfs", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample_name), path(report_json), path(bam), path(ref), val(doWeVarCall)

    when:
    doWeVarCall =~ /NOW\_VARCALL\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.bcftools.vcf"), emit: mpileup_vcf

    script:
    bcftools_vcf = "${sample_name}.bcftools.vcf"

    """
    bcftools mpileup -Ou -a 'INFO/AD' -f ${ref} ${bam} | bcftools call --threads ${task.cpus} -vm -O v -o ${bcftools_vcf}
    """

    stub:
    bcftools_vcf = "${sample_name}.bcftools.vcf"

    """
    touch ${bcftools_vcf}
    """
}

process getRefCortex {
    tag { sample_name }
    label 'clockwork'
    
    input:
    tuple val(sample_name), path(report_json), path(bam), path(ref), val(doWeVarCall)

    when:
    doWeVarCall =~ /NOW\_VARCALL\_${sample_name}/
    
    output:
    stdout
    
    script:
    """
    ref_dir=\$(jq -r '.top_hit.file_paths.clockwork_ref_dir' ${report_json})
    echo "\$ref_dir"
    """

    stub:
    """
    echo ${baseDir}/resources/tuberculosis
    """
    
    
}

process callVarsCortex {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'clockwork'


    publishDir "${params.output_dir}/$sample_name/output_vcfs", mode: 'copy', pattern: '*.vcf'
    
    input:
    tuple val(sample_name), path(report_json), path(bam), path(ref), val(doWeVarCall)
    path(ref_dir)

    when:
    doWeVarCall =~ /NOW\_VARCALL\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.cortex.vcf"), emit: cortex_vcf

    script:
    cortex_vcf = "${sample_name}.cortex.vcf"
    cortex_original = "cortex/cortex.out/vcfs/cortex_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf"
    
    """
    cp -r ${ref_dir}/* .

    clockwork cortex . ${bam} cortex ${sample_name}
    if [[ -f ${cortex_original} ]] ;
    then
        cp cortex/cortex.out/vcfs/cortex_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf ${cortex_vcf}
    else
        touch ${cortex_vcf}
    fi
    """

    stub:
    cortex_vcf = "${sample_name}.cortex.vcf"

    """
    touch ${cortex_vcf}
    """
}


process minos {
    /**
     * @QCcheckpoint check if top species is TB, if yes pass vcf to resistance profiling
     */

    tag { sample_name }
    label 'clockwork'

    publishDir "${params.output_dir}/$sample_name/output_vcfs", mode: 'copy', pattern: '*.vcf'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(report_json), path(bam), path(ref), val(doWeVarCall), path(cortex_vcf), path(bcftools_vcf)

    output:
    tuple val(sample_name), path(report_json), path(bam), path(ref), emit: minos_bam
    tuple val(sample_name), path("${sample_name}_allelic_depth.minos.vcf"), stdout, emit: minos_vcf
    tuple val(sample_name), path("${sample_name}_report.json"), emit: minos_report
    path "${sample_name}_err.json", emit: minos_log optional true

    script:
    minos_vcf = "${sample_name}.minos.vcf"
    error_log = "${sample_name}_err.json"

    """
    # 1) Generate a minimal reference FASTA for minos
    awk '{print \$1}' ${ref} > ref.fa

    # 2) Check number of variants in bcftools and cortex VCFs
    n_variants_bcf=\$(grep -i "^#" ${bcftools_vcf} | wc -l)
    n_variants_cortex=\$(grep -i "^#" ${cortex_vcf} | wc -l)

    if [[ \$n_variants_bcf == 0 ]] ; then
        grep "^#" ${cortex_vcf} > ${bcftools_vcf}
    elif [[ \$n_variants_cortex == 0 ]]; then
        grep "^#" ${bcftools_vcf} > ${cortex_vcf}
    fi

    # 3) Run minos adjudicate, create final combined VCF, and clean up
    minos adjudicate --force --reads ${bam} minos ref.fa ${bcftools_vcf} ${cortex_vcf}
    cp minos/final.vcf ${minos_vcf}
    rm -rf minos

    # 4) Index the main reference, create dict, etc.
    samtools faidx ${ref}
    samtools dict ${ref} -o ${ref.baseName}.dict
    mkdir tmp

    # 5) Fix fractional DP and AD in the minos VCF => cleaned.vcf
    awk -F '\\t' 'BEGIN { OFS="\\t" }
    !/^#/ {
        split(\$10, gt, ":");
        # DP – match with a string-based pattern
        if (gt[2] ~ "^[0-9]+\\\\.[0-9]+\$") gt[2] = int(gt[2]);
        # ALLELE_DP
        split(gt[3], ad, ",");
        for (i in ad) {
            if (ad[i] ~ "^[0-9]+\\\\.[0-9]+\$") ad[i] = int(ad[i]);
        }
        gt[3] = ad[1] "," ad[2];
        \$10 = gt[1] ":" gt[2] ":" gt[3] ":" gt[4] ":" gt[5] ":" gt[6] ":" gt[7] ":" gt[8];
    }
    { print }' ${minos_vcf} > cleaned.vcf

    # 6) Run GATK with cleaned VCF to annotate final VCF
    gatk VariantAnnotator \\
        -R ${ref} \\
        -I ${bam} \\
        -V cleaned.vcf \\
        -A DepthPerAlleleBySample \\
        -O ${sample_name}_allelic_depth.minos.vcf \\
        --tmp-dir tmp

    # 7) Determine if top species is TB and set output accordingly
    top_hit=\$(jq -r '.top_hit.file_paths.ref_fa' ${report_json})

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    if [[ \$top_hit =~ "/tuberculosis.fasta" ]]; then
        printf "CREATE_ANTIBIOGRAM_${sample_name}"
    else
        printf "CREATE_NTM_ANTIBIOGRAM_${sample_name}"
        echo '{"resistance-profiling-warning":"sample is not TB so cannot produce antibiogram using resistance profiling tools"}' \
          | jq '.' > ${error_log} &&
          jq -s ".[0] * .[1]" ${error_log} ${sample_name}_report_previous.json > ${report_json}
    fi
    """
}

process gvcf {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'clockwork'

    publishDir "${params.output_dir}/$sample_name/output_fasta", mode: 'copy', pattern: '*.fa'
    publishDir "${params.output_dir}/$sample_name/output_vcfs", mode: 'copy', pattern: '*.vcf.gz'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(report_json), path(bam), path(ref), val(doWeValCall), path(minos_vcf), val(isSampleTB)

    output:
    path("${sample_name}.gvcf.vcf.gz", emit: gvcf)
    path("${sample_name}.fa", emit: gvcf_fa)
    path "${sample_name}_err.json", emit: gvcf_log optional true
    path "${sample_name}_report.json", emit: gvcf_report optional true
    tuple val(sample_name), path(minos_vcf), path(report_json), emit: vcfmix_input
    tuple val(sample_name), path(minos_vcf), path(report_json), val(isSampleTB), emit: tbprofiler
    tuple val(sample_name), path(report_json), path(minos_vcf), val(isSampleTB), emit: gvcf_report_resistance

    script:
    gvcf = "${sample_name}.gvcf.vcf"
    gvcf_fa = "${sample_name}.fa"
    error_log = "${sample_name}_err.json"

    """
    awk '{print \$1}' ${ref} > ref.fa

    bcftools mpileup -Ou -f ref.fa ${bam} | bcftools call --threads ${task.cpus} -m -O v -o samtools_all_pos.vcf

    clockwork gvcf_from_minos_and_samtools ref.fa ${minos_vcf} samtools_all_pos.vcf ${gvcf}
    clockwork gvcf_to_fasta ${gvcf} ${gvcf_fa}

    rm samtools_all_pos.vcf
    gzip ${gvcf}

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    if [ ${params.resistance_profiler} == "none" ]; then echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${error_log} ${sample_name}_report_previous.json > ${report_json}; fi
    """

    stub:
    gvcf = "${sample_name}.gvcf.vcf.gz"
    gvcf_fa = "${sample_name}.fa"
    error_log = "${sample_name}_err.json"

    """
    touch ${gvcf}
    touch ${gvcf_fa}
    touch ${error_log}
    """
}

