// modules for the tb qc workflow (derived from preprocessingModules.nf)

process checkBamValidity {
    /**
    * @QCcheckpoint confirm that samtools validates bam
    */

    tag { bam_file.getBaseName() }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'
    
    publishDir "${params.output_dir}/${bam_file.getBaseName()}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple path(bam_file)

    output:
    tuple path(bam_file), stdout, emit: checkValidity_bam
    path "${bam_file.getBaseName()}_err.json", emit: checkValidityBam_log optional true
    path "${bam_file.getBaseName()}_report.json", emit: checkValidityBam_report optional true

    script:
    error_log = "${bam_file.getBaseName()}_err.json"
    report_json = "${bam_file.getBaseName()}_report.json"

    """
    is_ok=\$(samtools quickcheck $bam_file && echo 'OK' || echo 'FAIL' )

    if [ \$is_ok == 'OK' ]; then printf \$is_ok; else echo '{"error":"bam did not pass samtools quickcheck"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${error_log} > ${report_json}; fi
    """

    stub:
    error_log = "${bam_file.getBaseName()}_err.json"

    """
    touch ${error_log}
    printf ${params.checkBamValidity_isok}
    """
}

process checkFqValidity {
    /**
    * @QCcheckpoint confirm that fqtools validates both fastqs
    */

    tag { sample_name }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path "${sample_name}_err.json", emit: checkValidity_log optional true
    path "${sample_name}_report.json", emit: checkValidity_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then printf 'OK'; else echo '{"error":"sample did not pass fqtools validation check"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${error_log} > ${report_json}; fi
    """

    stub:
    error_log  = "${sample_name}_err.json"

    """
    printf ${params.checkFqValidity_isok}
    touch ${error_log}
    """
}


process bam2fastq {
    /**
    * @QCcheckpoint none
    */

    tag { bam_file.getBaseName() }
    label 'tb_qaat'
    label 'low_memory'
    label 'normal_cpu'

    input:
    tuple path(bam_file), val(is_ok)

    when:
    is_ok == 'OK'

    output:
    tuple val("${bam_file.getBaseName()}"), path("${bam_file.getBaseName()}_1.fq.gz"), path("${bam_file.getBaseName()}_2.fq.gz"), stdout, emit: bam2fastq_fqs

    script:
    """
    samtools sort -n $bam_file -o ${bam_file.getBaseName()}.sorted.bam

    bedtools bamtofastq -i ${bam_file.getBaseName()}.sorted.bam -fq ${bam_file.getBaseName()}_1.fq -fq2 ${bam_file.getBaseName()}_2.fq

    rm ${bam_file.getBaseName()}.sorted.bam

    gzip ${bam_file.getBaseName()}_1.fq || true
    gzip ${bam_file.getBaseName()}_2.fq || true

    printf 'OK'
    """

    stub:
    """
    touch ${bam_file.getBaseName()}_1.fq.gz
    touch ${bam_file.getBaseName()}_2.fq.gz
    printf 'OK'
    """
}

process countReads {
    /**
    * @QCcheckpoint fail sample if there are < 100k raw reads
    */

    tag { sample_name }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok)

    when:
    is_ok == 'OK'

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: countReads_fqs
    path "${sample_name}_err.json", emit: countReads_log optional true
    path "${sample_name}_report.json", emit: countReads_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > 100000 )); then printf "${sample_name}"; else jq -n --arg key "\$num_reads" '{"error": ("sample did not have > 100k pairs of raw reads it only contained " + \$key)}' > ${error_log} && printf "fail"  ${error_log} > ${report_json}; fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    printf ${params.countReads_runfastp}
    touch ${error_log}
    """
}

process fastp {
    /**
    * @QCcheckpoint confirm that there > 100k reads after cleaning with fastp
    */

    tag { sample_name }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy', pattern: '*_fastp.json'
    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz' // may be overwritten if unmixing needed
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp)

    when:
    run_fastp =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    path "${sample_name}_err.json", emit: fastp_log optional true
    path "${sample_name}_report.json", emit: fastp_report optional true

    script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    fastp -i $fq1 -I $fq2 -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

    rm -rf ${fastp_html}

    num_reads=\$(fqtools count ${clean_fq1} ${clean_fq2})

    if (( \$num_reads > 100000 )); then printf "${sample_name}"; else jq -n --arg key "\$num_reads" '{"error": ("after fastp sample did not have > 100k pairs of raw reads it only contained " + \$key)}' > ${error_log} && printf "fail" && jq -s ".[0] * .[1]"  ${error_log} > ${report_json}; fi
    """

    stub:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}_err.json"

    """
    printf ${params.fastp_enoughreads}
    touch ${error_log}
    touch ${clean_fq1}
    touch ${clean_fq2}
    touch ${fastp_json}
    touch ${fastp_html}
    """
}

process fastQC {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)

    output:
    path("*", emit: fastQC_all)

    script:
    """
    cat $fq1 $fq2 > ${sample_name}.fq.gz
    fastqc ${sample_name}.fq.gz
    rm ${sample_name}.fq.gz
    """

    stub:
    """
    touch ${sample_name}_fastqc.html
    touch ${sample_name}_fastqc.zip
    """
}
