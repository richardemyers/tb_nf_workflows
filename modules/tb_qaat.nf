// modules for the tb qc workflow

process checkFqValidity {
    /**
    * @QCcheckpoint confirm that fqtools validates both fastqs
    */

    tag { sample_name }
    label 'tb_qaat'
    label 'low_memory'
    label 'low_cpu'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log,_file}'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path "invalid_fastq_file", emit: checkValidity_invalid_log optional true
    path "${sample_name}_checkFqValidity_command.err", emit: checkFqValidity_err optional true
    path "${sample_name}_checkFqValidity_command.log", emit: checkFqValidity_log optional true

    script:
    invalid_err_file = "invalid_fastq_file"

    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then 
        printf 'OK'; 
    else
        echo "sample did not pass fqtools validation check. \$is_ok" > invalid_fastq_file
        print 'FAIL';
    fi

    cp .command.err "${sample_name}_checkFqValidity_command.err"
    cp .command.log "${sample_name}_checkFqValidity_command.log"
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

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_inadequate}'
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok)

    when:
    is_ok == 'OK'

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: countReads_fqs
    tuple val(sample_name), path("yield_raw_data_inadequate"), emit: countReads_yield_log, optional: true
    path "${sample_name}_countReads_command.err", emit: countReads_err optional true
    path "${sample_name}_countReads_command.log", emit: countReads_log optional true

    script:

    """
    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > ${params.min_reads} )); then 
        printf "${sample_name}"; 
    else
        echo "sample did not have > ${params.min_reads} pairs of raw reads, it only contained \$num_reads" > yield_raw_data_inadequate
        printf "FAIL";
    fi

    cp .command.err "${sample_name}_countReads_command.err"
    cp .command.log "${sample_name}_countReads_command.log"
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
    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*{.fq.gz,.txt}'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_inadequate}'
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp)

    when:
    run_fastp =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    tuple val(sample_name), path("yield_processed_inadequate"),  emit: fastp_yield_log, optional:true
    path "${sample_name}_fastp_command.err", emit: fastp_err optional true
    path "${sample_name}_fastp_command.log", emit: fastp_log optional true
    path("ComponentComplete.txt"), emit: fastp_component_complete optional true

    script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"

    """
    fastp -i $fq1 -I $fq2 -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --thread 1 --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

    rm -rf ${fastp_html}

    num_reads=\$(fqtools count ${clean_fq1} ${clean_fq2})

    if (( \$num_reads > ${params.min_reads} )); then 
        printf "${sample_name}";
    else
        echo "After fastp processing sample did not have > ${params.min_reads} pairs of raw reads, it only contained \$num_reads" > yield_processed_inadequate
        printf "FAIL";
    fi
    
    if [[ -f ${clean_fq1} && -f ${clean_fq2} ]]; then
         touch ComponentComplete.txt
    fi

    cp .command.err "${sample_name}_fastp_command.err"
    cp .command.log "${sample_name}_fastp_command.log"
    
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
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)

    output:
    path("*", emit: fastQC_all)
    path "${sample_name}_fastqc_command.err", emit: fastqc_err optional true
    path "${sample_name}_fastqc_command.log", emit: fastqc_log optional true

    script:
    """
    cat $fq1 $fq2 > ${sample_name}.fq.gz
    fastqc ${sample_name}.fq.gz
    rm ${sample_name}.fq.gz

    cp .command.err "${sample_name}_fastqc_command.err"
    cp .command.log "${sample_name}_fastqc_command.log"

    """
}
