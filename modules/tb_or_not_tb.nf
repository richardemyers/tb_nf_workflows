process tb_or_not_tb {

    tag { sample_name }
    label 'tb_or_not_tb'

    // keep your existing publish locations/patterns
    publishDir "${params.output_dir}/$sample_name/tb_or_not_tb_reports", mode: 'copy', pattern: '*{.json,.txt}'
    publishDir "${params.output_dir}/$sample_name/output_reads_tb_or_not_tb", mode: 'copy', pattern: '*{.fastq.gz,.fastq}'
    publishDir "${params.output_dir}/$sample_name/output_reads_tb_or_not_tb", mode: 'copy', pattern: 'ComponentComplete.txt'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_inadequate}'
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
    path(gem_index)

    when:
    enough_reads =~ /${sample_name}/

    output:
    // keep the legacy interface for downstream processes
    tuple val(sample_name), path("${sample_name}_myco_1.fastq.gz"), path("${sample_name}_myco_2.fastq.gz"), stdout, emit: tbntb_fqs
    path("${sample_name}_ledger.json"), emit: tbntb_json
    tuple path("${sample_name}_host_1.fastq"), path("${sample_name}_host_2.fastq"), path("${sample_name}_labels.txt"), emit: tbntb_cont_fqs
    path("ComponentComplete.txt"), emit: tb_or_not_component_complete optional true
    path "${sample_name}_tbrnot_command.err", emit: tbrnot_err optional true
    path "${sample_name}_tbrnot_command.log", emit: tbrnot_log optional true
    tuple val(sample_name), path("yield_processed_inadequate"), emit: tbrnot_yield_log, optional: true

    /*
     * v12D writes:
     *   ${sample}__TB-or-not-TB_keep_1.fastq
     *   ${sample}__TB-or-not-TB_keep_2.fastq
     *   ${sample}__TB-or-not-TB_disc_1.fastq
     *   ${sample}__TB-or-not-TB_disc_2.fastq
     *   ${sample}__TB-or-not-TB_labels.txt
     *   ${sample}__TB-or-not-TB_report.json
     *
     * We adapt to your legacy filenames:
     *   _keep_*   -> ${sample}_myco_*.fastq.gz  (gzip for compatibility)
     *   _disc_*   -> ${sample}_host_*.fastq     (kept uncompressed as before)
     *   labels    -> ${sample}_labels.txt
     *   report    -> ${sample}_ledger.json
     */

    script:
    """
    TB-or-not-TB \\
      -i ${gem_index} \\
      -1 ${fq1} \\
      -2 ${fq2} \\
      -o ${sample_name} \\
      -t ${task.cpus ?: 8}

    # map v12D outputs to legacy names
    # keep -> myco (gzip)
    gzip -c ${sample_name}*keep_1.fastq > ${sample_name}_myco_1.fastq.gz
    gzip -c ${sample_name}*keep_2.fastq > ${sample_name}_myco_2.fastq.gz

    # disc -> host (plain fastq)
    mv ${sample_name}*disc_1.fastq ${sample_name}_host_1.fastq
    mv ${sample_name}*disc_2.fastq ${sample_name}_host_2.fastq

    # labels + report (rename report->ledger for compatibility)
    mv ${sample_name}*labels.txt    ${sample_name}_labels.txt
    cp ${sample_name}*report.json   ${sample_name}_ledger.json

    num_reads=\$(python -c "import json; print(json.load(open('${sample_name}_ledger.json'))[0]['reads_kept'])")
  
    if (( \$num_reads > ${params.min_reads} )); then
        printf "${sample_name}";
        echo \$num_reads > number_of_reads.txt
    else
        echo "After tb_or_not_tb processing sample did not have > ${params.min_reads} pairs of raw reads, it only contained \$num_reads" > yield_processed_inadequate
        printf "FAIL";
    fi


    if [[ -f ${sample_name}_myco_1.fastq.gz && -f ${sample_name}_myco_2.fastq.gz ]]; then
      touch ComponentComplete.txt
    fi

    cp .command.err "${sample_name}_tbrnot_command.err"
    cp .command.log "${sample_name}_tbrnot_command.log"

    """

}

