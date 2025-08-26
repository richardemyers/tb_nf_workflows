process tb_or_not_tb {

    tag { sample_name }
    label 'tb_or_not_tb'


    publishDir "${params.output_dir}/$sample_name/tb_or_not_tb_reports", mode: 'copy', pattern: '*{.json,.txt}'
    publishDir "${params.output_dir}/$sample_name/output_reads_tb_or_not_tb", mode: 'copy', pattern: '*{.fastq.gz,.fastq}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
    path(gem_index)

    output:
    tuple val(sample_name), path("${sample_name}_myco_1.fastq.gz"), path("${sample_name}_myco_2.fastq.gz"), stdout, emit: tbntb_fqs
    path("${sample_name}_ledger.json"), emit: tbntb_json
    tuple path("${sample_name}_host_1.fastq"), path("${sample_name}_host_2.fastq"), path("${sample_name}_labels.txt"), emit: tbntb_cont_fqs


    script:
    clean_fq1  = "${sample_name}_myco_1.fastq"
    clean_fq2  = "${sample_name}_myco_2.fastq"

    """
    /opt/TB-or-not-TB -i $gem_index -1 $fq1 -2 $fq2 -u -o $sample_name
    printf ${sample_name}

    gzip ${clean_fq1}
    gzip ${clean_fq2}

    """

    stub:
    clean_fq1  = "${sample_name}_myco_1.fastq.gz"
    clean_fq2  = "${sample_name}_myco_2.fastq.gz"
    tbntb_json = "${sample_name}_ledger.json"

    """
    touch ${clean_fq1}
    touch ${clean_fq2}
    touch ${tbntb_json}
    """
}