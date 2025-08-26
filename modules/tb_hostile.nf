process hostile {

    tag { sample_name }
    label 'hostile'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name/hostile_reports", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/$sample_name/output_reads_hostile", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)


    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.clean_1.fastq.gz"), path("${sample_name}_cleaned_2.clean_2.fastq.gz"), stdout, emit: hostile_fqs
    path("${sample_name}_hostile.json", emit: hostile_json)


    script:
    clean_fq1  = "${sample_name}_cleaned_1.clean_1.fastq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.clean_2.fastq.gz"
    hostile_json = "${sample_name}_hostile.json"

    """
    hostile clean --fastq1 $fq1 --fastq2 $fq2 --index /mnt/tb/ukhsa_lodestone/hostile/human-t2t-hla > $hostile_json
    printf ${sample_name}
    """

    stub:
    clean_fq1  = "${sample_name}_cleaned_1.clean_1.fastq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.clean_2.fastq.gz"
    hostile_json = "${sample_name}_fastp.json"

    """
    touch ${clean_fq1}
    touch ${clean_fq2}
    touch ${hostile_json}
    """
}