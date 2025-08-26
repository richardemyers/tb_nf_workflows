// modules for the vcfpredict workflow

process tbprofiler_update_db {
    tag { sample_name }
    label 'tbprofiler'

    input:
    path(reference)

    script:
    """
    mkdir tmp
    tb-profiler update_tbdb --match_ref $reference --temp tmp
    """
}

//take vcf file from clockwork as input
process tbprofiler{
    tag { sample_name }
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.tbprofiler-out.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(minos_vcf), path(report_json), val(isSampleTB)

    output:
    tuple val(sample_name), path("${sample_name}.tbprofiler-out.json"), path("${sample_name}_report.json"), emit: tbprofiler_json
    path("${sample_name}/${sample_name}.results.json"), emit: collate_json
    tuple val(sample_name), path(minos_vcf), path(report_json), emit: vcfmix_in

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    script:
    error_log = "${sample_name}_err.json"
    tbprofiler_json = "${sample_name}.tbprofiler-out.json"
    
    """
    #keep the original vcf so we can collate the output and pass it down
    cp ${minos_vcf} tmp.vcf
    bgzip ${minos_vcf}
    mv tmp.vcf ${minos_vcf}
    
    mkdir tmp
    tb-profiler profile --vcf ${minos_vcf}.gz --threads ${task.cpus} --temp tmp --prefix ${sample_name}
    
    mv results ${sample_name}
    cp ${sample_name}/${sample_name}.results.json ${tbprofiler_json}
    
    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json  ${tbprofiler_json} > ${report_json}
    """

    stub:
    """
    mkdir ${sample_name}
    touch ${sample_name}.tbprofiler-out.json
    touch ${sample_name}_report.json
    touch ${sample_name}/${sample_name}.results.json
    """
}

//take fastq file as input
process tbprofiler_fastq{
    tag { sample_name }
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*{.tbprofiler-out.json,.bam,.bai,.gz,.csv}', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(report_json), path(top_hit_report_json), val(doWeAlign), path(species_data)

    output:
    tuple val(sample_name), path("${sample_name}.tbprofiler-out.json"), path("${sample_name}_report.json"), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), path("${sample_name}.targets.vcf.gz"),  emit: tbprofiler_json
    path("${sample_name}/${sample_name}.results.json"), emit: collate_json
    path("${sample_name}.csv"), emit: tbp_csv


    script:
    error_log = "${sample_name}_err.json"
    tbprofiler_json = "${sample_name}.tbprofiler-out.json"
    tbprofiler_bam = "${sample_name}.bam"
    tbprofiler_bai = "${sample_name}.bam.bai"
    tbprofiler_csv = "${sample_name}.csv"
    tbprofiler_vcf = "${sample_name}.targets.vcf.gz"
    
    """
    mkdir tmp
    tb-profiler profile -1 ${fq1} -2 ${fq2} --csv --prefix ${sample_name} --threads ${task.cpus} --temp tmp
    
    mv results ${sample_name}
    cp ${sample_name}/${sample_name}.results.json ${tbprofiler_json}
    cp ${sample_name}/${sample_name}.results.csv ${tbprofiler_csv}

    cp bam/${tbprofiler_bam} ${tbprofiler_bam}
    cp bam/${tbprofiler_bai} ${tbprofiler_bai}
    cp vcf/${tbprofiler_vcf} ${tbprofiler_vcf}
    
    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json  ${tbprofiler_json} > ${report_json}
    """

    stub:
    """
    mkdir ${sample_name}
    touch ${sample_name}.tbprofiler-out.json
    touch ${sample_name}_report.json
    touch ${sample_name}/${sample_name}.results.json
    """
}


process tbprofiler_collate{
    label 'tbprofiler'

    publishDir "${params.output_dir}", mode: 'copy', overwrite: 'true', pattern: 'tbprofiler.variants.csv'

    input:
    path(files), stageAs: "results/*"
    
    output:
    path("tbprofiler.variants.csv")

    script:
    """
    tb-profiler collate
    """
}

process tbp_parser{
    tag { sample_name }
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/tbp_parser", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/tbp_parser", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(tbprofiler_report_json), path(sample_report_json), path(bam_path), path(bai_path), path(tbprofiler_vcf)

    output:
    tuple val(sample_name), path("${sample_name}.looker_report.csv"), path("${sample_name}.lims_report.csv"), path("${sample_name}.laboratorian_report.csv"), path("${sample_name}.percent_gene_coverage.csv"), emit: tbp_parser_reports_csv

    script:
    error_log = "${sample_name}_err.json"
    tbp_parser_csv = "${sample_name}.looker_report.csv"
    
    """
    python3 /opt/tbp-parser/tbp_parser/tbp_parser.py ${tbprofiler_report_json} ${bam_path} -o '${sample_name}' --min_depth 12 --min_frequency 0.9 --sequencing_method 'Illumina NextSeq'

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    """

    stub:
    """
    mkdir ${sample_name}
    touch ${sample_name}/${sample_name}.looker_report.csv
    """
}


process fastlin{
    tag { sample_name }
    label 'fastlin'
    
    publishDir "${params.output_dir}/${sample_name}/fastlin", mode: 'copy', pattern: '*.tsv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/fastlin", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    val(sample_name)
    path(bam_path)
    path(bai_path)


    output:
    tuple val(sample_name), path("${sample_name}.fastlin.out.tsv"), emit: fastlin_csv

    script:
    error_log = "${sample_name}_err.json"
    fastlin_out_csv = "${sample_name}.fastlin.out.tsv"    
    
    """

    fastlin -d \$PWD -b ${params.fastlin_barcodes} -o ${fastlin_out_csv}
   
    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    """

    stub:
    """
    mkdir ${sample_name}
    touch ${sample_name}/${fastlin_out_csv}
    """
}


process snippy{
    tag { sample_name }
    label 'snippy'
    
    publishDir "${params.output_dir}/${sample_name}/snippy", mode: 'copy', pattern: '*{.aligned.fa,*.bam,*.bai,*.bed,*.vcf.gz}', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/snippy", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    val(sample_name)
    path(bam_path)
    path(bai_path)


    output:
    tuple val(sample_name), path("${sample_name}.fastlin.out.tsv"), emit: snippy_outputs

    script:
    error_log = "${sample_name}_err.json"
  
    
    """
    snippy --reference ${params.tbdb_default_reference_genome} --bam ${bam_path} --cpus ${task.cpus} --outdir \$PWD/snippy_results

    snippy-core --ref ${params.tbdb_default_reference_genome} --prefix ${sample_name} --gap-char - \$PWD/snippy_results
   
    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}
    """


}