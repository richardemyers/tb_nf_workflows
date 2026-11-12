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

//take fastq file as input
process tbprofiler_fastq{
    tag { sample_name }
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*{.tbprofiler-out.json,.bam,.bai,.gz,.csv,.txt}', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'
    publishDir "${params.output_dir}/${sample_name}/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:  
    tuple val(sample_name), path(fq1), path(fq2), path(report_json), path(species_data), val(doWeAlign)

    output:
    tuple val(sample_name), path("${sample_name}.tbprofiler-out.json"), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), path("${sample_name}.targets.vcf.gz"),  emit: tbprofiler_json
    path("${sample_name}/${sample_name}.results.json"), emit: collate_json
    path("${sample_name}.csv"), emit: tbp_csv
    path("ComponentComplete.txt"), emit: tbprofiler_component_complete optional true
    path "${sample_name}_tbprofiler_command.err", emit: tbprofiler_err optional true
    path "${sample_name}_tbprofiler_command.log", emit: tbprofiler_log optional true


    script:
    tbprofiler_json = "${sample_name}.tbprofiler-out.json"
    tbprofiler_bam = "${sample_name}.bam"
    tbprofiler_bai = "${sample_name}.bam.bai"
    tbprofiler_csv = "${sample_name}.csv"
    tbprofiler_vcf = "${sample_name}.targets.vcf.gz"
    
    """
    mkdir tmp
    tb-profiler profile -1 ${fq1} -2 ${fq2} --csv --prefix ${sample_name} --threads 2 --temp tmp --call-whole-genome --depth 5,5  --af 0.9,0.9  --strand 1,1  --csv --platform illumina --no_samclip
    
    mv results ${sample_name}
    cp ${sample_name}/${sample_name}.results.json ${tbprofiler_json}
    cp ${sample_name}/${sample_name}.results.csv ${tbprofiler_csv}

    cp bam/${tbprofiler_bam} ${tbprofiler_bam}
    cp bam/${tbprofiler_bai} ${tbprofiler_bai}
    cp vcf/${tbprofiler_vcf} ${tbprofiler_vcf}

    if [[ -f ${tbprofiler_bam} && -f ${tbprofiler_bai} && -f  ${tbprofiler_vcf} ]]; then
         touch ComponentComplete.txt
    fi    

    cp .command.err "${sample_name}_tbprofiler_command.err"
    cp .command.log "${sample_name}_tbprofiler_command.log"

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
    
    publishDir "${params.output_dir}/${sample_name}/tbp_parser", mode: 'copy', pattern: '*{.csv,.txt}', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/tbp_parser", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'
    publishDir "${params.output_dir}/${sample_name}/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(tbprofiler_report_json), path(bam_path), path(bai_path), path(tbprofiler_vcf)

    output:
    tuple val(sample_name), path("${sample_name}.looker_report.csv"), path("${sample_name}.lims_report.csv"), path("${sample_name}.laboratorian_report.csv"), path("${sample_name}.locus_coverage_report.csv"), emit: tbp_parser_reports_csv
    tuple val(sample_name), path(tbprofiler_vcf), path(bam_path), path(bai_path),  emit: tbp_vcf_bam_bai
    path("ComponentComplete.txt"), emit: tbp_parser_component_complete optional true
    path "${sample_name}_tbparser_command.err", emit: tbparser_err optional true
    path "${sample_name}_tbparser_command.log", emit: tbparser_log optional true
    
    script:
    
    """
    tbp-parser ${tbprofiler_report_json} ${bam_path} -o '${sample_name}' --min_frequency 0.9 --sequencing_method 'Illumina NextSeq' --min_depth 5 --min_percent_coverage 0.89 --min_read_support 1

    if [[ -f ${sample_name}.looker_report.csv && -f ${sample_name}.lims_report.csv && -f  ${sample_name}.laboratorian_report.csv ]]; then
         touch ComponentComplete.txt
    fi  

    cp .command.err "${sample_name}_tbparser_command.err"
    cp .command.log "${sample_name}_tbparser_command.log"

    """

}
