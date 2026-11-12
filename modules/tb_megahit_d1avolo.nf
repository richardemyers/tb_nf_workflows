
process megahit {
    /**
    * A simple-minded module to assemble Illumina sequencing data
    *
    */

    tag { sample_name }
    label 'megahit'

    publishDir "${params.output_dir}/$sample_name/megahit_assembly_reports", mode: 'copy', overwrite: 'true',pattern: '*{.json,.tsv,.csv,.log,.fasta,.txt}'
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_megahit)

    when:
    run_megahit =~ /${sample_name}/

    output:
    tuple val(sample_name), path(fq1), path(fq2),  path("${sample_name}__megahiTB.json"),  path("${sample_name}__megahiTB.fasta"), emit: megahi_assembly_fasta
    path("megahiTB.log"), emit: megahi_log optional true
    path("ComponentComplete.txt"), emit: megahit_component_complete optional true
    path "${sample_name}_megahit_command.err", emit: megahit_err optional true
    path "${sample_name}_megahit_command.log", emit: megahit_log optional true

    script:
    """
    /opt/megahiTB -1 $fq1 -2 $fq2 --output-prefix $sample_name > megahiTB.log

    if [[ -f ${sample_name}__megahiTB.fasta ]]; then
         touch ComponentComplete.txt
    fi

    cp .command.err "${sample_name}_megahit_command.err"
    cp .command.log "${sample_name}_megahit_command.log"
  
    """
}


process d1_avolo {
    /**
    * New speciation method developed internally by Paolo Ribeca
    */

    tag { sample_name }
    label 'd1avolo'

    publishDir "${params.output_dir}/$sample_name/d1avolo_speciation_reports", mode: 'copy', overwrite: 'true',pattern: '*{.fasta,.json,.tsv,.csv,.txt,.log,.pdf,*.tar.gz}'
    publishDir "${params.output_dir}/$sample_name/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(megahi_json), path(megahi_fasta)

    output:
    tuple val(sample_name), path(fq1), path(fq2), path("${sample_name}__d1avolo.json"), path("${sample_name}__d1avolo.txt"),  stdout,  emit: speciation_fqs
    path("${sample_name}__d1avolo_*.fasta"), emit: d1avolo_fasta_files optional true
    tuple path("${sample_name}__d1avolo_human.pdf"), path("${sample_name}__d1avolo_human.txt"), path("${sample_name}__d1avolo_human.json"), emit: d1avolo_pdf_report optional true
    path("${sample_name}__d1avolo_references-from-contigs.txt"), emit: d1avolo_ref_ctgs_file optional true
    path("${sample_name}__d1avolo_contigs-from-references.txt"), emit: d1avolo_ctgs_ref_file optional true
    path("ComponentComplete.txt"), emit: d1avolo_component_complete optional true
    path "${sample_name}_d1avolo_command.err", emit: d1avolo_err optional true
    path "${sample_name}_d1avolo_command.log", emit: d1avolo_log optional true

    script:
    d1avolo_json = "${sample_name}__d1avolo.json"
    error_log  = "${sample_name}_err.json"
    // Define the threshold in bytes (4MB)
    limit = 4 * 1000 * 1000

    """
    /opt/d1avolo -i $megahi_fasta -I  ${params.d1_ref_fasta}/index --output-prefix $sample_name

    if [[ -f ${sample_name}__d1avolo.json && -f ${sample_name}__d1avolo.txt ]]; then
         touch ComponentComplete.txt
    fi

    # Check: Does json file contain M.tuberculosis
    if [[ -f ${sample_name}__d1avolo.json ]] && jq -e '.result.is_failed == "N" and (.composition | tostring | contains("M.tuberculosis"))' ${sample_name}__d1avolo.json > /dev/null 2>&1; then
        printf "TB_SAMPLE_POSITIVE_${sample_name}"
    else
        printf "TB_SAMPLE_NEGATIVE_${sample_name}"
    fi


    cp .command.err "${sample_name}_d1avolo_command.err"
    cp .command.log "${sample_name}_d1avolo_command.log"
 
    """

}
