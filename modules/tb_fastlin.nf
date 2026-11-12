
process fastlin{
    tag { sample_name }
    label 'fastlin'
    
    publishDir "${params.output_dir}/${sample_name}/fastlin", mode: 'copy', pattern: '*{.tsv,.txt}', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/logs", mode: 'copy', overwrite: 'true', pattern: '*{.err,.log}'


    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastlin)

    when:
    run_fastlin =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.fastlin.out.tsv"), emit: fastlin_csv
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: fastlin_fq
    tuple val(sample_name), path("ComponentComplete.txt"), emit: fastlin_component_complete optional true
    path "${sample_name}_fastlin_command.err", emit: fastlin_err optional true
    path "${sample_name}_fastlin_command.log", emit: fastlin_log optional true



    script:
    fastlin_out_csv = "${sample_name}.fastlin.out.tsv"
    temp_dir =  "${sample_name}_temp"
    
    """
    mkdir -p $temp_dir
    cp $fq1 $temp_dir/
    cp $fq2 $temp_dir/

    for old_file in "$temp_dir"/*_myco_*.fastq.gz; do
        [ -e "\$old_file" ] || continue
        
        new_file=\${old_file/_myco_/_}
        
        mv "\$old_file" "\$new_file"
    done

    fastlin -d $temp_dir -b ${params.fastlin_barcodes} -o ${fastlin_out_csv}

    cp .command.err "${sample_name}_fastlin_command.err"
    cp .command.log "${sample_name}_fastlin_command.log"

    if [[ -f ${fastlin_out_csv}  ]]; then
         printf "${sample_name}";
         touch ComponentComplete.txt
    fi
        

    """
}

