
process mykrobe {
    /**
    * @QCcheckpoint confirm that the species is tb
    * parse_mykrobe_species.py prints tb, otherbug or FAILED to stdout
    */

    tag { sample_name }
    label 'mykrobe'

    publishDir "${params.output_dir}/$sample_name/mykrobe_speciation_reports", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json,_data.tsv,.csv}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_mykrobe)

    when:
    run_mykrobe =~ /${sample_name}/

    output:
    tuple val(sample_name), path(fq1), path(fq2),  path("${sample_name}_report.json"), path("${sample_name}_mykrobe_report.json"), stdout, path("species_data.tsv"), emit: speciation_fqs
    path("mykrobe.log"), emit: speciation_log optional true
    path("${sample_name}_mykrobe.tar.gz"), emit: speciation_tar optional true

    script:
    mykrobe_report_json = "${sample_name}_mykrobe_report.json"
    error_log  = "${sample_name}_err.json"
    //Add dummy top_hit species report for clockwork 
    report_json = "${sample_name}_report.json"

    """
    mykrobe predict --panel ${params.mykrobe_panel} --sample ${sample_name} --species tb --threads ${task.cpus} --format json --output ${mykrobe_report_json} -1 $fq1 $fq2
    is_tb=\$(python /opt/parse_mykrobe_species.py -j ${mykrobe_report_json} )

    if [[ \$is_tb == 'tb' ]];
    then
      printf "TB_SAMPLE_POSITIVE_${sample_name}";
      jq -n '{"top_hit": { "name": "${params.species_of_interest}",
                        "file_paths":{
                          "ref_fa": "${params.tb_ref_fasta}",
                          "clockwork_ref_dir": "${params.tb_ref_dir}"
                        }
                    }}' > ${report_json}

    else
      printf "TB_SAMPLE_NEGATIVE_${sample_name}";
      jq -n '{"top_hit": { "name": "",
                        "file_paths":{
                          "ref_fa": "",
                          "clockwork_ref_dir": ""
                        }
                    }}' > ${report_json}
      jq -n --arg key "\$is_tb" '{"error": ("Report do not have tb and have instead: " + \$key)}' > ${error_log}
    fi
 
    """

    stub:
    mykrobe_report_json = "${sample_name}_mykrobe_report.json"

    """
    touch ${mykrobe_report_json}
    printf ${sample_name}
    """
}
