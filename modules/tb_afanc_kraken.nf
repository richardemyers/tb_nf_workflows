
process afanc {
    /**
    * @QCcheckpoint confirm that the species is tb using afanc speciation
    * parse_afanc_json.py prints tb, otherbug or FAILED to stdout
    */

    tag { sample_name }
    label 'afanc'

    publishDir "${params.output_dir}/$sample_name/afanc_speciation_reports", mode: 'copy', overwrite: 'true',pattern: '*{.json,.tsv,.csv,.log,*.tar.gz}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_afanc)

    when:
    run_afanc =~ /${sample_name}/

    output:
    tuple val(sample_name), path(fq1), path(fq2),  path("${sample_name}_report.json"),  path("${sample_name}_afanc.json"), stdout, path("species_data.tsv"), emit: speciation_fqs
    path("afanc.log"), emit: speciation_log optional true
    path("${sample_name}_afanc.tar.gz"), emit: speciation_tar optional true

    script:
    afanc_res =  "${sample_name}_afanc"
    afanc_res_tar = "${afanc_res}.tar.gz"
    afanc_report_json = "${sample_name}_afanc.json"
    error_log  = "${sample_name}_err.json"
    //Add dummy top_hit species report for clockwork 
    top_hit_report_json = "${sample_name}_report.json"

    """
    afanc screen -o ${afanc_res}  ${params.afanc_myco_db}  $fq1  $fq2 > afanc.log

    is_tb=\$(python /opt/parse_afanc_json.py -j ${sample_name}_afanc/${afanc_report_json} )

    cp ${sample_name}_afanc/${afanc_report_json} .
    tar -czf ${afanc_res_tar} ${sample_name}_afanc


    if [[ \$is_tb == 'tb' ]];
    then
      printf "TB_SAMPLE_POSITIVE_${sample_name}";
      jq -n '{"top_hit": { "name": "${params.species_of_interest}",
                        "file_paths":{
                          "ref_fa": "${params.tb_ref_fasta}",
                          "clockwork_ref_dir": "${params.tb_ref_dir}"
                        }
                    }}' > ${top_hit_report_json}

    else
      printf "TB_SAMPLE_NEGATIVE_${sample_name}";
      jq -n '{"top_hit": { "name": "",
                        "file_paths":{
                          "ref_fa": "",
                          "clockwork_ref_dir": ""
                        }
                    }}' > ${top_hit_report_json}
      jq -n --arg key "\$is_tb" '{"error": ("Report do not have tb and have instead: " + \$key)}' > ${error_log}
    fi
 
    """

    stub:
    afanc_report_json = "${sample_name}_afanc.json"

    """
    touch ${afanc_report_json}
    printf ${sample_name}
    """
}



process kraken {
    /**
    * Run kraken for additional speciation confirmation
    */

    tag { sample_name }
    label 'afanc'

    publishDir "${params.output_dir}/$sample_name/kraken2_speciation_reports", mode: 'copy', overwrite: 'true',pattern: '*{.txt, .json,.tsv,.csv,.log,*.tar.gz}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_kraken)

    when:
    run_kraken =~ /${sample_name}/

    output:
    tuple val(sample_name), path(fq1), path(fq2),  path("${sample_name}_k2_output.txt"),  path("${sample_name}_k2_report.txt"), emit: k2_speciation_out

    script:
    k2_res =  "${sample_name}_k2"
    k2_rep =  "${k2_res}_output.txt"
    k2_out = "${k2_res}_report.txt"


    """
    kraken2 --db ${params.kraken2_db} --report $k2_rep --output $k2_out --paired $fq1  $fq2

    """

    stub:
    k2_res =  "${sample_name}_k2"
    k2_rep =  "${k2_res}_output.txt"
    k2_out = "${k2_res}_report.txt"

    """
    touch ${k2_rep}
    touch ${k2_out}
    """
}
