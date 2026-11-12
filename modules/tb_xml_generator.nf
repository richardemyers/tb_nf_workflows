process tb_xml_generator {

    label 'tb_xml_generator'

    input:
    val(sample_ids_list)

    output:
    val(sample_ids_list), emit: xml_sample_ids

    script:
  // Use id_list.join(' ') to convert the list ['ID1', 'ID2', 'ID3'] 
  // into the string "ID1 ID2 ID3" with space. This is what Bash needs for the 'for' loop.
    def sample_ids = sample_ids_list.join(' ')
    """
    components=("raw_read_QC_reports"
                "tb_or_not_tb_reports"
                "megahit_assembly_reports"
                "d1avolo_speciation_reports"
                "tbp_parser"
                "fastlin"
    )

    for sample_id in $sample_ids; do
      for component in "\${components[@]}"; do
        COMPONENT_PATH="${params.output_dir}/\${sample_id}/\${component}"

        if [ -d "\$COMPONENT_PATH" ]; then
          generate_xml_json.py\
            -i "${params.output_dir}/\${sample_id}"\
            -n "\${sample_id}"\
            -w tb-typing\
            -ct "\${component}"\
            --json_summary
        else
          echo "Skipping \${component} for \${sample_id}: Directory not found."
        fi
      done
    done

    """
}

process tb_combine_xml {

    label 'tb_combine_xml'

    input:
    val(sample_ids_list)

    output:
    val(sample_ids_list), emit: combine_xml_sample_ids

    script:
  // Use id_list.join(' ') to convert the list ['ID1', 'ID2', 'ID3']
  // into the string "ID1 ID2 ID3" with space. This is what Bash needs for the 'for' loop.
    def sample_ids = sample_ids_list.join(' ')
    """
    for sample_id in $sample_ids; do
      combine_xml.py -w ${params.workflow} -i "${params.output_dir}/\${sample_id}" -s "\${sample_id}"
    done

    """
    
}

process tb_generate_json {

    label 'tb_xml_generator'

    input:
    val(sample_ids_list)

    script:
  // Use id_list.join(' ') to convert the list ['ID1', 'ID2', 'ID3']
  // into the string "ID1 ID2 ID3" with space. This is what Bash needs for the 'for' loop.
    def sample_ids = sample_ids_list.join(' ')
    """
    for sample_id in $sample_ids; do
      generate_summary_json.py -i "${params.output_dir}/\${sample_id}" -s "\${sample_id}" -m ${params.module_version}
    done
    """
}
