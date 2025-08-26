// modules for the tb speciation workflow (derived from tb_speciationModules.nf)


process kraken2 {
    /**
    * @QCcheckpoint if Kraken's top family classification is NOT Mycobacteriaceae, sample will not proceed further than afanc
    */

    tag { sample_name }
    label 'tb_speciation'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_kraken_report.*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads), path(software_json)
    path(database)

    when:
    enough_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: kraken2_json
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, path(software_json), emit: kraken2_fqs
    path "${sample_name}_err.json", emit: kraken2_log optional true
    path "${sample_name}_report.json", emit: kraken2_report optional true

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_cleaned_1.fq"
    nonBac_depleted_reads_2 = "${sample_name}_cleaned_2.fq"
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2
    
    parse_kraken_report2.py ${kraken2_report} ${kraken2_json} ${params.percent_threshold} ${params.n_reads_threshold} ${params.permissive}

    extract_kraken_reads.py -k ${kraken2_read_classification} -r ${kraken2_report} -s $fq1 -s2 $fq2 -o ${nonBac_depleted_reads_1} -o2 ${nonBac_depleted_reads_2} --taxid 2 --include-children --fastq-output >/dev/null

    gzip -f ${nonBac_depleted_reads_1}
    gzip -f ${nonBac_depleted_reads_2}

    rm -rf ${sample_name}_read_classifications.txt

    run_afanc=\$(jq '.afanc' ${kraken2_json})

    if [ \$run_afanc == '\"true\"' ]; then printf "${sample_name}"; else echo '{"error":"Kraken's top family hit either wasn't Mycobacteriaceae, or there were < 100k Mycobacteriaceae reads. Sample will not proceed further than afanc."}' | jq '.' > ${error_log} && printf "no" && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_cleaned_1.fq.gz"
    nonBac_depleted_reads_2 = "${sample_name}_cleaned_2.fq.gz"
    error_log = "${sample_name}_err.json"

    """
    printf ${params.kraken2_runmykrobe}
    touch ${kraken2_report}
    touch ${kraken2_json}
    touch ${kraken2_read_classification}
    touch ${nonBac_depleted_reads_1}
    touch ${nonBac_depleted_reads_2}
    touch ${error_log}
    """
}


process afanc {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'tb_speciation'
    label 'normal_cpu'
    label 'high_memory'
    label 'retry_afanc'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_afanc*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_afanc), path(software_json), path(kraken_report), path(kraken_json)
    path(afanc_myco_db)
    val(resource_dir)
    path(refseq_path)

    output:
    tuple val(sample_name), path("${sample_name}_afanc_report.json"), stdout, emit: afanc_json
    path "${sample_name}_err.json", emit: afanc_log optional true
    path "${sample_name}_report.json", emit: afanc_report optional true
    path "${sample_name}_afanc_original.json", emit: afanc_original optional true

    script:
    afanc_report = "${sample_name}_afanc_report.json"
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    if [[ ${run_afanc} =~ ${sample_name} ]]
    then
	afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 5.0 -n 1000 -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv > afanc.log
        cp ${sample_name}/${sample_name}.json ${sample_name}_afanc_original.json
	reformat_afanc_json.py ${sample_name}/${sample_name}.json
	printf ${sample_name}
    else
	afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 2.0 -n 500 -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv > afanc.log
        cp ${sample_name}/${sample_name}.json ${sample_name}_afanc_original.json
	reformat_afanc_json.py ${sample_name}/${sample_name}.json

	identify_tophit_and_contaminants2.py ${afanc_report} ${kraken_json} $refseq_path ${params.species} ${params.unmix_myco} $resource_dir null ${params.permissive} 0
        mv "${sample_name}"_species_in_sample_pass_0.json "${sample_name}"_species_in_sample.json 

	echo '{"error":"Kraken's top family hit either wasn't Mycobacteriaceae, or there were < 100k Mycobacteriaceae reads. Sample will not proceed further than afanc."}' | jq '.' > ${error_log} && printf "no" && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${sample_name}_species_in_sample.json > ${report_json}

    fi

    """

    stub:
    afanc_report = "${sample_name}_afanc_report.json"

    """
    touch ${afanc_report}
    printf ${sample_name}
    """
}


process bowtie2 {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'tb_speciation'
    label 'normal_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_myco_reads), path(software_json)
    path(index)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), path(software_json), emit: bowtie2_fqs
    path(software_json), emit: software_json

    script:
    bam = "${sample_name}.bam"
    humanfree_fq1 = "${sample_name}_cleaned_1.fq"
    humanfree_fq2 = "${sample_name}_cleaned_2.fq"

    """
    bowtie2 --very-sensitive -p ${task.cpus} -x ./${params.bowtie_index_name} -1 $fq1 -2 $fq2 | samtools view -f 4 -Shb - > ${bam}
    samtools fastq -1 ${humanfree_fq1} -2 ${humanfree_fq2} -s singleton.fq ${bam}

    rm -rf ${bam}
    rm -rf singleton.fq

    gzip -f ${humanfree_fq1}
    gzip -f ${humanfree_fq2}
    """

    stub:
    humanfree_fq1 = "${sample_name}_cleaned_1.fq"
    humanfree_fq2 = "${sample_name}_cleaned_2.fq"

    """
    touch ${humanfree_fq1}.gz
    touch ${humanfree_fq2}.gz
    """
}

