process MULTIQC_DIFFERENTIAL {
    label 'process_single'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path multiqc_config
    val differential_dir

    output:
    path "differential_multiqc_report.html", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc -f --config ${multiqc_config} --filename differential_multiqc_report.html ${differential_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
