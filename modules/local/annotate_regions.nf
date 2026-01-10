process ANNOTATE_REGIONS {
    tag "${results_tsv.simpleName}"
    label 'process_medium'

    publishDir = [
        path: { meta?.publish_dir ?: task.ext.publish_dir },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        enabled: { meta?.publish_dir ?: task.ext.publish_dir ? true : false }
    ]

    conda "bioconda::bedtools=2.31.1 bioconda::bedops=2.4.41 conda-forge::python=3.11 conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(results_tsv)
    val gene_bed
    val gtf

    output:
    path "*.annotated.tsv", emit: annotated
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: results_tsv.baseName
    def gene_arg = gene_bed ? "--features ${gene_bed}" : ''
    def gtf_arg = (!gene_bed && gtf) ? "--gtf ${gtf}" : ''
    """
    python ${projectDir}/bin/annotate_regions.py \
        --input ${results_tsv} \
        --output ${prefix}.annotated.tsv \
        ${gene_arg} \
        ${gtf_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """
}
