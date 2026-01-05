process DIFFBIND_RUN {
    tag "${group}.${caller}"
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/01_diffbind/${caller}/${group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bioconductor-diffbind bioconda::bioconductor-deseq2 bioconda::bioconductor-edger bioconda::bioconductor-genomicranges conda-forge::r-base=4.2.3 conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-jsonlite conda-forge::r-yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-diffbind:3.14.0--r42hdfd78af_0' :
        'biocontainers/bioconductor-diffbind:3.14.0--r42hdfd78af_0' }"

    input:
    path samplesheet
    val treated
    val control
    val group
    val caller
    val use_spikein
    val fdr
    val lfc
    val min_overlap
    val summits
    val backend
    val extra_params
    val output_prefix

    output:
    tuple val(group), val(caller), path("${output_prefix}.results.tsv")          , emit: results
    tuple val(group), val(caller), path("${output_prefix}.significant.bed")      , emit: significant
    tuple val(group), val(caller), path("${output_prefix}.significant_up.bed")   , emit: significant_up
    tuple val(group), val(caller), path("${output_prefix}.significant_down.bed") , emit: significant_down
    tuple val(group), val(caller), path("${output_prefix}.summary.tsv")          , emit: summary
    tuple val(group), val(caller), path("plots")                                 , emit: plots
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_arg = extra_params ? "--extra_params ${extra_params}" : ''
    """
    mkdir -p plots
    Rscript ${projectDir}/bin/diffbind_run.R \
        --samplesheet ${samplesheet} \
        --treated ${treated} \
        --control ${control} \
        --group ${group} \
        --caller ${caller} \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --min_overlap ${min_overlap} \
        --summits ${summits} \
        --backend ${backend} \
        --use_spikein ${use_spikein} \
        --prefix ${output_prefix} \
        ${extra_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diffbind: \$(Rscript -e 'packageVersion("DiffBind")' 2>/dev/null | tr -d '[]')
        r-base: \$(R --version 2>&1 | head -n 1 | sed -e 's/.*R version //; s/ .*//')
    END_VERSIONS
    """
}
