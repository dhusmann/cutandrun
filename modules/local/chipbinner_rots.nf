process CHIPBINNER_ROTS {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename == 'versions.yml') { return null }
            if (filename == 'plots') { return 'plots' }
            if (filename.startsWith('plots')) { return "plots/${filename}" }
            if (filename.endsWith('.bed')) { return "bed/${filename}" }
            if (filename.contains('differential')) { return "differential/${filename}" }
            return filename
        }
    ]

    conda "bioconda::r-rots conda-forge::r-base=4.2.3 conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-jsonlite"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-rots:1.19.0--r42hdfd78af_0' :
        'biocontainers/r-rots:1.19.0--r42hdfd78af_0' }"

    input:
    path matrix
    path clusters
    val samples_json
    val treated
    val control
    val group
    val fdr
    val lfc
    val bootstrap
    val k_value

    output:
    tuple val(group), path("chipbinner.differential.tsv")     , emit: results
    tuple val(group), path("chipbinner.significant.bed")      , emit: significant
    tuple val(group), path("chipbinner.significant_up.bed")   , emit: up
    tuple val(group), path("chipbinner.significant_down.bed") , emit: down
    tuple val(group), path("chipbinner.summary.tsv")          , emit: summary
    tuple val(group), path("plots")                           , emit: plots
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    mkdir -p plots
    python - <<'PY'
    import json
    import csv
    samples = json.loads(r'''${samples_json}''')
    with open('chipbinner.samplesheet.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample_id', 'condition'])
        for row in samples:
            writer.writerow([row['sample_id'], row['condition']])
    PY

    Rscript ${projectDir}/bin/chipbinner_rots.R \
        --matrix ${matrix} \
        --clusters ${clusters} \
        --samplesheet chipbinner.samplesheet.csv \
        --treated ${treated} \
        --control ${control} \
        --group ${group} \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --bootstrap ${bootstrap} \
        --k_value ${k_value}

    cat chipbinner.significant_up.bed chipbinner.significant_down.bed | awk 'NF' | sort -k1,1 -k2,2n | uniq > chipbinner.significant.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed -e 's/.*R version //; s/ .*//')
    END_VERSIONS
    """.stripIndent()
}
