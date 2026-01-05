process SPAN_COMPARE {
    tag "${meta.id}"
    label 'process_long'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/${meta.group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::openjdk=21.0.2 bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://eclipse-temurin:21-jre' :
        'eclipse-temurin:21-jre' }"

    input:
    tuple val(meta), path(treated_bams), path(control_bams)
    path chrom_sizes
    path omnipeaks_jar
    val gap
    val bin
    val fdr
    val java_heap

    output:
    tuple val(meta), path("span.differential.tsv")     , emit: tsv
    tuple val(meta), path("span.differential.bed")     , emit: bed
    tuple val(meta), path("span.significant_up.bed")   , emit: up
    tuple val(meta), path("span.significant_down.bed") , emit: down
    tuple val(meta), path("span.summary.tsv")          , emit: summary
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.span_compare"
    def fdr_arg = fdr ? "--fdr ${fdr}" : ''
    def bin_arg = bin ? "--bin ${bin}" : ''
    """
    treated_list=(${treated_bams})
    control_list=(${control_bams})

    treated_input=${treated_bams}
    control_input=${control_bams}
    if [ \${#treated_list[@]} -gt 1 ]; then
        samtools merge -f treated_pool.bam ${treated_bams}
        samtools index treated_pool.bam
        treated_input=treated_pool.bam
    fi
    if [ \${#control_list[@]} -gt 1 ]; then
        samtools merge -f control_pool.bam ${control_bams}
        samtools index control_pool.bam
        control_input=control_pool.bam
    fi

    java -Xmx${java_heap} -jar ${omnipeaks_jar} compare \
        -t ${treated_input} \
        -c ${control_input} \
        --cs ${chrom_sizes} \
        --gap ${gap} \
        ${bin_arg} \
        ${fdr_arg} \
        -p ${prefix}

    diff_file=\$(ls ${prefix}*.peak 2>/dev/null | head -n 1 || true)
    if [ -z "\$diff_file" ]; then
        diff_file=\$(ls ${prefix}*.bed 2>/dev/null | head -n 1 || true)
    fi
    if [ -z "\$diff_file" ]; then
        echo "SPAN compare produced no output" >&2
        exit 1
    fi

    cp "\$diff_file" span.differential.tsv
    awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3}' "\$diff_file" > span.differential.bed

    # Attempt to split by sign using column 5 if numeric
    awk 'BEGIN{OFS="\t"} {if (\$5 ~ /^-?[0-9.]+\$/) {if (\$5>0) print \$1,\$2,\$3}}' "\$diff_file" > span.significant_up.bed
    awk 'BEGIN{OFS="\t"} {if (\$5 ~ /^-?[0-9.]+\$/) {if (\$5<0) print \$1,\$2,\$3}}' "\$diff_file" > span.significant_down.bed

    total=\$(grep -v '^#' "\$diff_file" | wc -l | awk '{print \$1}')
    n_up=\$(wc -l < span.significant_up.bed | awk '{print \$1}')
    n_down=\$(wc -l < span.significant_down.bed | awk '{print \$1}')
    n_fdr_pass=\$((n_up + n_down))
    cat <<-END_SUMMARY > span.summary.tsv
    method\tgroup\tcaller\ttreated_condition\tcontrol_condition\tn_tested\tn_fdr_pass\tn_up\tn_down
    span\t${meta.group}\tNA\t${meta.treated ?: 'NA'}\t${meta.control ?: 'NA'}\t\${total}\t\${n_fdr_pass}\t\${n_up}\t\${n_down}
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/\"//g')
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
