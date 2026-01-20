process SPAN_OMNIPEAKS_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::openjdk=21.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://eclipse-temurin:21-jdk' :
        'eclipse-temurin:21-jdk' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    path  chrom_sizes
    path  omnipeaks_jar
    val   gap
    val   fdr
    val   java_heap

    output:
    tuple val(meta), path("*.peak"), emit: peaks
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def raw_prefix = task.ext.prefix ?: "${meta.id}_${meta.caller}"
    def prefix = raw_prefix.endsWith('.peak') ? raw_prefix : "${raw_prefix}.peak"
    def fdr_arg = fdr ? "--fdr ${fdr}" : ''
    """
    java --add-modules jdk.incubator.vector -Xmx${java_heap} -jar ${omnipeaks_jar} analyze \
        -t ${treatment_bam} \
        -c ${control_bam} \
        --cs ${chrom_sizes} \
        --gap ${gap} \
        ${fdr_arg} \
        -p ${prefix} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/"//g')
    END_VERSIONS
    """
}
