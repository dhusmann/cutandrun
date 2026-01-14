process EXPORT_DIFFBIND_SAMPLESHEETS {
    label 'process_single'

    input:
    tuple val(group), val(caller), path(samplesheet)

    output:
    tuple val(group), val(caller), path("${group}.csv"), emit: sheet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp ${samplesheet} ${group}.csv
    """
}
