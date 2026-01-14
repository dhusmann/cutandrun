process RENAME_SUMMARY {
    label 'process_single'

    input:
    tuple val(prefix), path(summary)

    output:
    path "${prefix}.summary.tsv", emit: summary

    script:
    """
    cp ${summary} ${prefix}.summary.tsv
    """
}
