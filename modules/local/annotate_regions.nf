process ANNOTATE_REGIONS {
    label 'process_single'

    conda "bioconda::bedtools=2.31.0 conda-forge::python=3.8.3"
    container "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_1"

    input:
    tuple val(method), val(group), val(caller), path(regions), path(gene_bed), val(output_name)

    output:
    tuple val(method), val(group), val(caller), path("${output_name}"), emit: annotated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    annotate_regions.py \
        --regions ${regions} \
        --gene-bed ${gene_bed} \
        --out ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
