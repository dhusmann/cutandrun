process ANNOTATE_REGIONS {
    label 'process_single'

    conda "bioconda::bedtools=2.31.0"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"

    input:
    tuple val(method), val(group), val(caller), path(regions), path(gene_bed), val(output_name)

    output:
    tuple val(method), val(group), val(caller), path("${output_name}"), emit: annotated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    out_dir=\$(dirname "${output_name}")
    mkdir -p "\$out_dir"
    first_line=\$(head -n 1 ${regions} || true)
    if [[ -z "\$first_line" ]]; then
        : > ${output_name}
    else
        if echo "\$first_line" | awk -F'\t' 'NF<3{exit 0} \$2 ~ /^[0-9]+${'$'}/ {exit 1} {exit 0}'; then
            header_line="\$first_line"
            tail -n +2 ${regions} > regions.data
        else
            header_line=""
            cp ${regions} regions.data
        fi

        if [[ ! -s regions.data ]]; then
            if [[ -n "\$header_line" ]]; then
                printf "%s\tnearest_feature_id\tnearest_gene_name\tdistance_to_feature\n" "\$header_line" > ${output_name}
            else
                : > ${output_name}
            fi
        else
            bcols=\$(awk -F'\t' 'NF>0{print NF; exit}' ${gene_bed} || true)
            bcols=\${bcols:-0}
            awk -F'\t' -v OFS='\t' '{print \$0, NR-1}' regions.data > regions.with_idx
            bedtools closest -d -a regions.with_idx -b ${gene_bed} > regions.closest

            awk -F'\t' -v OFS='\t' -v bcols="\$bcols" -v header="\$header_line" 'BEGIN{ if (header != "") print header, "nearest_feature_id", "nearest_gene_name", "distance_to_feature" } {
                a_cols = NF - bcols - 1
                if (a_cols < 1) next
                out = \$1
                for (i=2; i<=a_cols-1; i++) out = out OFS \$i
                gene_id = "NA"
                gene_name = "NA"
                if (bcols >= 4) gene_id = \$(a_cols + 4)
                if (bcols >= 5) gene_name = \$(a_cols + 5)
                if (gene_id == "." || gene_id == "") gene_id = "NA"
                if (gene_name == "." || gene_name == "") gene_name = "NA"
                distance = \$(NF)
                print out, gene_id, gene_name, distance
            }' regions.closest > ${output_name}
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
