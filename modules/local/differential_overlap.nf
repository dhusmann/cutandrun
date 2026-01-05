process DIFFERENTIAL_OVERLAP {
    tag "overlap"
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/04_cross_comparison" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bedtools=2.31.1 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path callers_tsv
    path methods_tsv

    output:
    path "overlap_callers.tsv" , emit: callers
    path "overlap_methods.tsv" , emit: methods
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    set -euo pipefail

    printf "group\tcaller_a\tcaller_b\tintersection\tunion\tjaccard\n" > overlap_callers.tsv
    if [[ -s "${callers_tsv}" ]]; then
        awk -F'\t' 'BEGIN{OFS="\t"} { g=\$1; c=\$2; p=\$3; n[g]++; callers[g,n[g]]=c; paths[g,n[g]]=p }
            END { for (g in n) { for (i=1; i<=n[g]; i++) { for (j=i+1; j<=n[g]; j++) { print g, callers[g,i], paths[g,i], callers[g,j], paths[g,j] } } } }' "${callers_tsv}" \
            | while IFS=\$'\\t' read -r group caller_a path_a caller_b path_b; do
                jaccard_line=\$(bedtools jaccard -a "\$path_a" -b "\$path_b" | awk 'NR==2 {print \$2"\\t"\$3"\\t"\$4}')
                if [[ -z "\$jaccard_line" ]]; then
                    inter=0
                    union=0
                    jacc=0
                else
                    IFS=\$'\\t' read -r inter union jacc <<< "\$jaccard_line"
                fi
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "\$group" "\$caller_a" "\$caller_b" "\$inter" "\$union" "\$jacc" >> overlap_callers.tsv
            done
    fi

    printf "group\tmethod_a\tmethod_b\tcaller\tintersection\tunion\tjaccard\n" > overlap_methods.tsv
    if [[ -s "${methods_tsv}" ]]; then
        awk -F'\t' 'BEGIN{OFS="\t"} { g=\$1; m=\$2; c=\$3; p=\$4; n[g]++; methods[g,n[g]]=m; callers[g,n[g]]=c; paths[g,n[g]]=p }
            END { for (g in n) { for (i=1; i<=n[g]; i++) { for (j=i+1; j<=n[g]; j++) { print g, methods[g,i], callers[g,i], paths[g,i], methods[g,j], callers[g,j], paths[g,j] } } } }' "${methods_tsv}" \
            | while IFS=\$'\\t' read -r group method_a caller_a path_a method_b caller_b path_b; do
                jaccard_line=\$(bedtools jaccard -a "\$path_a" -b "\$path_b" | awk 'NR==2 {print \$2"\\t"\$3"\\t"\$4}')
                if [[ -z "\$jaccard_line" ]]; then
                    inter=0
                    union=0
                    jacc=0
                else
                    IFS=\$'\\t' read -r inter union jacc <<< "\$jaccard_line"
                fi
                caller="NA"
                if [[ "\$method_a" == "diffbind" ]]; then
                    caller="\$caller_a"
                elif [[ "\$method_b" == "diffbind" ]]; then
                    caller="\$caller_b"
                fi
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "\$group" "\$method_a" "\$method_b" "\$caller" "\$inter" "\$union" "\$jacc" >> overlap_methods.tsv
            done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """.stripIndent()
}
