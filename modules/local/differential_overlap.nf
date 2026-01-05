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

    declare -A bp_cache

    calc_bp() {
        local bed="\$1"
        if [[ ! -s "\$bed" ]]; then
            echo 0
            return
        fi
        bedtools sort -i "\$bed" | bedtools merge -i - | awk '{sum += (\$3 - \$2)} END {print sum + 0}'
    }

    get_bp() {
        local bed="\$1"
        if [[ -n "\${bp_cache[\$bed]+x}" ]]; then
            echo "\${bp_cache[\$bed]}"
            return
        fi
        local value
        value=\$(calc_bp "\$bed")
        bp_cache["\$bed"]="\$value"
        echo "\$value"
    }

    printf "group\tcaller_a\tcaller_b\tintersection\tunion\tjaccard\n" > overlap_callers.tsv
    if [[ -s "${callers_tsv}" ]]; then
        awk -F'\\t' 'BEGIN{OFS="\\t"} { g=\$1; c=\$2; p=\$3; n[g]++; callers[g,n[g]]=c; paths[g,n[g]]=p }
            END { for (g in n) { for (i=1; i<=n[g]; i++) { for (j=i+1; j<=n[g]; j++) { print g, callers[g,i], paths[g,i], callers[g,j], paths[g,j] } } } }' "${callers_tsv}" \
            | while IFS=\$'\\t' read -r group caller_a path_a caller_b path_b; do
                jaccard_line=\$(bedtools jaccard -a "\$path_a" -b "\$path_b" | awk 'NR==2 {print \$1"\\t"\$2"\\t"\$3}')
                if [[ -z "\$jaccard_line" ]]; then
                    inter=0
                    union=0
                    jacc=0
                else
                    IFS=\$'\\t' read -r inter union jacc <<< "\$jaccard_line"
                fi
                printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\$group" "\$caller_a" "\$caller_b" "\$inter" "\$union" "\$jacc" >> overlap_callers.tsv
            done
    fi

    printf "group\tcaller\tmethod_a\tmethod_b\tn_a\tn_b\tn_overlap\tjaccard\tspan_mode_used\n" > overlap_methods.tsv
    if [[ -s "${methods_tsv}" ]]; then
        awk -F'\\t' 'BEGIN{OFS="\\t"} { g=\$1; m=\$2; c=\$3; p=\$4; s=\$5; n[g]++; methods[g,n[g]]=m; callers[g,n[g]]=c; paths[g,n[g]]=p; spans[g,n[g]]=s }
            END { for (g in n) { for (i=1; i<=n[g]; i++) { for (j=i+1; j<=n[g]; j++) { if (methods[g,i] == methods[g,j]) { next } print g, methods[g,i], callers[g,i], paths[g,i], spans[g,i], methods[g,j], callers[g,j], paths[g,j], spans[g,j] } } } }' "${methods_tsv}" \
            | while IFS=\$'\\t' read -r group method_a caller_a path_a span_a method_b caller_b path_b span_b; do
                n_a=\$(get_bp "\$path_a")
                n_b=\$(get_bp "\$path_b")
                jaccard_line=\$(bedtools jaccard -a "\$path_a" -b "\$path_b" | awk 'NR==2 {print \$1"\\t"\$3}')
                if [[ -z "\$jaccard_line" ]]; then
                    n_overlap=0
                    jacc=0
                else
                    IFS=\$'\\t' read -r n_overlap jacc <<< "\$jaccard_line"
                fi
                caller="NA"
                if [[ "\$method_a" == "diffbind" ]]; then
                    caller="\$caller_a"
                elif [[ "\$method_b" == "diffbind" ]]; then
                    caller="\$caller_b"
                fi
                span_mode="NA"
                if [[ "\$method_a" == "span" ]]; then
                    span_mode="\$span_a"
                elif [[ "\$method_b" == "span" ]]; then
                    span_mode="\$span_b"
                fi
                printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\$group" "\$caller" "\$method_a" "\$method_b" "\$n_a" "\$n_b" "\$n_overlap" "\$jacc" "\$span_mode" >> overlap_methods.tsv
            done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """.stripIndent()
}
