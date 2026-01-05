process DIFFBIND_RUN {
    label 'process_diffbind'

    conda "conda-forge::r-base=4.2.3 bioconda::bioconductor-diffbind conda-forge::r-jsonlite conda-forge::r-yaml"
    container "quay.io/biocontainers/bioconductor-diffbind:3.14.0--r42_0"

    input:
    tuple val(group), val(caller), path(records)
    val contrast
    val use_spikein
    val fdr
    val lfc
    val min_overlap
    val backend
    val recenter
    val summits
    val norm_method
    val extra_params
    val export_sheets
    val allow_partial

    output:
    tuple val(group), val(caller), path("diffbind.results.tsv", optional: true), emit: results
    tuple val(group), val(caller), path("diffbind.significant.bed", optional: true), emit: bed
    tuple val(group), val(caller), path("diffbind.significant_up.bed", optional: true), emit: bed_up
    tuple val(group), val(caller), path("diffbind.significant_down.bed", optional: true), emit: bed_down
    tuple val(group), val(caller), path("diffbind.summary.tsv"), emit: summary
    tuple val(group), val(caller), path("diffbind.samplesheet.csv"), emit: samplesheet
    tuple val(group), val(caller), path("diffbind.normalization_factors.tsv", optional: true), emit: norm_factors_out
    tuple val(group), val(caller), path("diffbind.dba.rds", optional: true), emit: dba
    tuple val(group), val(caller), path("diffbind.error.txt", optional: true), emit: error
    tuple val(group), val(caller), path("plots", optional: true), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_arg = extra_params ? "--extra_params ${extra_params}" : ''
    """
    diffbind_run.R \
        --records ${records} \
        --outdir . \
        --contrast '${contrast}' \
        --caller '${caller}' \
        --group '${group}' \
        --use_spikein ${use_spikein} \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --min_overlap ${min_overlap} \
        --backend ${backend} \
        --recenter ${recenter} \
        --summits ${summits} \
        --norm_method ${norm_method} \
        ${extra_arg} \
        --export_sheets ${export_sheets} \
        --allow_partial ${allow_partial}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p plots
    cat <<-EOF > diffbind.samplesheet.csv
SampleID,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller,Tissue
stub,${group},stub,1,stub.bam,stub.peak,${caller},CUTRUN
EOF

    cat <<-EOF > diffbind.results.tsv
chr\tstart\tend\tlog2FC\tpval\tFDR
chr1\t1\t2\t0.0\t1.0\t1.0
EOF

    printf "chr1\t1\t2\n" > diffbind.significant.bed
    printf "chr1\t1\t2\n" > diffbind.significant_up.bed
    printf "chr1\t1\t2\n" > diffbind.significant_down.bed

    cat <<-EOF > diffbind.summary.tsv
caller\tgroup\ttreated\tcontrol\tn_tested\tn_fdr_pass\tn_up\tn_down\tstatus\treason
${caller}\t${group}\tstub\tstub\t1\t0\t0\t0\tRUN\tok
EOF

    cat <<-EOF > diffbind.normalization_factors.tsv
sample_id\tsize_factor
stub\t1.0
EOF

    echo "stub" > diffbind.dba.rds
    echo "stub" > plots/PCA.pdf
    echo "stub" > plots/correlation_heatmap.pdf
    echo "stub" > plots/MA.pdf
    echo "stub" > plots/volcano.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
