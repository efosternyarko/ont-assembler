// ── ONT_PLOT_METRICS: read metrics visualisation across all samples ────────
// Merges per-sample read stats TSVs and generates an 8-panel figure:
//   Row 1: Read N50 histogram | Read N50 boxplot | Depth histogram | Depth boxplot
//   Row 2: Total reads hist   | Total reads box  | Mean length hist | Mean length box
//
// A vertical dashed line on the depth histogram marks --min_read_depth.

process ONT_PLOT_METRICS {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0'

    input:
    path(stats_files)   // collected list of *_read_stats.tsv files
    val(min_depth)

    output:
    path("ont_read_metrics.png"),          emit: plot
    path("ont_read_metrics_summary.tsv"),  emit: summary

    script:
    """
    plot_ont_metrics.py \\
        --stats ${stats_files} \\
        --min_depth ${min_depth} \\
        --output ont_read_metrics.png \\
        --summary ont_read_metrics_summary.tsv
    """
}
