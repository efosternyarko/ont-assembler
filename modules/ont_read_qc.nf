// ── ONT_READ_QC: per-sample ONT read statistics ────────────────────────────
// Computes read length and quality metrics from ONT FASTQ using seqkit stats.
// Estimated depth = total bases / genome_size_estimate.
//
// Output TSV columns:
//   sample_id  num_reads  total_bases  read_N50  mean_length
//   mean_quality  gc_pct  estimated_depth

process ONT_READ_QC {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/seqkit.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(sample_id), path(reads)
    val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}_read_stats.tsv"), emit: stats

    script:
    """
    seqkit stats -a -T -j ${task.cpus} ${reads} > seqkit_raw.txt
    parse_read_stats.py ${sample_id} ${genome_size} seqkit_raw.txt \\
        > ${sample_id}_read_stats.tsv
    """
}
