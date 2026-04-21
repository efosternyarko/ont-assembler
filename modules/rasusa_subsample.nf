// ── RASUSA_SUBSAMPLE: subsample ONT reads to a target depth ────────────────
// Runs only on samples whose estimated depth exceeds params.max_depth.
// Subsampling is reproducible via a fixed seed (42).

process RASUSA_SUBSAMPLE {
    tag "${sample_id} (depth > ${max_depth}×)"
    label 'medium'

    conda "${projectDir}/envs/rasusa.yml"

    input:
    tuple val(sample_id), path(reads)
    val(genome_size)
    val(max_depth)

    output:
    tuple val(sample_id), path("${sample_id}_subsampled.fastq.gz"), emit: reads

    script:
    """
    rasusa reads \\
        --input   ${reads} \\
        --coverage ${max_depth} \\
        --genome-size ${genome_size} \\
        --output  ${sample_id}_subsampled.fastq.gz \\
        --seed    42
    """
}
