// ── DRAGONFLYE_ASSEMBLY: ONT long-read assembly with DragonFlye ───────────
// DragonFlye wraps Flye with additional QC steps (read trimming, length
// filtering, optional short-read polishing).

process DRAGONFLYE_ASSEMBLY {
    tag "${sample_id}"
    label 'high'
    errorStrategy 'ignore'   // failed assembly skips the sample; never kills the whole run

    conda     "${projectDir}/envs/dragonflye.yml"
    container 'quay.io/biocontainers/dragonflye:1.2.0--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads)
    val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: assembly
    path("${sample_id}_dragonflye.log"),               emit: log
    path("${sample_id}_timing.tsv"),                   emit: timing

    script:
    """
    _start=\$(date +%s)

    dragonflye \\
        --reads ${reads} \\
        --outdir ${sample_id}_dragonflye \\
        --gsize ${genome_size} \\
        --cpus ${task.cpus} \\
        --minreadlen 1000 \\
        --minquality 8

    cp ${sample_id}_dragonflye/contigs.fa ${sample_id}.fasta

    cp ${sample_id}_dragonflye/dragonflye.log ${sample_id}_dragonflye.log 2>/dev/null \\
        || echo "sample: ${sample_id}" > ${sample_id}_dragonflye.log

    _end=\$(date +%s)
    _elapsed=\$(( _end - _start ))
    _hms=\$(printf '%02d:%02d:%02d' \$(( _elapsed/3600 )) \$(( (_elapsed%3600)/60 )) \$(( _elapsed%60 )))
    printf 'sample\tassembler\tduration_seconds\tduration_hms\n' > ${sample_id}_timing.tsv
    printf '%s\tdragonflye\t%d\t%s\n' "${sample_id}" "\$_elapsed" "\$_hms" >> ${sample_id}_timing.tsv
    """
}
