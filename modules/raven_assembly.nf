// ── RAVEN_ASSEMBLY: ONT long-read assembly with Raven ──────────────────────
// Raven is a de novo assembler for long uncorrected reads.
// Output is FASTA written to stdout — redirected to file here.

process RAVEN_ASSEMBLY {
    tag "${sample_id}"
    label 'high'
    errorStrategy 'ignore'   // failed assembly skips the sample; never kills the whole run

    conda     "${projectDir}/envs/raven.yml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: assembly
    path("${sample_id}_raven.gfa"),                   emit: gfa
    path("${sample_id}_timing.tsv"),                  emit: timing

    script:
    """
    _start=\$(date +%s)

    raven \\
        --threads ${task.cpus} \\
        --graphical-fragment-assembly ${sample_id}_raven.gfa \\
        ${reads} \\
        > ${sample_id}.fasta

    _end=\$(date +%s)
    _elapsed=\$(( _end - _start ))
    _hms=\$(printf '%02d:%02d:%02d' \$(( _elapsed/3600 )) \$(( (_elapsed%3600)/60 )) \$(( _elapsed%60 )))
    printf 'sample\tassembler\tduration_seconds\tduration_hms\n' > ${sample_id}_timing.tsv
    printf '%s\traven\t%d\t%s\n' "${sample_id}" "\$_elapsed" "\$_hms" >> ${sample_id}_timing.tsv
    """
}
