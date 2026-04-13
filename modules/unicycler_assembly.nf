// ── UNICYCLER_ASSEMBLY: ONT long-read assembly with Unicycler ─────────────
// Runs Unicycler in long-read only mode (--long).
// For hybrid assembly (Illumina + ONT) use the standard short+long flags instead.

process UNICYCLER_ASSEMBLY {
    tag "${sample_id}"
    label 'high'
    errorStrategy 'ignore'   // failed assembly skips the sample; never kills the whole run

    conda     "${projectDir}/envs/unicycler.yml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: assembly
    path("${sample_id}_unicycler.log"),                emit: log
    path("${sample_id}_timing.tsv"),                   emit: timing

    script:
    """
    _start=\$(date +%s)

    unicycler \\
        --long ${reads} \\
        --out ${sample_id}_unicycler \\
        --threads ${task.cpus}

    cp ${sample_id}_unicycler/assembly.fasta ${sample_id}.fasta

    cp ${sample_id}_unicycler/unicycler.log ${sample_id}_unicycler.log 2>/dev/null \\
        || echo "sample: ${sample_id}" > ${sample_id}_unicycler.log

    _end=\$(date +%s)
    _elapsed=\$(( _end - _start ))
    _hms=\$(printf '%02d:%02d:%02d' \$(( _elapsed/3600 )) \$(( (_elapsed%3600)/60 )) \$(( _elapsed%60 )))
    printf 'sample\tassembler\tduration_seconds\tduration_hms\n' > ${sample_id}_timing.tsv
    printf '%s\tunicycler\t%d\t%s\n' "${sample_id}" "\$_elapsed" "\$_hms" >> ${sample_id}_timing.tsv
    """
}
