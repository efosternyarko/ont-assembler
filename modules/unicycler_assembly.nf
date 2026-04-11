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

    script:
    """
    unicycler \\
        --long ${reads} \\
        --out ${sample_id}_unicycler \\
        --threads ${task.cpus}

    cp ${sample_id}_unicycler/assembly.fasta ${sample_id}.fasta

    cp ${sample_id}_unicycler/unicycler.log ${sample_id}_unicycler.log 2>/dev/null \\
        || echo "sample: ${sample_id}" > ${sample_id}_unicycler.log
    """
}
