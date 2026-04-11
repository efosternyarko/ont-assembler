// ── FLYE_ASSEMBLY: ONT long-read assembly with Flye ───────────────────────
// Uses --nano-hq mode (Guppy 5+ / Dorado Q20+ reads — the current standard).
// Override to --nano-raw for older/lower-accuracy basecalling if needed.

process FLYE_ASSEMBLY {
    tag "${sample_id}"
    label 'high'
    errorStrategy 'ignore'   // failed assembly skips the sample; never kills the whole run

    conda     "${projectDir}/envs/flye.yml"
    container 'quay.io/biocontainers/flye:2.9.4--py39h6935b12_1'

    input:
    tuple val(sample_id), path(reads)
    val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: assembly
    path("${sample_id}_flye_info.txt"),                        emit: info

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --out-dir ${sample_id}_flye \\
        --genome-size ${genome_size} \\
        --threads ${task.cpus}

    cp ${sample_id}_flye/assembly.fasta ${sample_id}.fasta

    # Save basic assembly info
    cp ${sample_id}_flye/assembly_info.txt ${sample_id}_flye_info.txt 2>/dev/null \\
        || echo "sample: ${sample_id}" > ${sample_id}_flye_info.txt
    """
}
