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
    path("${sample_id}_flye_info.txt"),                emit: info
    path("${sample_id}_flye.gfa"),                     emit: gfa
    path("${sample_id}_timing.tsv"),                   emit: timing

    script:
    """
    _start=\$(date +%s)

    flye \\
        --nano-hq ${reads} \\
        --out-dir ${sample_id}_flye \\
        --genome-size ${genome_size} \\
        --threads ${task.cpus}

    cp ${sample_id}_flye/assembly.fasta ${sample_id}.fasta

    # Save assembly info and graph
    cp ${sample_id}_flye/assembly_info.txt ${sample_id}_flye_info.txt 2>/dev/null \\
        || echo "sample: ${sample_id}" > ${sample_id}_flye_info.txt
    cp ${sample_id}_flye/assembly_graph.gfa ${sample_id}_flye.gfa 2>/dev/null || true

    _end=\$(date +%s)
    _elapsed=\$(( _end - _start ))
    _hms=\$(printf '%02d:%02d:%02d' \$(( _elapsed/3600 )) \$(( (_elapsed%3600)/60 )) \$(( _elapsed%60 )))
    printf 'sample\tassembler\tduration_seconds\tduration_hms\n' > ${sample_id}_timing.tsv
    printf '%s\tflye\t%d\t%s\n' "${sample_id}" "\$_elapsed" "\$_hms" >> ${sample_id}_timing.tsv
    """
}
