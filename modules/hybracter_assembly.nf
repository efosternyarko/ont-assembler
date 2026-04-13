// ── HYBRACTER_ASSEMBLY: ONT long-read assembly with Hybracter ─────────────
// Runs hybracter long (ONT-only mode).
// Complete assemblies: chromosome + plasmid FASTAs are concatenated.
// Incomplete assemblies: the draft contig FASTA is used as-is.
//
// Two modes (controlled by params.hybracter_auto):
//   false (default) — 3-column CSV: sample, reads, min_chromosome_length
//   true            — 2-column CSV + --auto flag: hybracter estimates chromosome length

process HYBRACTER_ASSEMBLY {
    tag "${sample_id}"
    label 'high'
    errorStrategy 'ignore'   // failed assembly skips the sample; never kills the whole run

    conda     "${projectDir}/envs/hybracter.yml"

    input:
    tuple val(sample_id), path(reads)
    val(chromosome_size)
    val(auto_mode)
    val(no_medaka)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"),          emit: assembly
    path("${sample_id}_hybracter_summary.tsv"),                 emit: summary
    path("${sample_id}_chromosome.fasta"),  optional: true,     emit: chromosome
    path("${sample_id}_plasmid.fasta"),     optional: true,     emit: plasmid
    path("${sample_id}_incomplete.fasta"),  optional: true,     emit: incomplete
    path("${sample_id}_flye.gfa"),          optional: true,     emit: gfa
    path("${sample_id}_timing.tsv"),                            emit: timing

    script:
    def auto_flag     = auto_mode ? '--auto' : ''
    def no_medaka_flag = no_medaka ? '--no_medaka' : ''
    """
    # Hybracter 0.12+ requires a CSV input (no header).
    # --auto mode: 2 columns (sample, reads) — hybracter estimates chromosome length.
    # Manual mode: 3 columns (sample, reads, min_chromosome_length).
    if ${auto_mode}; then
        echo "${sample_id},${reads}" > hybracter_input.csv
    else
        echo "${sample_id},${reads},${chromosome_size}" > hybracter_input.csv
    fi

    _start=\$(date +%s)

    # Force osx-64 for internal Snakemake conda envs (needed on Apple Silicon;
    # has no effect on Linux/x86_64).  --conda-prefix reuses envs already built
    # during 'hybracter install' so arm64 builds don't need to solve again.
    export CONDA_SUBDIR=osx-64
    hybracter long \\
        -i hybracter_input.csv \\
        -o ${sample_id}_hybracter \\
        -t ${task.cpus} \\
        --conda-prefix ${projectDir}/.snakemake/conda \\
        ${no_medaka_flag} \\
        ${auto_flag}

    final="${sample_id}_hybracter/FINAL"
    chrom="\${final}/${sample_id}_chromosome.fasta"
    plasmid="\${final}/${sample_id}_plasmid.fasta"
    incomplete="\${final}/${sample_id}_incomplete.fasta"

    if [ -f "\${chrom}" ]; then
        cat "\${chrom}" > ${sample_id}.fasta
        # Append plasmid contigs if file is non-empty
        [ -s "\${plasmid}" ] && cat "\${plasmid}" >> ${sample_id}.fasta || true
    elif [ -f "\${incomplete}" ]; then
        cp "\${incomplete}" ${sample_id}.fasta
    else
        # Fallback: first FASTA found anywhere in the output tree
        fallback=\$(find "${sample_id}_hybracter" -name "*.fasta" | head -1)
        if [ -n "\${fallback}" ]; then
            cp "\${fallback}" ${sample_id}.fasta
        else
            echo "ERROR: hybracter produced no assembly for ${sample_id}" >&2; exit 1
        fi
    fi

    # Preserve individual chromosome / plasmid / incomplete FASTAs
    [ -f "\${chrom}" ]      && cp "\${chrom}"      ${sample_id}_chromosome.fasta  || true
    [ -s "\${plasmid}" ]    && cp "\${plasmid}"    ${sample_id}_plasmid.fasta     || true
    [ -f "\${incomplete}" ] && cp "\${incomplete}" ${sample_id}_incomplete.fasta  || true

    # Copy hybracter summary table (one row per sample)
    summary_src="${sample_id}_hybracter/hybracter_summary.tsv"
    if [ -f "\${summary_src}" ]; then
        cp "\${summary_src}" ${sample_id}_hybracter_summary.tsv
    else
        printf 'sample\\tstatus\\n${sample_id}\\tunknown\\n' > ${sample_id}_hybracter_summary.tsv
    fi

    # Copy Flye assembly graph — buried at processing/assemblies/{sample}/assembly_graph.gfa
    cp "${sample_id}_hybracter/processing/assemblies/${sample_id}/assembly_graph.gfa" \
        "${sample_id}_flye.gfa" 2>/dev/null || true

    _end=\$(date +%s)
    _elapsed=\$(( _end - _start ))
    _hms=\$(printf '%02d:%02d:%02d' \$(( _elapsed/3600 )) \$(( (_elapsed%3600)/60 )) \$(( _elapsed%60 )))
    printf 'sample\tassembler\tduration_seconds\tduration_hms\n' > ${sample_id}_timing.tsv
    printf '%s\thybracter\t%d\t%s\n' "${sample_id}" "\$_elapsed" "\$_hms" >> ${sample_id}_timing.tsv
    """
}
