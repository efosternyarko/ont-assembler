#!/usr/bin/env nextflow
// ── assemble.nf — standalone ONT long-read assembly pipeline ──────────────
//
// Workflow:
//   1. Read QC   — seqkit stats → per-sample TSV → 8-panel metrics plot
//   2. Depth filter — drop samples below --min_read_depth
//   3. Assembly  — hybracter (default) | flye | dragonflye | unicycler
//
// Input (one of):
//   --samplesheet  CSV with columns:  id,reads
//   --input_dir    Directory of FASTQ / FASTQ.gz files (sample ID = filename stem)
//
// Key parameters:
//   --assembler          hybracter (default) | flye | dragonflye | unicycler
//   --genome_size        Expected genome size for depth calculation (default: 5m)
//   --min_read_depth     Drop samples below this estimated depth × (default: 20)
//   --chromosome_size    Hybracter: min contig length to call a chromosome (bp) (default: 2500000)
//                        Ignored when --hybracter_auto true
//   --hybracter_auto     Let hybracter estimate chromosome size (default: true)
//   --hybracter_no_medaka  Skip medaka polishing — required on macOS ARM (default: false)
//   --outdir             Output directory (default: assembly_results)
//
// Usage examples:
//   # Hybracter (default), Linux/HPC:
//   nextflow run assemble.nf -profile conda \
//       --input_dir /path/to/fastq \
//       --genome_size 5m
//
//   # macOS Apple Silicon:
//   CONDA_SUBDIR=osx-64 nextflow run assemble.nf -profile conda,arm64 \
//       --input_dir /path/to/fastq \
//       --hybracter_no_medaka true
//
//   # Flye assembler:
//   nextflow run assemble.nf -profile conda \
//       --samplesheet samples.csv \
//       --assembler flye

nextflow.enable.dsl = 2

include { ONT_READ_QC       } from './modules/ont_read_qc'
include { ONT_PLOT_METRICS  } from './modules/ont_plot_metrics'
include { HYBRACTER_ASSEMBLY  } from './modules/hybracter_assembly'
include { FLYE_ASSEMBLY       } from './modules/flye_assembly'
include { DRAGONFLYE_ASSEMBLY } from './modules/dragonflye_assembly'
include { UNICYCLER_ASSEMBLY  } from './modules/unicycler_assembly'

// ── Input channel helpers ─────────────────────────────────────────────────

def samplesheet_channel(csv_path) {
    Channel.fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def id    = row.id    ?: row.sample_id ?: row.sample
            def reads = file(row.reads, checkIfExists: true)
            tuple(id, reads)
        }
}

def input_dir_channel(dir_path) {
    Channel.fromPath("${dir_path}/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
        .map { f ->
            // Strip .fastq / .fastq.gz / .fq / .fq.gz suffix to get sample ID
            def id = f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(id, f)
        }
}

// ── Workflow ──────────────────────────────────────────────────────────────

workflow {

    // 1. Build input channel
    if (!params.samplesheet && !params.input_dir) {
        error "Provide --samplesheet <csv> or --input_dir <dir>"
    }
    ch_reads = params.samplesheet
        ? samplesheet_channel(params.samplesheet)
        : input_dir_channel(params.input_dir)

    // 2. Read QC — per-sample stats
    ONT_READ_QC(ch_reads, params.genome_size)

    // 3. Depth filter — parse estimated_depth from stats TSV, drop below threshold
    //    TSV columns: sample_id  num_reads  total_bases  read_N50  mean_length
    //                 mean_quality  gc_pct  estimated_depth   (index 7)
    ch_passing = ONT_READ_QC.out.stats
        .map { id, tsv ->
            def lines = tsv.text.readLines()
            def depth = lines.size() > 1 ? lines[1].split('\t')[7].toDouble() : 0.0
            tuple(id, depth, tsv)
        }
        .filter { id, depth, tsv ->
            if (depth < params.min_read_depth as double) {
                log.warn "Skipping ${id}: estimated depth ${depth}× < ${params.min_read_depth}×"
                return false
            }
            return true
        }

    // 4. Metrics plot — passing samples only (depth ≥ min_read_depth)
    ONT_PLOT_METRICS(
        ch_passing.map { _id, _depth, tsv -> tsv }.collect(),
        params.min_read_depth
    )

    // Re-join passing IDs with original reads channel to get the reads path
    ch_reads_with_depth = ch_passing
        .map { id, _depth, _tsv -> tuple(id, 'pass') }
        .join(ch_reads, by: 0)
        .map { id, _flag, reads -> tuple(id, reads) }

    // 5. Assembly — dispatch to selected assembler
    def asm = params.assembler.toLowerCase()

    if (asm == 'hybracter') {
        HYBRACTER_ASSEMBLY(
            ch_reads_with_depth,
            params.chromosome_size,
            params.hybracter_auto,
            params.hybracter_no_medaka
        )
        ch_assemblies = HYBRACTER_ASSEMBLY.out.assembly
    } else if (asm == 'flye') {
        FLYE_ASSEMBLY(ch_reads_with_depth, params.genome_size)
        ch_assemblies = FLYE_ASSEMBLY.out.assembly
    } else if (asm == 'dragonflye') {
        DRAGONFLYE_ASSEMBLY(ch_reads_with_depth, params.genome_size)
        ch_assemblies = DRAGONFLYE_ASSEMBLY.out.assembly
    } else if (asm == 'unicycler') {
        UNICYCLER_ASSEMBLY(ch_reads_with_depth)
        ch_assemblies = UNICYCLER_ASSEMBLY.out.assembly
    } else {
        error "Unknown assembler '${params.assembler}'. Choose: hybracter | flye | dragonflye | unicycler"
    }

    // Emit assembled FASTAs (publishDir in assemble.config copies them to outdir/assemblies/)
    ch_assemblies.view { id, fasta -> "Assembled: ${id}  →  ${fasta}" }
}
