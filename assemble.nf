#!/usr/bin/env nextflow
// ── assemble.nf — standalone ONT long-read assembly pipeline ──────────────
//
// Workflow:
//   1. Read QC   — seqkit stats → per-sample TSV → 8-panel metrics plot
//   2. Depth filter — drop samples below --min_read_depth
//   3. Assembly  — hybracter (default) | flye | raven
//
// Input (one of):
//   --samplesheet  CSV with columns:  id,reads
//   --input_dir    Directory of FASTQ / FASTQ.gz files (sample ID = filename stem)
//
// Key parameters:
//   --assembler          hybracter (default) | flye | raven
//   --genome_size        Expected genome size for depth calculation (default: 5m)
//   --min_read_depth     Drop samples below this estimated depth × (default: 20)
//   --chromosome_size    Hybracter: min contig length to call a chromosome (bp) (default: 2500000)
//                        Ignored when --hybracter_auto true
//   --hybracter_auto     Let hybracter estimate chromosome size (default: true)
//   --hybracter_no_medaka  Skip medaka polishing — required on macOS ARM (default: false)
//   --outdir             Output directory (default: assembly_results_<assembler>)
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
include { RAVEN_ASSEMBLY      } from './modules/raven_assembly'

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

// ── Help ──────────────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║                         ont-assembler                                   ║
    ║   ONT long-read assembly with read QC and depth filtering               ║
    ╚══════════════════════════════════════════════════════════════════════════╝

    USAGE
      nextflow run assemble.nf -c assemble.config -profile <profile> [options]

    INPUT (one required)
      --input_dir   <dir>   Directory of FASTQ / FASTQ.gz files
                            (sample ID = filename stem)
      --samplesheet <csv>   CSV with columns: id,reads

    COMMON OPTIONS
      --outdir              <dir>    Output directory (default: assembly_results_<assembler>)
      --assembler           <name>   hybracter (default) | flye | raven
      --genome_size         <size>   Expected genome size (default: 5m)
      --min_read_depth      <n>      Min estimated depth to assemble (default: 20)
      --hybracter_no_medaka          Skip medaka polishing (required on macOS ARM)
      --max_cpus            <n>      Max CPUs per process (default: 16)

    EXAMPLES
      # Hybracter (default) — Linux / HPC
      nextflow run assemble.nf -c assemble.config -profile conda \\
          --input_dir /path/to/fastq/

      # macOS Apple Silicon (M1 and above)
      CONDA_SUBDIR=osx-64 nextflow run assemble.nf -c assemble.config -profile conda,arm64 \\
          --input_dir /path/to/fastq/ \\
          --hybracter_no_medaka true

      # HPC (SLURM + Singularity)
      nextflow run assemble.nf -c assemble.config -profile singularity,slurm \\
          --input_dir /path/to/fastq/

      # Alternative assembler
      nextflow run assemble.nf -c assemble.config -profile conda \\
          --input_dir /path/to/fastq/ \\
          --assembler flye

      # Resume after interruption
      nextflow run assemble.nf -c assemble.config -profile conda \\
          --input_dir /path/to/fastq/ -resume

    NOTE  Assembled FASTAs (hybracter) are written to:
            assembly_results/all_assemblies/
          Pass this directory to enteric-typer with --input_dir.

    PROFILES
      conda           Local — conda environments
      mamba           Local — mamba (faster env solving)
      arm64           Add on Apple Silicon (combine with conda/mamba)
      docker          Local — Docker
      singularity     HPC — Singularity/Apptainer
      slurm           HPC — SLURM executor
      pbs             HPC — PBS/Torque executor

    Full documentation: README.md
    """.stripIndent()
    exit 0
}

// Resolve output directory before any process runs so publishDir closures in
// assemble.config pick up the correct value.
params.outdir = params.outdir ?: "assembly_results_${params.assembler}"

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
        ch_timing     = HYBRACTER_ASSEMBLY.out.timing
    } else if (asm == 'flye') {
        FLYE_ASSEMBLY(ch_reads_with_depth, params.genome_size)
        ch_assemblies = FLYE_ASSEMBLY.out.assembly
        ch_timing     = FLYE_ASSEMBLY.out.timing
    } else if (asm == 'raven') {
        RAVEN_ASSEMBLY(ch_reads_with_depth)
        ch_assemblies = RAVEN_ASSEMBLY.out.assembly
        ch_timing     = RAVEN_ASSEMBLY.out.timing
    } else {
        error "Unknown assembler '${params.assembler}'. Choose: hybracter | flye | raven"
    }

    // Emit assembled FASTAs (publishDir in assemble.config copies them to outdir/assemblies/)
    ch_assemblies.view { id, fasta -> "Assembled: ${id}  →  ${fasta}" }

    // Collate per-sample timing TSVs into a single summary (header from first file, data rows only from rest)
    ch_timing
        .collectFile(
            name:     "assembly_timing_summary.tsv",
            storeDir: "${params.outdir}",
            keepHeader: true,
            sort:     true
        )
}
