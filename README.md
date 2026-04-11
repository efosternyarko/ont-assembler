# ont-assembler

A Nextflow pipeline for assembling Oxford Nanopore (ONT) long reads into genome assemblies.

**Workflow steps:**
1. **Read QC** — per-sample metrics (read N50, depth, total reads, mean length) via `seqkit stats`
2. **Depth filter** — samples below `--min_read_depth` (default 20×) are skipped with a warning
3. **Read metrics plot** — 8-panel figure (histograms + boxplots) for passing samples
4. **Assembly** — [Hybracter](https://github.com/gbouras13/hybracter) by default, with [Flye](https://github.com/fenderglass/Flye), [DragonFlye](https://github.com/rpetit3/dragonflye), and [Unicycler](https://github.com/rrwick/Unicycler) as alternatives

Assembled FASTAs are ready to feed directly into [enteric-typer](https://github.com/efosternyarko/enteric-typer).

---

## Installation

### Requirements

- [Nextflow](https://nextflow.io) ≥ 23.04
- [Conda](https://conda-forge.org/miniforge/) or Mamba (recommended)
- Java 17+ (`java -version`)

### Clone

```bash
git clone https://github.com/efosternyarko/ont-assembler
cd ont-assembler
```

### Hybracter: one-time setup

Hybracter downloads its internal databases on first use. Run this once before the first assembly:

```bash
# Linux / Intel Mac
conda run --prefix .snakemake/conda/hybracter-env hybracter install

# macOS Apple Silicon
CONDA_SUBDIR=osx-64 conda run --prefix .snakemake/conda/hybracter-env hybracter install
```

> Nextflow creates the conda environments on the first run. If `hybracter install`
> fails because the env does not exist yet, run the pipeline once with a small test
> file first (`errorStrategy = 'ignore'` prevents failure), then run `hybracter install`.

---

## Quick start

```bash
# Hybracter (default) — Linux / HPC
nextflow run assemble.nf -c assemble.config -profile conda \
    --input_dir /path/to/fastq/

# macOS Apple Silicon (M1/M2/M3/M4)
CONDA_SUBDIR=osx-64 nextflow run assemble.nf -c assemble.config -profile conda,arm64 \
    --input_dir /path/to/fastq/ \
    --hybracter_no_medaka true

# Samplesheet input (CSV with columns: id,reads)
nextflow run assemble.nf -c assemble.config -profile conda \
    --samplesheet samples.csv

# Alternative assembler
nextflow run assemble.nf -c assemble.config -profile conda \
    --input_dir /path/to/fastq/ \
    --assembler flye
```

**Feed assembled FASTAs into enteric-typer** (hybracter default):

```bash
nextflow run /path/to/enteric-typer/main.nf -profile conda \
    --input_dir assembly_results/all_assemblies/ \
    --outdir    typing_results/
```

---

## Assemblers

| `--assembler` | Tool | Notes |
|---|---|---|
| `hybracter` | [Hybracter](https://github.com/gbouras13/hybracter) | **Default.** Circularises chromosomes and plasmids. Uses `--auto` to estimate chromosome size automatically. |
| `flye` | [Flye](https://github.com/fenderglass/Flye) | `--nano-hq` mode (Guppy 5+ / Dorado Q20+ reads). |
| `dragonflye` | [DragonFlye](https://github.com/rpetit3/dragonflye) | Flye wrapper with read trimming and length filtering. |
| `unicycler` | [Unicycler](https://github.com/rrwick/Unicycler) | Long-read only mode (`--long`). Circularises assemblies where possible. |

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | `null` | Directory of FASTQ / FASTQ.gz files. Sample ID = filename stem (strip extension). |
| `--samplesheet` | `null` | CSV with columns `id,reads` |
| `--outdir` | `assembly_results` | Output directory |
| `--assembler` | `hybracter` | Assembly tool: `hybracter` \| `flye` \| `dragonflye` \| `unicycler` |
| `--genome_size` | `5m` | Expected genome size for depth calculation (e.g. `5m`, `4500000`) |
| `--min_read_depth` | `20` | Minimum estimated depth (×). Shallower samples are skipped and excluded from the plot. |
| `--chromosome_size` | `2500000` | Hybracter: minimum contig length (bp) to call a chromosome. Ignored when `--hybracter_auto true`. |
| `--hybracter_auto` | `true` | Let hybracter estimate chromosome size automatically (`--auto`). |
| `--hybracter_no_medaka` | `false` | Skip medaka polishing. **Required on macOS Apple Silicon** (OpenSSL conflict). Leave `false` on Linux/HPC. |
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |

---

## Output files

```
assembly_results/
│
│  ── Read QC ──────────────────────────────────────────────────────────────────
│
├── ont_read_qc/
│   └── {sample}_read_stats.tsv
│         Per-sample read metrics (tab-separated, one data row):
│           sample_id       — sample name
│           num_reads       — total read count
│           total_bases     — total sequenced bases
│           read_N50        — read N50 (bp): half of all bases are in reads ≥ this length
│           mean_length     — mean read length (bp)
│           mean_quality    — mean Phred quality score
│           gc_pct          — GC content (%)
│           estimated_depth — total_bases / genome_size (×)
│
├── ont_plot_metrics/
│   ├── ont_read_metrics.png
│   │     8-panel figure for samples passing the depth filter:
│   │       A  Read N50 histogram      B  Read N50 boxplot
│   │       C  Depth histogram         D  Depth boxplot
│   │           (vertical dashed line marks --min_read_depth threshold)
│   │       E  Total reads histogram   F  Total reads boxplot
│   │       G  Mean read length histogram   H  Mean read length boxplot
│   └── ont_read_metrics_summary.tsv
│         Merged TSV of all passing samples' read_stats rows
│
│  ── Assembly outputs (hybracter — default) ────────────────────────────────────
│
├── all_assemblies/                          ← pass this to enteric-typer
│   └── {sample}.fasta
│         Chromosome + plasmid sequences merged into a single FASTA per sample.
│         Complete assemblies: chromosome sequence(s) then plasmid contigs (if any).
│         Incomplete assemblies: draft contig set.
│
├── hybracter_output/
│   ├── {sample}_chromosome.fasta
│   │     Chromosome sequence(s) only. Present for complete (circularised) assemblies;
│   │     absent for incomplete assemblies.
│   ├── {sample}_plasmid.fasta
│   │     Plasmid contig sequences. Present only when plasmids were detected and
│   │     circularised; absent when no plasmids were found.
│   ├── {sample}_incomplete.fasta
│   │     Draft contig set for samples where the chromosome could not be circularised.
│   │     Present only for incomplete assemblies.
│   └── {sample}_hybracter_summary.tsv
│         Assembly outcome table (one row per sample):
│           sample            — sample name
│           complete          — TRUE / FALSE (chromosome circularised)
│           chromosome_length — assembled chromosome length (bp)
│           plasmid_count     — number of plasmid contigs detected
│
│  ── Assembly outputs (flye / dragonflye / unicycler) ──────────────────────────
│
├── assemblies/
│   └── {sample}.fasta        ← assembled contigs (pass to enteric-typer)
│
├── flye_assembly/            (if --assembler flye)
│   └── {sample}_flye_info.txt    — per-contig assembly statistics from Flye
│
├── dragonflye_assembly/      (if --assembler dragonflye)
│   └── {sample}_dragonflye.log   — DragonFlye run log
│
├── unicycler_assembly/       (if --assembler unicycler)
│   └── {sample}_unicycler.log    — Unicycler run log
│
│  ── Pipeline info ─────────────────────────────────────────────────────────────
│
└── pipeline_info/
    ├── timeline.html   — per-task runtime and resource usage
    ├── report.html     — summary execution report
    └── dag.svg         — workflow directed acyclic graph
```

> **Hybracter: complete vs incomplete assemblies.**
> The `{sample}_hybracter_summary.tsv` in `hybracter_output/` reports whether each
> assembly is complete (chromosome circularised) or incomplete (draft contigs). For
> complete assemblies, chromosome and plasmid FASTAs are saved separately alongside
> the merged `all_assemblies/{sample}.fasta`. For incomplete assemblies, only the
> draft contig FASTA is retained.

---

## Execution profiles

| Profile | Use case |
|---|---|
| `conda` | Local workstation with conda/mamba |
| `mamba` | Same as conda but with faster env solving |
| `arm64` | **Add on Apple Silicon (M1/M2/M3/M4)** — forces osx-64 conda envs via Rosetta 2 |
| `docker` | Local with Docker Desktop |
| `singularity` | HPC with Singularity/Apptainer |
| `slurm` | SLURM HPC executor (combine with another profile: `-profile conda,slurm`) |
| `pbs` | PBS/Torque HPC executor |

### HPC example

```bash
nextflow run assemble.nf -c assemble.config -profile singularity,slurm \
    --input_dir /path/to/fastq/ \
    --outdir    assembly_results/
```

### macOS Apple Silicon

```bash
CONDA_SUBDIR=osx-64 nextflow run assemble.nf -c assemble.config -profile conda,arm64 \
    --input_dir /path/to/fastq/ \
    --hybracter_no_medaka true
```

> **`--hybracter_no_medaka true` is required on macOS Apple Silicon.**
> Hybracter uses [medaka](https://github.com/nanoporetech/medaka) for neural network
> consensus polishing as its final assembly step. Medaka has a hard dependency on
> **OpenSSL 1.1.x**, which conflicts with the OpenSSL 3.x libraries present in macOS
> conda environments — even under Rosetta 2 (osx-64) emulation. Conda cannot resolve
> an environment that satisfies both constraints, so medaka's internal conda environment
> fails to build. Passing `--hybracter_no_medaka true` skips the polishing step:
> hybracter still runs Flye assembly, Plassembler plasmid recovery, and chromosome
> circularisation — you simply do not get the final medaka polishing pass.
>
> **If medaka polishing is required**, run the pipeline on Linux or an HPC cluster
> (omit `--hybracter_no_medaka` — it defaults to `false`). With high-accuracy Dorado
> basecalling (Q20+ reads), the impact of skipping medaka polishing on consensus
> accuracy is generally minimal.

> `CONDA_SUBDIR=osx-64` is required on Apple Silicon. Nextflow's conda integration
> does not consistently pass `--platform` to libmamba, so setting this environment
> variable ensures Rosetta 2 (x86_64) conda environments are created. Environments
> are cached after the first run.

---

## Cleaning up

```bash
# Remove Nextflow temporary files (safe once you are happy with results)
rm -rf work/ .nextflow/ .nextflow.log*
```

Keep `work/` if you want to use `-resume` to restart from a checkpoint.

---

## Citation

If you use ont-assembler, please also cite the underlying tools:

- **Hybracter**: Bouras et al. (2024) Microbial Genomics 10(5)
- **Flye**: Kolmogorov et al. (2019) Nature Biotechnology 37:540–546
- **DragonFlye**: github.com/rpetit3/dragonflye
- **Unicycler**: Wick et al. (2017) PLOS Computational Biology 13(6)
- **seqkit**: Shen et al. (2016) PLOS ONE 11(10):e0163962
- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology 35:316–319
