# ont-assembler

🌐 **Language / Langue / Língua:** [English](README.md) | [Français](README_fr.md) | Português

---

Um pipeline Nextflow para montagem de leituras longas Oxford Nanopore (ONT) em genomas.

**Etapas do pipeline:**
1. **Controlo de qualidade das leituras** — métricas por amostra (N50, profundidade, total de leituras, comprimento médio) via `seqkit stats`
2. **Filtro de profundidade** — amostras abaixo de `--min_read_depth` (predefinição: 20×) são ignoradas com um aviso
3. **Figura de métricas de leitura** — figura com 8 painéis (histogramas + diagramas de caixa) para as amostras aprovadas
4. **Montagem** — [Hybracter](https://github.com/gbouras13/hybracter) por predefinição, com [Flye](https://github.com/fenderglass/Flye) e [Raven](https://github.com/lbcb-sci/raven) como alternativas

Os FASTAs montados estão prontos para serem transmitidos directamente ao [enteric-typer](https://github.com/efosternyarko/enteric-typer).

### Exemplo de figura de métricas de leitura (n = 7 amostras de E. coli)

![Exemplo de métricas ONT](assets/ont_read_metrics_example.png)

*Figura com 8 painéis: N50 das leituras, profundidade média de leitura, total de leituras e comprimento médio das leituras — cada um apresentado como histograma (esquerda) e diagrama de caixa (direita). A linha vermelha a tracejado no painel C assinala o limiar `--min_read_depth` (predefinição: 20×); as amostras abaixo deste valor são excluídas da montagem.*

---

## Instalação

### Requisitos

- [Nextflow](https://nextflow.io) ≥ 23.04
- [Conda](https://conda-forge.org/miniforge/) ou Mamba (recomendado)
- Java 17+ (`java -version`)

### Clonar o repositório

```bash
git clone https://github.com/efosternyarko/ont-assembler
cd ont-assembler
```

### Hybracter: configuração única

O Hybracter descarrega as suas bases de dados internas (~2 GB) na primeira utilização. O Nextflow cria o ambiente conda do Hybracter durante a primeira execução (as montagens falhadas são ignoradas via `errorStrategy = 'ignore'`). Após essa primeira execução, localize o ambiente e instale as bases de dados:

```bash
# 1. Encontrar o ambiente Hybracter criado pelo Nextflow
HYBRACTER_ENV=$(for e in work/conda/env-*/; do [ -f "${e}bin/hybracter" ] && echo "$e" && break; done)

# 2. Instalar as bases de dados — Linux / Mac Intel
conda run --prefix "$HYBRACTER_ENV" hybracter install

# 2. Instalar as bases de dados — macOS Apple Silicon (M1 e superior)
CONDA_SUBDIR=osx-64 conda run --prefix "$HYBRACTER_ENV" hybracter install
```

Depois, volte a executar o pipeline — as bases de dados ficam guardadas no ambiente entre execuções.

---

## Início rápido

```bash
# Hybracter (predefinição) — Linux / HPC / WSL2
nextflow run assemble.nf -c assemble.config -profile conda \
    --input_dir /caminho/para/fastq/ \
    --outdir    assembly_results_hybracter/

# macOS Apple Silicon (M1 e superior)
CONDA_SUBDIR=osx-64 nextflow run assemble.nf -c assemble.config -profile conda,arm64 \
    --input_dir /caminho/para/fastq/ \
    --outdir    assembly_results_hybracter/ \
    --hybracter_no_medaka true

# Ficheiro de amostras (CSV com colunas: id,reads)
nextflow run assemble.nf -c assemble.config -profile conda \
    --samplesheet samples.csv \
    --outdir      assembly_results_hybracter/

# Montadores alternativos
nextflow run assemble.nf -c assemble.config -profile conda \
    --input_dir /caminho/para/fastq/ \
    --assembler flye \
    --outdir    assembly_results_flye/

nextflow run assemble.nf -c assemble.config -profile conda \
    --input_dir /caminho/para/fastq/ \
    --assembler raven \
    --outdir    assembly_results_raven/
```

> Cada montador tem o seu próprio directório de saída para que execuções paralelas nunca se sobreponham.
> Se `--outdir` for omitido, o directório predefinido é `assembly_results_<montador>`.

**Transmitir os FASTAs montados ao enteric-typer** (com Hybracter por predefinição):

```bash
nextflow run /caminho/para/enteric-typer/main.nf -profile conda \
    --input_dir assembly_results_hybracter/all_assemblies/ \
    --outdir    typing_results/
```

---

## Montadores

| `--assembler` | Ferramenta | Notas |
|---|---|---|
| `hybracter` | [Hybracter](https://github.com/gbouras13/hybracter) | **Predefinição.** Circulariza cromossomas e plasmídeos. Utiliza `--auto` para estimar automaticamente o tamanho do cromossoma. |
| `flye` | [Flye](https://github.com/fenderglass/Flye) | Modo `--nano-hq` (leituras Guppy 5+ / Dorado Q20+). |
| `raven` | [Raven](https://github.com/lbcb-sci/raven) | Montador rápido para leituras longas não corrigidas. Não requer tamanho de genoma. Produz também um grafo de montagem GFA. |

---

## Parâmetros

| Parâmetro | Predefinição | Descrição |
|---|---|---|
| `--input_dir` | `null` | Directório de ficheiros FASTQ / FASTQ.gz. O identificador da amostra corresponde ao nome do ficheiro sem extensão. |
| `--samplesheet` | `null` | CSV com colunas `id,reads` |
| `--outdir` | `assembly_results_<montador>` | Directório de saída (o nome do montador é incluído para evitar conflitos) |
| `--assembler` | `hybracter` | Ferramenta de montagem: `hybracter` \| `flye` \| `raven` |
| `--genome_size` | `5m` | Tamanho estimado do genoma para cálculo de profundidade (ex. `5m`, `4500000`) |
| `--min_read_depth` | `20` | Profundidade mínima estimada (×). Amostras menos profundas são ignoradas e excluídas da figura. |
| `--chromosome_size` | `2500000` | Hybracter: comprimento mínimo de contig (pb) para ser considerado cromossoma. Ignorado se `--hybracter_auto true`. |
| `--hybracter_auto` | `true` | Deixar o Hybracter estimar o tamanho do cromossoma automaticamente (`--auto`). |
| `--hybracter_no_medaka` | `false` | Ignorar o polimento medaka. **Obrigatório no macOS Apple Silicon** (conflito OpenSSL). Manter `false` em Linux/HPC. |
| `--max_cpus` | `16` | Número máximo de CPUs por processo |
| `--max_memory` | `128.GB` | Memória máxima por processo |

---

## Ficheiros de saída

```
assembly_results_<montador>/
│
│  ── Controlo de qualidade das leituras ────────────────────────────────────────
│
├── ont_read_qc/
│   └── {amostra}_read_stats.tsv
│         Métricas de leitura por amostra (TSV, uma linha de dados):
│           sample_id       — nome da amostra
│           num_reads       — número total de leituras
│           total_bases     — total de bases sequenciadas
│           read_N50        — N50 das leituras (pb)
│           mean_length     — comprimento médio das leituras (pb)
│           mean_quality    — pontuação de qualidade Phred média
│           gc_pct          — conteúdo GC (%)
│           estimated_depth — total_bases / tamanho_genoma (×)
│
├── ont_plot_metrics/
│   ├── ont_read_metrics.png
│   │     Figura com 8 painéis para amostras que passam o filtro de profundidade:
│   │       A  Histograma N50          B  Diagrama de caixa N50
│   │       C  Histograma profundidade D  Diagrama de caixa profundidade
│   │           (linha vertical a tracejado = limiar --min_read_depth)
│   │       E  Histograma total leituras   F  Diagrama de caixa total leituras
│   │       G  Histograma comprimento médio   H  Diagrama de caixa comprimento médio
│   └── ont_read_metrics_summary.tsv
│         TSV combinado de todas as linhas read_stats das amostras aprovadas
│
│  ── Montagem (hybracter — predefinição) ───────────────────────────────────────
│
├── all_assemblies/                          ← transmitir ao enteric-typer
│   └── {amostra}.fasta
│
├── hybracter_output/
│   ├── {amostra}_chromosome.fasta
│   ├── {amostra}_plasmid.fasta
│   ├── {amostra}_incomplete.fasta
│   ├── {amostra}_hybracter_summary.tsv
│   └── {amostra}_flye.gfa             — grafo de montagem Flye (visualizar no Bandage)
│
│  ── Montagem (flye / raven) ───────────────────────────────────────────────────
│
├── assemblies/
│   └── {amostra}.fasta
│
├── flye_assembly/            (se --assembler flye)
│   ├── {amostra}_flye_info.txt
│   └── {amostra}_flye.gfa
│
├── raven_assembly/           (se --assembler raven)
│   └── {amostra}_raven.gfa
│
├── assembly_timing/
│   └── {amostra}_timing.tsv
│
└── assembly_timing_summary.tsv
│
│  ── Informações do pipeline ───────────────────────────────────────────────────
│
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

> **Hybracter: montagens completas e incompletas.**
> O ficheiro `{amostra}_hybracter_summary.tsv` indica se cada montagem é completa
> (cromossoma circularizado) ou incompleta (contigs em rascunho).

---

## Perfis de execução

| Perfil | Caso de utilização |
|---|---|
| `conda` | Estação de trabalho local com conda/mamba |
| `mamba` | Igual ao conda mas com resolução de ambiente mais rápida |
| `arm64` | **Adicionar no Apple Silicon (M1 e superior)** — força ambientes conda osx-64 via Rosetta 2 |
| `docker` | Local com Docker Desktop |
| `singularity` | Cluster HPC com Singularity/Apptainer |
| `slurm` | Executor HPC SLURM (combinar com outro perfil: `-profile conda,slurm`) |
| `pbs` | Executor HPC PBS/Torque |

### Exemplo HPC

```bash
nextflow run assemble.nf -c assemble.config -profile singularity,slurm \
    --input_dir /caminho/para/fastq/
```

### macOS Apple Silicon

```bash
CONDA_SUBDIR=osx-64 nextflow run assemble.nf -c assemble.config -profile conda,arm64 \
    --input_dir /caminho/para/fastq/ \
    --hybracter_no_medaka true
```

> **`--hybracter_no_medaka true` é obrigatório no macOS Apple Silicon.**
> O Hybracter utiliza o [medaka](https://github.com/nanoporetech/medaka) para polimento por
> consenso com rede neuronal na etapa final de montagem. O medaka depende do **OpenSSL 1.1.x**,
> o que entra em conflito com as bibliotecas OpenSSL 3.x presentes nos ambientes conda do macOS
> — mesmo sob emulação Rosetta 2 (osx-64). Passar `--hybracter_no_medaka true` ignora esta etapa
> de polimento: o Hybracter continua a efectuar a montagem Flye, a recuperação de plasmídeos com
> o Plassembler e a circularização do cromossoma.
>
> **Se o polimento medaka for necessário**, execute o pipeline em Linux ou num cluster HPC.
> Com basecalling Dorado de alta precisão (leituras Q20+), o impacto de omitir o polimento
> medaka na precisão do consenso é geralmente mínimo.

> `CONDA_SUBDIR=osx-64` é necessário no Apple Silicon para garantir a criação de ambientes
> conda Rosetta 2 (x86_64). Os ambientes ficam em cache após a primeira execução.

---

## Limpeza

```bash
# Remover ficheiros temporários do Nextflow (seguro após verificar os resultados)
rm -rf work/ .nextflow/ .nextflow.log*
```

Mantenha `work/` se pretender utilizar `-resume` para retomar a partir de um ponto de controlo.

---

## Citação

Se utilizar o ont-assembler, cite também as ferramentas subjacentes:

- **Hybracter**: Bouras et al. (2024) Microbial Genomics 10(5)
- **Flye**: Kolmogorov et al. (2019) Nature Biotechnology 37:540–546
- **Raven**: Vaser & Šikić (2021) Nature Computational Science 1:332–336
- **seqkit**: Shen et al. (2016) PLOS ONE 11(10):e0163962
- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology 35:316–319
