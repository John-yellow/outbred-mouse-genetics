# Configuration

The tracked configuration files describe the workflow schema. They are not populated with study-specific sample identifiers or server paths. Copy `config/pipeline.env` to an ignored local file for a private run:

```bash
cp config/pipeline.env config/pipeline.local.env
export CONFIG_FILE="$PWD/config/pipeline.local.env"
source "$CONFIG_FILE"
```

The shell scripts source the file selected by `CONFIG_FILE`; if it is unset, they use `config/pipeline.env`. Relative paths are resolved from the repository root by `scripts/lib/common.sh`. An absolute local configuration path is recommended when running from another directory.

## Pipeline paths and files

| Variable | Type | Default | Purpose |
|---|---|---|---|
| `RAW_DIR` | directory | `data/raw` | Per-sample paired-end FASTQ directories |
| `CLEAN_DIR` | directory | `data/clean` | fastp cleaned reads and reports |
| `BAM_DIR` | directory | `data/bam` | Sorted/indexed BAM files |
| `VCF_DIR` | directory | `data/vcf` | Raw, filtered, and annotated VCF files |
| `WORK_DIR` | directory | `work` | Temporary lists and stage intermediates |
| `GENETICS_DIR` | directory | `results/genetics` | PLINK, PCA, FST, and ADMIXTURE outputs |
| `FIGURE_DIR` | directory | `results/figures` | R-generated figures and tables |
| `POPULATION_DIR` | directory | `results/genetics/populations` | Population-specific PLINK text outputs |
| `MAP_DIR` | directory | `results/genetics/populations` | Map-file input for `08.den.R` |
| `PCA_PREFIX` | file prefix | `results/genetics/obm_pca` | PLINK PCA eigenvector/eigenvalue prefix |
| `Q_PREFIX` | file prefix | `results/genetics/obm_pruned_dataset.3` | ADMIXTURE Q-matrix prefix for `10.GBC.R` |
| `FAM_FILE` | file | `results/genetics/obm_pruned_dataset.fam` | PLINK FAM file for `10.GBC.R` |

The default data, results, reference, and work paths are local workspace conventions only. They do not imply that those directories or their contents are distributed by GitHub.

## Sample and reference inputs

| Variable | Type | Default | Purpose |
|---|---|---|---|
| `SAMPLE_MAP` | file | `config/sample_map.tsv` | Private two-column raw-to-clean sample map |
| `SAMPLE_GROUPS` | file | `config/sample_groups.tsv` | Private three-column sample-group table |
| `REFERENCE_FASTA` | file | `reference/Mus_musculus.GRCm39.dna.primary_assembly.fa` | GRCm39 FASTA for alignment, calling, and normalization |
| `VEP_CACHE_DIR` | directory | `reference/vep_cache` | Offline mouse VEP cache, release 115 |
| `VEP_FASTA` | file | `reference/Mus_musculus.GRCm39.dna.primary_assembly.fa` | FASTA used by VEP |
| `RAW_VCF` | file | `data/vcf/obm.raw.vcf.gz` | Stage 3 raw VCF |
| `FILTERED_VCF` | file | `data/vcf/obm_snps.vcf.gz` | Stage 4 filtered SNP VCF |
| `ANNOTATED_VCF` | file | `data/vcf/obm2026_snps_annotated.vcf.gz` | Stage 5 annotated VCF |
| `RENAMED_VCF` | file | `results/genetics/obm2026_snps_renamedID.vcf.gz` | Stage 6 variant-ID-renamed VCF |

The reference FASTA, FASTA indexes, BWA indexes, annotation database, and VEP cache must be obtained separately and configured locally. They are not included in this repository.

## Sample-table schemas

The tracked `config/sample_map.tsv` and `config/sample_groups.tsv` files contain comments describing their schemas but no rows. Private tables should be tab-delimited and contain one sample per row.

### `SAMPLE_MAP`

Use two columns without a header:

```text
raw_name    clean_name
```

`raw_name` is the directory name under `RAW_DIR` containing the paired FASTQ files. `clean_name` becomes the cleaned-read and downstream sample name.

### `SAMPLE_GROUPS`

Use three columns without a header:

```text
sample_id    population    sex
```

`sample_id` must match the sample identifier in the VCF/PLINK files. The current workflow expects the population labels `ICR1`, `ICR2`, and `KM`. The `sex` field is retained as table metadata; the downstream grouping uses the population field.

Study-specific identifiers must remain in the private table and must not be copied into the tracked templates.

## Chromosomes and resource settings

| Variable | Purpose |
|---|---|
| `TARGET_CHROMS` | Chromosomes retained during BAM filtering and variant calling |
| `QC_CHROMS` | Chromosomes retained by the VCF QC stage |
| `FASTP_THREADS`, `FASTP_JOBS` | fastp per-job threads and parallel sample jobs |
| `BWA_THREADS`, `VIEW_THREADS`, `SORT_THREADS`, `MAPPING_JOBS` | Alignment, BAM filtering/sorting, and mapping parallelism |
| `BCFTOOLS_THREADS`, `CALL_JOBS`, `REGION_SIZE` | Variant calling threads, region jobs, and window size |
| `VEP_FORKS` | VEP worker count |
| `ADMIXTURE_K`, `ADMIXTURE_CV`, `ADMIXTURE_THREADS` | ADMIXTURE model, cross-validation setting, and thread count |

The resource settings control execution cost. The analysis filters and population-genetic procedures remain encoded in the scripts and are documented in [`workflow.md`](workflow.md).

## R-stage paths

The R scripts do not source the shell configuration automatically. Pass the configured paths explicitly, as in the README:

```bash
Rscript scripts/08.den.R --input-dir "$MAP_DIR" --output-dir "$FIGURE_DIR"
Rscript scripts/09.PCA.R --prefix "$PCA_PREFIX" --output-dir "$FIGURE_DIR"
Rscript scripts/10.GBC.R --q-prefix "$Q_PREFIX" --fam "$FAM_FILE" --output-dir "$FIGURE_DIR"
```

Each R script also has a documented relative default and creates its output directory when needed.
