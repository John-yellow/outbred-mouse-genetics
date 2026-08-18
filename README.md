# Outbred Mouse Genetics

Reproducible scripts for whole-genome resequencing, variant discovery, variant annotation, and population-genetic characterization of three commercial outbred mouse stocks: ICR1, ICR2, and KM.

This repository is the public code companion for the manuscript *Whole-genome resequencing dataset for characterizing genomic variation and population structure in three commercial outbred mouse stocks in China*. The GitHub repository is `outbred-mouse-genetics`:

<https://github.com/John-yellow/outbred-mouse-genetics>

The workflow follows the ten scripts used for the study:

| Stage | Script | Main result |
|---|---|---|
| 1 | `scripts/01.run_fastp.sh` | Quality-controlled paired-end FASTQ files and fastp reports |
| 2 | `scripts/02.run_bwa_s2b_sorted.sh` | Coordinate-sorted, indexed BAM files |
| 3 | `scripts/03.bcftools_SNPcalling.sh` | Joint raw VCF and index |
| 4 | `scripts/04.vcf_qc.sh` | Filtered high-confidence SNP VCF and index |
| 5 | `scripts/05.Vep_annotate.sh` | VEP-annotated SNP VCF |
| 6 | `scripts/06.Genetic_Val.sh` | LD-pruned genotypes, PCA, pairwise FST, and ADMIXTURE input/output |
| 7 | `scripts/07.split_fam.sh` | Population-specific PLINK text datasets |
| 8 | `scripts/08.den.R` | SNP-density plots |
| 9 | `scripts/09.PCA.R` | PCA plot |
| 10 | `scripts/10.GBC.R` | ADMIXTURE bar plot and group-level ancestry table |

## Important data policy

The repository contains code, configuration templates, documentation, and small example tables only. It does **not** contain the study's complete FASTQ, BAM, VCF/BCF, reference genome, VEP cache, PLINK binary files, or large result files. These files are intentionally excluded by `.gitignore` and must be obtained from the appropriate data or reference repositories before a full run.

The `data/demo/` directory is a miniature, non-production example for checking file formats and understanding the expected inputs. It is not a substitute for the 90-sample study data and is not intended to reproduce the manuscript's numerical results. Do not replace the tracked demo sample tables with real sample identifiers before uploading this repository.

## Quick start for the full analysis

From the repository root:

```bash
mamba env create -f envs/conda.yml
mamba activate outbred-mouse-genetics
```

1. Download or otherwise provide the GRCm39 primary-assembly FASTA and its index files in private storage.
2. Install or provide the offline VEP mouse cache (release 115) and set `VEP_CACHE_DIR` in the local configuration.
3. Keep the real FASTQ files and sample tables in private storage. Create a local, ignored configuration file that points to them; do not commit real sample identifiers or machine-specific paths.
4. Review CPU and memory settings in the local configuration.
5. Run the stages from the repository root:

```bash
bash scripts/01.run_fastp.sh
bash scripts/02.run_bwa_s2b_sorted.sh
bash scripts/03.bcftools_SNPcalling.sh
bash scripts/04.vcf_qc.sh
bash scripts/05.Vep_annotate.sh
bash scripts/06.Genetic_Val.sh
bash scripts/07.split_fam.sh

Rscript scripts/08.den.R
Rscript scripts/09.PCA.R
Rscript scripts/10.GBC.R
```

For a private run, a convenient pattern is:

```bash
cp config/pipeline.env config/pipeline.local
# Edit config/pipeline.local: set SAMPLE_MAP, SAMPLE_GROUPS,
# REFERENCE_FASTA, VEP_FASTA, and VEP_CACHE_DIR to private paths.
export CONFIG_FILE="$PWD/config/pipeline.local"
source "$CONFIG_FILE"
```

The file `config/pipeline.local` is ignored by Git. The tracked `config/pipeline.env`, `config/sample_map.tsv`, and `config/sample_groups.tsv` remain public templates only.

Every shell script resolves paths relative to the repository root, so it can also be called from another working directory. The scripts source the file selected by `CONFIG_FILE` (default: `config/pipeline.env`); use the ignored `config/pipeline.local` for machine-specific paths.

## Expected input conventions

Raw reads are arranged as one directory per sample:

```text
data/raw/
└── SAMPLE_ID/
    ├── SAMPLE_ID_R1.fq.gz
    └── SAMPLE_ID_R2.fq.gz
```

`config/sample_map.tsv` contains two whitespace-separated columns, `raw_name` and `clean_name`. `config/sample_groups.tsv` contains `sample_id`, `population`, and `sex`; the `sample_id` values must match the sample IDs in the VCF/PLINK files. The default population labels are `ICR1`, `ICR2`, and `KM`.

The reference FASTA must use the GRCm39 chromosome names expected by the scripts: `1`–`19`, `X`, `Y`, and `MT`. The variant-calling and filtering thresholds are documented in `docs/workflow.md` and can be adjusted in the relevant script/configuration when a different study design requires it.

## Output layout

Large generated files are written below ignored directories:

```text
data/clean/       fastp output
data/bam/         alignment output
data/vcf/         raw, filtered, and annotated VCF output
results/genetics/ PLINK, PCA, FST, and ADMIXTURE output
results/figures/  R-generated figures
work/             temporary lists and QC intermediates
```

Do not force-add these directories to Git. Deposit complete sequencing and variant files in the public domain-specific repositories described by the manuscript and record their final accessions in the manuscript/data-availability documentation.

## Reproducibility notes

The software versions in `envs/conda.yml` mirror the versions reported in the manuscript where practical. VEP cache data are distributed separately from the code environment and are not bundled here. For a publication release, create a Git tag after the accession numbers and final configuration have been confirmed.

## Citation and code availability

The ready-to-paste code-availability wording is in `docs/manuscript-code-availability.txt`. Please cite the associated manuscript and the data-accession records when reusing the workflow or the resulting genomic resource.

## Upload to GitHub

After reviewing the templates, without adding private sample tables or large data files:

```bash
git add .
git commit -m "Add outbred mouse genetics workflow"
git push -u origin main
```

The repository is initialized on the `main` branch locally and uses the manuscript's GitHub URL.
