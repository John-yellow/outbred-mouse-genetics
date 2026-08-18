# Outbred Mouse Genetics

This repository contains the computational workflow and analysis code used for whole-genome resequencing and population-genetic analyses of three commercial outbred mouse populations.

The workflow covers read quality control, alignment to the *Mus musculus* GRCm39 reference genome, SNP calling and quality control, variant annotation, and downstream population-genetic analyses. This is a code-and-workflow release: genomic data and large analysis outputs are not distributed through GitHub.

## Overview

The repository provides:

- the ten shell and R scripts used by the analysis workflow;
- a portable configuration framework and sample-table schemas;
- a Conda environment definition for the command-line tools and R packages;
- documentation of the workflow, inputs, outputs, and reproducibility requirements.

It does not provide raw reads, processed genomic data, reference resources, or complete analysis results. Users must provide their own appropriately formatted input files and external reference resources before running the workflow.

## Study context

| Field | Value |
|---|---|
| Organism | *Mus musculus* |
| Population type | Commercial outbred mouse populations |
| Reference assembly | GRCm39 |
| BioProject | CNCB-NGDC BioProject accession **PRJCA069117** |

## Workflow

The operational order is:

```text
paired-end FASTQ
    -> 01.run_fastp.sh
    -> 02.run_bwa_s2b_sorted.sh
    -> 03.bcftools_SNPcalling.sh
    -> 04.vcf_qc.sh
    -> 05.Vep_annotate.sh
    -> 06.Genetic_Val.sh
    -> 07.split_fam.sh
    -> 08.den.R / 09.PCA.R / 10.GBC.R
```

The script-level inputs and outputs are documented in [`docs/workflow.md`](docs/workflow.md). The configuration variables and sample-table schemas are documented in [`docs/configuration.md`](docs/configuration.md).

## Repository structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   ├── pipeline.env
│   ├── sample_map.tsv
│   └── sample_groups.tsv
├── docs/
│   ├── code-availability.md
│   ├── configuration.md
│   └── workflow.md
├── envs/
│   └── conda.yml
└── scripts/
    ├── 01.run_fastp.sh
    ├── 02.run_bwa_s2b_sorted.sh
    ├── 03.bcftools_SNPcalling.sh
    ├── 04.vcf_qc.sh
    ├── 05.Vep_annotate.sh
    ├── 06.Genetic_Val.sh
    ├── 07.split_fam.sh
    ├── 08.den.R
    ├── 09.PCA.R
    ├── 10.GBC.R
    └── lib/common.sh
```

## Requirements

The command-line tools and R packages used by the scripts are listed in [`envs/conda.yml`](envs/conda.yml). Create the environment with:

```bash
conda env create -f envs/conda.yml
conda activate outbred-mouse-genetics
```

The workflow also requires external resources that are not bundled here:

- paired-end FASTQ files;
- the *Mus musculus* GRCm39 primary-assembly FASTA and the indexes required by the tools;
- the offline mouse VEP cache (release 115) and its FASTA resource;
- private sample metadata and population-group tables.

## Configuration

The tracked files in `config/` are portable templates and contain no individual-level sample rows. Before a private run, copy the pipeline configuration to an ignored local file and edit the paths and resource settings:

```bash
cp config/pipeline.env config/pipeline.local.env
export CONFIG_FILE="$PWD/config/pipeline.local.env"
source "$CONFIG_FILE"
```

Set `SAMPLE_MAP` and `SAMPLE_GROUPS` to private tab-delimited tables. The expected columns are described in [`docs/configuration.md`](docs/configuration.md). Keep study-specific sample identifiers and machine-specific paths outside the tracked files.

## Running the workflow

Run the shell stages in order from the repository root after configuring `CONFIG_FILE`:

```bash
bash scripts/01.run_fastp.sh
bash scripts/02.run_bwa_s2b_sorted.sh
bash scripts/03.bcftools_SNPcalling.sh
bash scripts/04.vcf_qc.sh
bash scripts/05.Vep_annotate.sh
bash scripts/06.Genetic_Val.sh
bash scripts/07.split_fam.sh
```

Run the R stages after the corresponding PLINK and ADMIXTURE outputs exist:

```bash
Rscript scripts/08.den.R \
  --input-dir "$MAP_DIR" --output-dir "$FIGURE_DIR"
Rscript scripts/09.PCA.R \
  --prefix "$PCA_PREFIX" --output-dir "$FIGURE_DIR"
Rscript scripts/10.GBC.R \
  --q-prefix "$Q_PREFIX" --fam "$FAM_FILE" --output-dir "$FIGURE_DIR"
```

The shell scripts resolve relative paths from the repository root. The R scripts accept the command-line paths shown above and otherwise use the documented relative defaults. The full workflow requires substantial compute and storage; this repository is not intended to run without the corresponding input data and external reference resources.

## Input data expectations

The first stage expects paired-end FASTQ files arranged by sample, for example:

```text
<raw-data-directory>/SAMPLE_ID/SAMPLE_ID_R1.fq.gz
<raw-data-directory>/SAMPLE_ID/SAMPLE_ID_R2.fq.gz
```

The private sample map is a two-column, tab-delimited table with one row per sample:

```text
raw_name    clean_name
```

The private sample-group table is a three-column, tab-delimited table with one row per sample:

```text
sample_id    population    sex
```

Sample IDs must match the VCF/PLINK identifiers used by the downstream scripts. The population labels expected by the current workflow are `ICR1`, `ICR2`, and `KM`. The tracked configuration files document these schemas without containing sample records.

## Outputs

The stages create local, ignored outputs including cleaned reads and reports, sorted/indexed BAM files, raw and filtered VCF files, annotated variants, PLINK datasets, PCA/FST/ADMIXTURE outputs, population-specific PLINK text files, and R-generated figures. Output locations are controlled by `config/pipeline.env` or its ignored local copy; none of these generated files are part of the public repository.

## Data note

This repository contains source code and workflow documentation only. Raw sequencing data, reference files, intermediate files, and full analysis outputs are not distributed through this Git repository.

The sequencing study associated with this workflow is registered under CNCB-NGDC BioProject accession **PRJCA069117**.

## Reproducibility

The repository provides the analysis workflow and software environment needed to reproduce the computational procedures when the corresponding input data and external reference resources are available. The repository alone is therefore not sufficient to reproduce numerical results. Tool versions and R package dependencies are recorded in [`envs/conda.yml`](envs/conda.yml), while workflow assumptions and stage details are recorded in [`docs/workflow.md`](docs/workflow.md).

## Citation

If you use this workflow, please cite the associated publication when available and reference this repository:

<https://github.com/John-yellow/outbred-mouse-genetics>

A final publication citation and DOI are not specified in this repository.

## License

This repository is released under the [MIT License](LICENSE).
