# Data and repository policy

This code repository is deliberately separate from the large genomic data products.

## Included

- Shell and R scripts used by the analysis workflow.
- Small configuration templates and example tables.
- Documentation of file formats, parameters, and expected software.

## Not included

- Complete raw or cleaned FASTQ files.
- BAM/CRAM alignment files and indexes.
- Raw, filtered, or annotated VCF/BCF files and indexes.
- The GRCm39 reference genome and indexes.
- The offline VEP cache.
- PLINK binary datasets, ADMIXTURE matrices, and generated figures/results.

The `.gitignore` file protects these classes of files by default. Do not use `git add -f` on a large data file merely to make a local run convenient. Keep full data in the appropriate public archive or in a separately managed storage location.

The example files under `data/demo/` are tiny format examples and are not the study data. They must not be interpreted as a validation of the manuscript's reported values.

Do not overwrite the tracked demo rows in `config/sample_map.tsv` or `config/sample_groups.tsv` with real sample identifiers before pushing the repository. For a private analysis, keep those tables outside the public repository (or in an ignored local file) and point to them through `config/pipeline.local`/`CONFIG_FILE`.

The current project-level data record is CNCB-NGDC BioProject `PRJCA069117`. This accession identifies the BioProject only. The GSA and GVM accessions are not yet available and are not guessed in this repository.
