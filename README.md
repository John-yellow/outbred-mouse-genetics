# Whole-genome resequencing of three commercial outbred mouse stocks

This repository contains the whole-genome resequencing variant calling, annotation, and quality-control workflow for 90 commercial outbred mice representing ICR1, ICR2, and KM (30 mice per stock; 15 males and 15 females). Quality-controlled Illumina NovaSeq X Plus paired-end reads are archived in NGDC GSA (`CRA046460`, BioProject `PRJCA069117`), and the resulting biallelic SNP callsets (candidate set: 16,529,182 SNPs; main set: 13,076,115 SNPs; BGZF-compressed VCF with CSI indexes) are publicly available in NGDC GVM under accession `GVM001647` (`https://ngdc.cncb.ac.cn/gvm/getProjectDetail?project=GVM001647`).

## Workflow

The workflow is documented in `docs/workflow.md`. It performs BWA-MEM2 alignment (`scripts/01_align.sh`), FastDup mark-only duplicate processing (`scripts/02_mark_duplicates.sh`), interval-based `bcftools` calling (`scripts/03_call_windows.py`), per-genotype DP/GQ masking and candidate/main SNP filtering (`scripts/04_filter_sets.sh`), Ensembl VEP 115 annotation (`scripts/05_annotate_vep.sh`), sample/variant QC and plotting (`scripts/06_qc_summary.py`, `scripts/07_known_novel_spectrum.py`), and release VCF.gz + CSI conversion and validation (`scripts/08_convert_release_vcfs.sh`, `scripts/validate_release_vcf.py`). Earlier v1 scripts remain in `legacy/v1/` for reference.

## Getting started

1. Install the software and obtain the GRCm39 reference and matching Ensembl VEP 115 cache described in envs/README.md.
2. Copy config/pipeline.env.example to config/pipeline.env and edit the paths and resource allocations.
3. Copy config/sample_manifest.example.tsv to config/sample_manifest.tsv. Supply one row per sample with stable sample ID, stock/group, sex, and paired QC FASTQ paths. Keep row order stable.
4. From the repository root, run scripts/01_align.sh and scripts/02_mark_duplicates.sh.
5. Run scripts/03_call_windows.py with the required input and resource arguments shown in docs/workflow.md, then run scripts/04_filter_sets.sh and scripts/05_annotate_vep.sh.
6. For QC summaries, prepare the work-directory layout in docs/workflow.md and run scripts/06_qc_summary.py. Run scripts/07_known_novel_spectrum.py on the annotated main BCF to generate the dbSNP 150 match summary and allele-frequency spectrum plot.

The production thread and memory allocations are recorded in docs/workflow.md. Resource settings can be adjusted for a different machine; the evidence-backed quality thresholds, sample order, ploidy, chromosome scope, and retained fields should remain fixed to reproduce this dataset.

## Inputs and outputs

Reads, the GRCm39 reference, VEP cache files, callsets, and large QC tables are supplied by their respective data sources, not Git. Example manifests and configuration are provided here. Never commit local sequencing paths or credentials.

`scripts/08_convert_release_vcfs.sh` converts the annotated BCF callsets to BGZF-compressed `.vcf.gz` with `.csi` indexes and runs `scripts/validate_release_vcf.py` to verify sample order, indexed counts, full record and genotype equivalence, X-ploidy concordance, `INFO/CSQ` headers, CSI region queries, and SHA-256 checksums.

Code is licensed under MIT. A separate license must be selected for genomic data before deposit.
