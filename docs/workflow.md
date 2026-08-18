# Workflow

This document describes the execution order and the principal file types used by the ten scripts. Paths are read from the shell configuration selected through `CONFIG_FILE`; the R scripts also accept explicit input and output paths.

## Execution order

```text
paired-end FASTQ
    ↓
01.run_fastp.sh
    ↓
clean FASTQ and QC reports
    ↓
02.run_bwa_s2b_sorted.sh
    ↓
sorted/indexed BAM
    ↓
03.bcftools_SNPcalling.sh
    ↓
raw variant VCF
    ↓
04.vcf_qc.sh
    ↓
filtered SNP VCF
    ↓
05.Vep_annotate.sh
    ↓
annotated SNP VCF
    ↓
06.Genetic_Val.sh and 07.split_fam.sh
    ↓
PLINK, PCA, FST, ADMIXTURE, and population-specific files
    ↓
08.den.R, 09.PCA.R, and 10.GBC.R
    ↓
population-genetic plots and summary tables
```

## Script reference

| Step | Script | Main input | Main output |
|---|---|---|---|
| 01 | `scripts/01.run_fastp.sh` | Paired-end FASTQ directories and the private sample map | Cleaned paired-end FASTQ files plus fastp HTML/JSON reports |
| 02 | `scripts/02.run_bwa_s2b_sorted.sh` | Cleaned FASTQ and GRCm39 FASTA | Coordinate-sorted, indexed BAM files |
| 03 | `scripts/03.bcftools_SNPcalling.sh` | Sorted BAM files and GRCm39 FASTA | Joint raw VCF and index |
| 04 | `scripts/04.vcf_qc.sh` | Raw VCF and GRCm39 FASTA | Filtered biallelic SNP VCF and index |
| 05 | `scripts/05.Vep_annotate.sh` | Filtered VCF, GRCm39 FASTA, and offline VEP cache | VEP-annotated VCF and optional index |
| 06 | `scripts/06.Genetic_Val.sh` | Annotated VCF and private population table | PLINK datasets, LD-pruned genotypes, PCA, pairwise FST, and ADMIXTURE files |
| 07 | `scripts/07.split_fam.sh` | Renamed VCF and private population table | Population-specific PLINK text datasets |
| 08 | `scripts/08.den.R` | Population-specific PLINK `.map` files | Chromosome-density PDF plots |
| 09 | `scripts/09.PCA.R` | PLINK eigenvector/eigenvalue files | PCA PDF and R workspace |
| 10 | `scripts/10.GBC.R` | ADMIXTURE Q matrix and PLINK FAM file | ADMIXTURE PDF and group-level ancestry tables |

## Stages 1–5: sequence processing and variants

1. `01.run_fastp.sh` discovers paired reads in one directory per sample and runs fastp. The private sample map connects each raw directory name to the cleaned-read sample name.
2. `02.run_bwa_s2b_sorted.sh` aligns reads to the GRCm39 primary assembly with BWA-MEM, retains the configured chromosome set, and creates sorted/indexed BAM files.
3. `03.bcftools_SNPcalling.sh` reads the BAM header from the first sorted BAM, creates non-overlapping windows, jointly calls variants with bcftools, and concatenates the interval VCFs.
4. `04.vcf_qc.sh` selects the configured chromosomes, normalizes records against the reference, recalculates `AN`, `AC`, and `MAF`, removes unused alternate alleles, and applies the filters encoded in the script:

   ```text
   QUAL >= 30 && MAF > 0.01 && F_MISSING == 0 && DP >= 5 && GQ >= 20
   ```

   The output is restricted to biallelic SNPs.
5. `05.Vep_annotate.sh` performs offline mouse VEP annotation for assembly GRCm39 using cache release 115.

## Stages 6–10: population genetics

6. `06.Genetic_Val.sh` assigns variant IDs as `CHROM:POS`, creates PLINK datasets, updates population labels from the private group table, performs LD pruning with `50 5 0.2`, runs PCA, estimates pairwise FST, and runs ADMIXTURE with the configured `K` value.
7. `07.split_fam.sh` uses the private group table to make population-specific PLINK text datasets for `ICR1`, `ICR2`, and `KM`, applying the MAF setting encoded in the script.
8. `08.den.R` reads population-specific PLINK map files and creates chromosome-density plots with CMplot.
9. `09.PCA.R` reads PLINK eigenvectors and eigenvalues and creates the PC1/PC2 plot.
10. `10.GBC.R` reads the ADMIXTURE Q matrix and PLINK FAM file, orders samples by population, and writes the ADMIXTURE plot and group-level ancestry tables.

## Input and output assumptions

- The reference is *Mus musculus* GRCm39 with chromosome names compatible with the configured target chromosome list.
- The sample map has two tab-delimited columns without a header: `raw_name` and `clean_name`.
- The population table has three tab-delimited columns without a header: `sample_id`, `population`, and `sex`. The sample IDs must match VCF/PLINK identifiers, and the population labels must be `ICR1`, `ICR2`, or `KM` for the current scripts.
- Intermediate and final files are written to the configured local data, results, reference, and work paths. Those paths are intentionally not represented by tracked directories in this repository.
- The workflow is computationally intensive and requires the external input data, reference indexes, and VEP cache described in the README and configuration guide.
