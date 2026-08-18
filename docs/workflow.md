# Workflow notes

The workflow is based on the ten analysis scripts supplied with the study and follows the methods described in the manuscript.

## Stages 1–5: sequence processing and variants

1. `01.run_fastp.sh` discovers paired-end reads in one directory per sample and runs fastp. It writes cleaned reads and HTML/JSON reports.
2. `02.run_bwa_s2b_sorted.sh` aligns cleaned reads to the GRCm39 primary assembly with BWA-MEM, retains chromosomes `1`–`19`, `X`, `Y`, and `MT`, and creates sorted/indexed BAM files.
3. `03.bcftools_SNPcalling.sh` builds non-overlapping windows from the first BAM header, jointly calls variants across all BAMs with bcftools, and concatenates interval VCFs.
4. `04.vcf_qc.sh` retains chromosomes `1`–`19` and `X`, normalizes records, recalculates `AN`, `AC`, and `MAF`, trims unused alleles, and applies the study filters:

   ```text
   QUAL >= 30 && MAF > 0.01 && F_MISSING == 0 && DP >= 5 && GQ >= 20
   ```

5. `05.Vep_annotate.sh` performs offline mouse VEP annotation using assembly GRCm39 and cache release 115.

## Stages 6–10: population genetics

6. `06.Genetic_Val.sh` assigns variant IDs as `CHROM:POS`, performs LD pruning with `50 5 0.2`, creates a PLINK binary dataset, runs PCA, estimates pairwise FST, and runs ADMIXTURE at `K = 3`.
7. `07.split_fam.sh` makes population-specific PLINK text datasets at `MAF >= 0.01` using the sample-group table.
8. `08.den.R` reads population `.map` files and creates chromosome-density plots with CMplot.
9. `09.PCA.R` reads PLINK eigenvectors/eigenvalues and creates the PC1/PC2 plot.
10. `10.GBC.R` reads the ADMIXTURE Q matrix and PLINK FAM file, orders individuals by population, and writes the ADMIXTURE plot and group-level ancestry table.

## Important assumptions

- VCF sample IDs match `config/sample_groups.tsv`.
- Population labels are exactly `ICR1`, `ICR2`, and `KM`.
- The first two PLINK FAM columns are FID and IID; the scripts use the FID as the population label after the optional ID update in stage 6.
- The full workflow is computationally intensive. The demo directory is not designed to run all ten stages.
