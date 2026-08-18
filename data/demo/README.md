# Miniature demo files

These files are deliberately tiny, synthetic format examples. They are included to show how the downstream tables are shaped without distributing any complete FASTQ, BAM, VCF, reference, or study result file.

- `sample_map.tsv.example` shows the stage-1 sample map.
- `raw.vcf.example` and `reference.fa.example` illustrate the text formats only; they are not GRCm39 data and are not sufficient to run the full pipeline.
- `pca.*.example`, `admixture.*.example`, and `*.map.example` are miniature inputs for understanding the R scripts' expected downstream files.

The full workflow still requires the study data, the GRCm39 reference, VEP cache, and the software listed in `envs/conda.yml`.
