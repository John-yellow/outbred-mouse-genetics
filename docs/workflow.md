# Workflow

The workflow begins with provider-supplied quality-controlled paired-end FASTQ files (NGDC GSA `CRA046460`) listed in `config/sample_manifest.tsv`. `legacy/v1/scripts_legacy/01.run_fastp.sh` is retained as an optional preprocessing script for raw reads.

## Run order

Run from the repository root after copying config/pipeline.env.example to config/pipeline.env and editing it:

1. scripts/01_align.sh — BWA-MEM2 v2.3 alignment to the configured GRCm39 reference (checks for `${REFERENCE_FASTA}.bwt.2bit.64` and builds the BWA-MEM2 index if missing); coordinate sorting and indexing.
2. scripts/02_mark_duplicates.sh — FastDup 1.0.0 marks duplicates and retains them in BAM with the duplicate flag.
3. Run interval calling with the configured paths and resources:

       source config/pipeline.env
       python3 scripts/03_call_windows.py          --reference "$REFERENCE_FASTA"          --manifest "$SAMPLE_MANIFEST"          --bam-dir "$MARKED_DIR"          --ploidy config/ploidy.tsv          --output-dir "$CALL_DIR"          --bcftools "$BCFTOOLS"          --samtools "$SAMTOOLS"          --jobs "$CALL_JOBS"          --call-threads "$CALL_THREADS"

   The script expects files named SAMPLE.marked.bam and their indexes under MARKED_DIR. It writes 10-Mb windows for autosomes and X non-PAR and concatenates them in reference order.
4. scripts/04_filter_sets.sh — normalize, mask low-quality sample genotypes, recalculate allele statistics, and write candidate and main sets.
5. scripts/05_annotate_vep.sh — annotate the candidate set once with VEP 115.2/cache 115 and join CSQ to candidate and main callsets by CHROM/POS/REF/ALT.
6. Prepare the QC work directory described below, then run scripts/06_qc_summary.py WORKDIR.
7. Run scripts/07_known_novel_spectrum.py --bcf ANNOTATED_MAIN_BCF --output-dir OUTPUT_DIR to summarize exact dbSNP 150 matches and plot the allele-frequency spectrum.
8. `scripts/08_convert_release_vcfs.sh` — convert the annotated BCF callsets to BGZF-compressed `VCF.gz` with `CSI` indexes and validate output integrity via `scripts/validate_release_vcf.py`.

## QC work directory

scripts/06_qc_summary.py expects:

- state/run_config.json, based on config/qc_run_config.example.json;
- inputs/sample_manifest.tsv, inputs/ploidy.tsv, and inputs/reference.fa.fai;
- results/main_variants_vep115_GRCm39.bcf and its index;
- state/, logs/, stats/, and plots/ directories.

The run configuration records the expected main record count, sample count, and tool executable names. Paths in the configuration are relative to WORKDIR unless absolute.

## Production settings recorded for this dataset

- Reference: GRCm39 RefSeq GCF_000001635.27; FASTA SHA-256 ca6f2442484586b1cd8627e69394eabfe1fdb01f748bd1f089058ff3c0fa1b47.
- Alignment: BWA-MEM2 v2.3 release binary, -t 32 -K 50000000; samtools 1.10 sort -@ 2 -m 768M -l 3. BAM @PG records the upstream 2.2.1 version string.
- Duplicate marking: FastDup 1.0.0, upstream commit 7b3f62587283a257fd38c5c20ceeda0364c285ef; local patch and sanitized build record are in patches/. Production settings: SUM_OF_BASE_QUALITIES, 32 threads, read-name regex disabled, BAI index; duplicates are marked, not removed.
- Calling: bcftools/htslib 1.10.2; -q 20 -Q 20 -E -d 500 --ff UNMAP,SECONDARY,QCFAIL,DUP -p -m 3 -F 0.25; -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD; bcftools call -mAv -f GQ,GP with sample-sex and ploidy files. GRCm39 PAR X:168,752,755–169,376,592 (1-based inclusive) is excluded; it is not diploid-called.
- Filtering: mask a GT when FORMAT/DP or FORMAT/GQ is absent, DP<5, or GQ<20; recalculate AC/AN/AF and missingness. Candidate: biallelic SNP, QUAL≥30, 0<AC<AN. Main adds F_MISSING≤0.05. Neither set has a MAF floor.
- Annotation: VEP 115.2, cache 115, GRCm39, GENCODE M38, dbSNP 150. Original candidate and main callset fields are retained.

These thread and memory values describe the production run. For other machines, resource allocations can change while evidence-backed study thresholds and sample/ploidy settings remain fixed.
