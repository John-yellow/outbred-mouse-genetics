# Data-release checklist

Verification checklist for releasing the 90-sample GRCm39 biallelic SNP callsets:

- [x] Convert annotated candidate and main BCF sources to `candidate_variants_vep115_GRCm39.vcf.gz` and `main_variants_vep115_GRCm39.vcf.gz` with `.csi` indexes (`scripts/08_convert_release_vcfs.sh`).
- [x] Validate both `.vcf.gz` + `.csi` pairs via `scripts/validate_release_vcf.py` (90-sample order, indexed counts, byte-exact full record/genotype equivalence, X non-PAR ploidy concordance, `INFO/CSQ` header, non-empty Chr1/ChrX CSI region queries, and SHA-256 checksums).
- [x] Publish both full-genotype callsets in NGDC GVM under accession `GVM001647` (`https://ngdc.cncb.ac.cn/gvm/getProjectDetail?project=GVM001647`, Open Access) alongside the quality-controlled sequencing reads in NGDC GSA (`CRA046460`, BioProject `PRJCA069117`).
