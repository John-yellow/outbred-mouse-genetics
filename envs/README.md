# Software installation and versions

Install Python 3.8 or newer and record the following software versions when reproducing the accepted workflow. The included Python scripts otherwise use only the standard library:

| Tool | Production version/settings |
|---|---|
| BWA-MEM2 | v2.3 release binary; BAM `@PG` records the upstream 2.2.1 version string |
| samtools | 1.10 |
| FastDup | 1.0.0, commit `7b3f62587283a257fd38c5c20ceeda0364c285ef`, with `patches/FastDup-local.patch` |
| bcftools / htslib | 1.10.2 |
| Ensembl VEP | 115.2, offline cache release 115 for `mus_musculus`, GRCm39 |
| PLINK 2 | 2.0.0-a.7LM |
| R | 4.x with base graphics; no extra package required by included plotting scripts |

Reference asset: GRCm39 RefSeq assembly `GCF_000001635.27`. The production FASTA SHA-256 is recorded in `docs/workflow.md`. Obtain the exact reference FASTA and VEP cache from their upstream distributors and verify checksums before use. Neither is bundled here.

FastDup is built from its upstream source commit with the patch in `patches/`. The sanitized build toolchain and the SHA-256 of the production binary are recorded in `patches/fastdup_build_provenance.json`; the binary itself is not distributed.

The files in `legacy/v1/envs_legacy/` describe an older workflow and should not be used to install the current analysis environment.
