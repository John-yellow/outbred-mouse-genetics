#!/usr/bin/env python3
"""Summarize main-set AF by exact VEP Existing_variation rs match."""
import argparse
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

AUTOSOMES = {str(i) for i in range(1, 20)}
X_NONPAR = ((1, 168752754), (169376593, 169476592))
BINS = [0.0, 0.01, 0.05] + [i / 100 for i in range(10, 101, 5)]
CLASSES = ("dbSNP_rs_match", "candidate_novel_vs_dbSNP", "unclassified_no_CSQ")
DBSNP_ID = re.compile(r"(^|[&,;])rs[0-9]+($|[&,;])")


def bin_labels():
    result=[]
    for i,(lo,hi) in enumerate(zip(BINS[:-1],BINS[1:])):
        result.append("[%.2f, %.2f%s"%(lo,hi,"]" if i==len(BINS)-2 else ")"))
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--bcf", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--plot-script", type=Path, default=Path(__file__).with_name("plot_known_novel_spectrum.R"))
    args=ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    stats=args.output_dir/"stats"; plots=args.output_dir/"plots"
    stats.mkdir(parents=True,exist_ok=True); plots.mkdir(parents=True,exist_ok=True)
    header=subprocess.run([args.bcftools,"view","-h",str(args.bcf)],check=True,text=True,stdout=subprocess.PIPE).stdout
    csq_line=next((x for x in header.splitlines() if x.startswith("##INFO=<ID=CSQ,")),None)
    if not csq_line:
        raise SystemExit("Input header has no VEP CSQ declaration")
    match=re.search(r"Format: ([^\"]+)",csq_line)
    if not match:
        raise SystemExit("Could not parse CSQ Format from VCF header")
    csq_fields=match.group(1).split("|")
    if "Existing_variation" not in csq_fields:
        raise SystemExit("CSQ header has no Existing_variation field")
    existing_i=csq_fields.index("Existing_variation")
    fmt=r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\n"
    query=subprocess.Popen([args.bcftools,"query","-f",fmt,str(args.bcf)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    counts=Counter(); ac1=Counter(); hist=Counter(); denominators=Counter(); total=0
    for line in query.stdout:
        chrom,pos_s,ref,alt,ac_s,an_s,csq=line.rstrip("\n").split("\t")
        pos=int(pos_s); ac=int(ac_s.split(",")[0]); an=int(an_s)
        if an<=0 or not (0<ac<an): raise RuntimeError(f"Invalid AC/AN at {chrom}:{pos}: {ac}/{an}")
        if chrom in AUTOSOMES: genomic="autosome"
        elif chrom=="X" and any(lo<=pos<=hi for lo,hi in X_NONPAR): genomic="X_nonPAR"
        else: raise RuntimeError(f"Unexpected site outside autosome/X non-PAR scope: {chrom}:{pos}")
        seen_rs=False
        if csq not in ("", "."):
            for row in csq.split(","):
                fields=row.split("|")
                if len(fields)>existing_i and DBSNP_ID.search(fields[existing_i]):
                    seen_rs=True; break
            variant="dbSNP_rs_match" if seen_rs else "candidate_novel_vs_dbSNP"
        else:
            variant="unclassified_no_CSQ"
        af=ac/an
        bi=next((i for i,(lo,hi) in enumerate(zip(BINS[:-1],BINS[1:])) if lo<=af<hi),None)
        if af==1.0: bi=len(BINS)-2
        if bi is None: raise RuntimeError(f"AF outside bins: {af}")
        key=(variant,genomic)
        counts[key]+=1; hist[(variant,genomic,bi)]+=1; denominators[key]+=1
        ac1[key]+=int(ac==1); total+=1
    query.stdout.close(); err=query.stderr.read(); rc=query.wait()
    if rc: raise RuntimeError("bcftools query failed: "+err)
    labels=bin_labels()
    with (stats/"known_candidate_novel_classification_counts.tsv").open("w",newline="") as handle:
        w=csv.writer(handle,delimiter="\t",lineterminator="\n")
        w.writerow(["variant_class","genomic_class","variant_count","within_class_denominator","proportion_of_genomic_class","AC1_count","classification_database_source"])
        for variant in CLASSES:
            for genomic in ("autosome","X_nonPAR"):
                n=counts[(variant,genomic)]; denom=sum(counts[(v,genomic)] for v in CLASSES)
                w.writerow([variant,genomic,n,denom,("%.10g"%(n/denom) if denom else "NA"),ac1[(variant,genomic)],"VEP Existing_variation; dbSNP 150"])
    with (stats/"known_candidate_novel_allele_frequency_spectrum.tsv").open("w",newline="") as handle:
        w=csv.writer(handle,delimiter="\t",lineterminator="\n")
        w.writerow(["variant_class","genomic_class","AF_bin","variant_count","within_class_denominator","proportion_within_class","AC1_total_in_class","frequency_definition","classification_database_source"])
        for variant in CLASSES:
            for genomic in ("autosome","X_nonPAR"):
                denom=counts[(variant,genomic)]
                for bi,label in enumerate(labels):
                    n=hist[(variant,genomic,bi)]
                    w.writerow([variant,genomic,label,n,denom,("%.10g"%(n/denom) if denom else "NA"),ac1[(variant,genomic)],"AF=AC/AN from effective GTs","VEP Existing_variation; dbSNP 150"])
    expected_classes=sum(counts.values())
    if expected_classes!=total: raise RuntimeError("Internal classification count mismatch")
    subprocess.run([args.rscript,str(args.plot_script),str(args.output_dir)],check=True)
    print(f"COMPLETE main_records={total} output={stats/'known_candidate_novel_allele_frequency_spectrum.tsv'}")


if __name__=="__main__":
    main()
