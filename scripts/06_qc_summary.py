#!/usr/bin/env python3
import csv, gzip, json, math, os, re, statistics, subprocess, sys, time, shutil
from bisect import bisect_right
from datetime import datetime
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else REPO/'work'
CFG=json.loads((BASE/'state/run_config.json').read_text())
def cfg_path(name, default):
    p=Path(CFG.get(name, default))
    return p if p.is_absolute() else BASE/p
MAIN=cfg_path('main_bcf','results/main_variants_vep115_GRCm39.bcf')
MAN=cfg_path('sample_manifest','inputs/sample_manifest.tsv')
PLID=cfg_path('ploidy_tsv','inputs/ploidy.tsv')
FAI=cfg_path('reference_fai','inputs/reference.fa.fai')
BCFTOOLS=shutil.which(CFG.get('bcftools','bcftools')) or CFG.get('bcftools','bcftools')
PLINK2=shutil.which(CFG.get('plink2','plink2')) or CFG.get('plink2','plink2')
RSCRIPT=shutil.which(CFG.get('rscript','Rscript')) or CFG.get('rscript','Rscript')
PROGRESS=BASE/'state/qc_progress.json'
COMMANDS=BASE/'state/commands.jsonl'
LOG=BASE/'logs/qc/qc_master.log'
LOG.parent.mkdir(parents=True,exist_ok=True)
SAMPLES={}; GROUPS={}; SEXES={}
with MAN.open() as f:
    reader=csv.DictReader(f,delimiter='\t')
    for row in reader:
        SAMPLES[row['sample_id']]=row
        GROUPS[row['sample_id']]=row['group']
        SEXES[row['sample_id']]=row['sex']
SAMPLE_ORDER=list(SAMPLES)
AUTOSOMES=[str(i) for i in range(1,20)]
WIDTH=1000000
START=time.time()

def now():
    return datetime.now().astimezone().isoformat()

def atomic_json(path,obj):
    tmp=Path(str(path)+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
    os.replace(str(tmp),str(path))

def write_qc_acceptance_summary():
    qcp=BASE/'state/qc_progress.json'; annp=BASE/'state/annotation_pipeline_summary.json'
    sitep=BASE/'state/site_metrics_summary.json'; plinkp=BASE/'state/plink_qc_summary.json'
    qc_complete=(sitep.exists() and plinkp.exists() and
                 json.loads(sitep.read_text()).get('status')=='COMPLETE' and
                 json.loads(plinkp.read_text()).get('status')=='COMPLETE')
    ann_complete=annp.exists() and json.loads(annp.read_text()).get('status')=='COMPLETE'
    if qc_complete and ann_complete: state_name='COMPLETE'
    elif qc_complete: state_name='QC_COMPLETE_ANNOTATION_PENDING'
    elif ann_complete: state_name='ANNOTATION_COMPLETE_QC_PENDING'
    else: state_name='RUNNING'
    atomic_json(BASE/'state/qc_acceptance_summary.json',{
        'status':state_name,'updated_at':now(),'qc_complete':qc_complete,'annotation_complete':ann_complete,
        'qc_progress':str(qcp),'site_metrics_summary':str(sitep),'plink_qc_summary':str(plinkp),
        'annotation_pipeline_summary':str(annp),
        'annotation_records':json.loads(annp.read_text()).get('candidate_records') if ann_complete else None,
        'main_annotation_records':json.loads(annp.read_text()).get('main_records') if ann_complete else None})

def proc_start_ticks(pid):
    try:
        txt=Path('/proc/%s/stat'%pid).read_text()
        tail=txt[txt.rfind(')')+2:].split()
        return int(tail[19])
    except Exception:
        return None

def status(stage, state='RUNNING', **kw):
    obj={'status':state,'stage':stage,'updated_at':now(),'elapsed_seconds':round(time.time()-START,2),
         'main_bcf':str(MAIN),'expected_main_records':CFG['expected_main_records'],
         'sample_count':CFG['sample_count']}
    obj.update(kw); atomic_json(PROGRESS,obj)

def record_cmd(stage,cmd,log_path,pid,start_ticks,started,finished,code,**kw):
    row={'stage':stage,'command':cmd,'log':str(log_path),'pid':pid,'pid_start_ticks':start_ticks,
         'started_at':started,'finished_at':finished,'exit_code':code}
    row.update(kw)
    with COMMANDS.open('a') as f:
        f.write(json.dumps(row,ensure_ascii=False)+'\n')
    return row

def run_cmd(stage,cmd,log_path):
    started=now(); log_path=Path(log_path); log_path.parent.mkdir(parents=True,exist_ok=True)
    with log_path.open('w') as log:
        log.write('COMMAND_JSON='+json.dumps(cmd,ensure_ascii=False)+'\n')
        log.flush()
        p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=False)
        st=proc_start_ticks(p.pid)
        status(stage,child_pid=p.pid,child_start_ticks=st,command=cmd,log=str(log_path),command_started_at=started)
        code=p.wait()
    finished=now(); row=record_cmd(stage,cmd,log_path,p.pid,st,started,finished,code)
    if code!=0: raise RuntimeError('%s failed with exit %d; see %s'%(stage,code,log_path))
    return row

def parse_fai(path):
    ans={}
    for line in path.read_text().splitlines():
        if line:
            p=line.split('\t'); ans[p[0]]=int(p[1])
    return ans

def x_intervals():
    out=[]
    for line in PLID.read_text().splitlines():
        p=line.split()
        if len(p)>=5 and p[0]=='X' and p[1]!='*' and p[2]!='*' and p[3]=='M' and p[4]=='1':
            out.append((int(p[1]),int(p[2])))
    return out

def windows_for(chrom,length,intervals=None):
    out={}
    for st in range(1,length+1,WIDTH):
        en=min(length,st+WIDTH-1)
        bp=en-st+1
        if intervals is not None:
            bp=sum(max(0,min(en,b)-max(st,a)+1) for a,b in intervals)
            if bp==0: continue
        out[(chrom,st)]={'chrom':chrom,'start':st,'end':en,'bp':bp,'n':0}
    return out

def write_tsv(path,header,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(header); w.writerows(rows)

def make_metrics():
    lengths=parse_fai(FAI)
    intervals=x_intervals()
    if not intervals: raise RuntimeError('No male haploid X intervals in ploidy.tsv')
    windows={}
    for c in AUTOSOMES:
        if c not in lengths: raise RuntimeError('Missing autosome %s in reference FAI'%c)
        windows.update(windows_for(c,lengths[c]))
    if 'X' not in lengths: raise RuntimeError('Missing X in reference FAI')
    windows.update(windows_for('X',lengths['X'],intervals))
    stats={'main_records':0,'autosomal_records':0,'X_records':0,'X_nonPAR_records':0,'X_PAR_records':0,
           'unexpected_contig_records':0,'ac1_autosome':0,'ac1_X_nonPAR':0,'ac1_X_PAR':0}
    bins=[0.0,0.01,0.05]+[i/100.0 for i in range(10,101,5)]
    labels=[]
    for i in range(len(bins)-1):
        labels.append('[%.2f, %.2f%s'%(bins[i],bins[i+1],']' if i==len(bins)-2 else ')'))
    classes=['autosome','X_nonPAR','X_PAR']
    hist={c:[0]*(len(bins)-1) for c in classes}
    totals={c:0 for c in classes}; ac1={c:0 for c in classes}
    base_out=BASE/'stats'; base_out.mkdir(exist_ok=True)
    site_path=BASE/'stats/main_snp_site_metrics.tsv.gz'
    fmt=r'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n'
    cmd=[BCFTOOLS,'query','-f',fmt,str(MAIN)]
    log_path=BASE/'logs/qc/main_site_metrics.stderr.log'
    started=now()
    err=log_path.open('w')
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1,start_new_session=False)
    child_start=proc_start_ticks(p.pid)
    status('MAIN_SITE_METRICS',child_pid=p.pid,child_start_ticks=child_start,command=cmd,log=str(log_path),command_started_at=started,records_seen=0)
    seen=0; last_flush=0
    with gzip.open(site_path,'wt',compresslevel=5) as gz:
        gz.write('CHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n')
        for line in p.stdout:
            f=line.rstrip('\n').split('\t')
            if len(f)!=6: raise RuntimeError('Unexpected bcftools query row width %d at record %d'%(len(f),seen+1))
            chrom,pos,ref,alt,ac_s,an_s=f
            ac=int(ac_s.split(',')[0]); an=int(an_s); pos=int(pos)
            if an<=0 or not (0<ac<an): raise RuntimeError('Invalid AC/AN at %s:%d: %s/%s'%(chrom,pos,ac_s,an_s))
            af=ac/an; seen+=1; stats['main_records']+=1
            gz.write('%s\t%d\t%s\t%s\t%d\t%d\t%.10g\n'%(chrom,pos,ref,alt,ac,an,af))
            if chrom in AUTOSOMES:
                stats['autosomal_records']+=1; cls='autosome'
                key=(chrom,((pos-1)//WIDTH)*WIDTH+1)
                if key in windows: windows[key]['n']+=1
            elif chrom=='X':
                stats['X_records']+=1
                nonpar=any(a<=pos<=b for a,b in intervals)
                if nonpar:
                    cls='X_nonPAR'; stats['X_nonPAR_records']+=1
                    key=('X',((pos-1)//WIDTH)*WIDTH+1)
                    if key in windows: windows[key]['n']+=1
                else:
                    cls='X_PAR'; stats['X_PAR_records']+=1
            else:
                stats['unexpected_contig_records']+=1
                continue
            totals[cls]+=1
            ix=bisect_right(bins,af)-1
            if af==1.0: ix=len(bins)-2
            if ix<0 or ix>=len(labels): raise RuntimeError('AF outside bins: %s'%af)
            hist[cls][ix]+=1
            if ac==1: ac1[cls]+=1
            if seen-last_flush>=250000:
                last_flush=seen
                status('MAIN_SITE_METRICS',child_pid=p.pid,child_start_ticks=child_start,command=cmd,log=str(log_path),
                       records_seen=seen,expected_records=CFG['expected_main_records'],fraction=seen/CFG['expected_main_records'])
    p.stdout.close(); err.close(); code=p.wait(); finished=now()
    record_cmd('MAIN_SITE_METRICS',cmd,log_path,p.pid,child_start,started,finished,code,records_seen=seen,output=str(site_path))
    if code!=0: raise RuntimeError('bcftools query failed (%d); see %s'%(code,log_path))
    if seen!=CFG['expected_main_records']: raise RuntimeError('Main record count mismatch: %d != %d'%(seen,CFG['expected_main_records']))
    if stats['unexpected_contig_records']: raise RuntimeError('Unexpected contig records present: %d'%stats['unexpected_contig_records'])
    window_rows=[]
    by_chrom={}
    for (chrom,st),x in sorted(windows.items(),key=lambda z:((int(z[0][0]) if z[0][0].isdigit() else 23),z[0][1])):
        rate=x['n']*1e6/x['bp']
        window_rows.append([chrom,st,x['end'],x['bp'],x['n'],'%.8f'%rate])
        by_chrom.setdefault(chrom,[]).append((x['n'],x['bp'],rate))
    write_tsv(base_out/'window_variant_density_1Mb.tsv',
              ['chrom','window_start_1based','window_end_1based_inclusive','included_reference_bp','snp_count','snps_per_Mb'],window_rows)
    summary=[]
    for chrom,vals in by_chrom.items():
        counts=[v[0] for v in vals]; rates=[v[2] for v in vals]; bpsum=sum(v[1] for v in vals)
        summary.append([chrom,len(vals),sum(v==0 for v in counts),sum(counts),bpsum,
                        '%.8f'%statistics.mean(rates),'%.8f'%statistics.median(rates),
                        '%.8f'%min(rates),'%.8f'%max(rates)])
    autos=[v for c,vs in by_chrom.items() if c in AUTOSOMES for v in vs]
    auto_rates=[v[2] for v in autos]; auto_counts=[v[0] for v in autos]
    summary.append(['AUTOSOME_1_19',len(autos),sum(v==0 for v in auto_counts),sum(auto_counts),sum(v[1] for v in autos),
                    '%.8f'%statistics.mean(auto_rates),'%.8f'%statistics.median(auto_rates),
                    '%.8f'%min(auto_rates),'%.8f'%max(auto_rates)])
    write_tsv(base_out/'window_density_summary.tsv',
              ['chrom_or_scope','window_count','zero_snp_windows','snp_count','included_reference_bp','mean_window_snps_per_Mb','median_window_snps_per_Mb','min_window_snps_per_Mb','max_window_snps_per_Mb'],summary)
    freq_rows=[]
    outclasses=['autosome','X_nonPAR','X_PAR','X_total']
    for cls in outclasses:
        denom=totals['X_nonPAR']+totals['X_PAR'] if cls=='X_total' else totals.get(cls,0)
        ac1_n=ac1['X_nonPAR']+ac1['X_PAR'] if cls=='X_total' else ac1.get(cls,0)
        for i,label in enumerate(labels):
            n=(hist['X_nonPAR'][i]+hist['X_PAR'][i]) if cls=='X_total' else hist.get(cls,[0]*len(labels))[i]
            freq_rows.append(['unclassified_pending_allele_database',cls,label,n,denom,('%.10g'%(n/denom) if denom else 'NA'),ac1_n,'AF=AC/AN'])
    write_tsv(base_out/'allele_frequency_spectrum.tsv',
              ['variant_class','genomic_class','AF_bin','variant_count','within_class_denominator','proportion_within_class','AC1_total_in_class','frequency_definition'],freq_rows)
    write_tsv(base_out/'ac1_site_counts.tsv',['genomic_class','AC1_variant_count','main_variant_count','proportion'],
              [[c,ac1.get(c,ac1['X_nonPAR']+ac1['X_PAR'] if c=='X_total' else 0),
                totals.get(c,totals['X_nonPAR']+totals['X_PAR'] if c=='X_total' else 0),
                ('%.10g'%((ac1.get(c,ac1['X_nonPAR']+ac1['X_PAR'] if c=='X_total' else 0))/(totals.get(c,totals['X_nonPAR']+totals['X_PAR'] if c=='X_total' else 0))) if totals.get(c,totals['X_nonPAR']+totals['X_PAR'] if c=='X_total' else 0) else 'NA')]
               for c in outclasses])
    atomic_json(BASE/'state/site_metrics_summary.json',{'status':'COMPLETE','completed_at':now(),'input':str(MAIN),
        'records':stats,'frequency_bins':bins,'frequency_classification':'genomic classes; known/candidate-novel class is summarized by the separate annotation spectrum script',
        'ac1_counts':ac1,'genomic_class_counts':totals,'x_nonpar_intervals':intervals,'query_command':cmd,
        'site_metrics_output':str(site_path)})
    status('MAIN_SITE_METRICS',child_pid=p.pid,child_start_ticks=child_start,records_seen=seen,expected_records=seen,output=str(site_path),stage_exit_code=0)
    return stats,intervals,by_chrom

def parse_het(path):
    with Path(path).open() as f:
        header=None
        for line in f:
            if line.startswith('#') or line.strip():
                header=line.strip().lstrip('#').split()
                break
        if not header: raise RuntimeError('Empty PLINK2 .het: '+str(path))
        upper=[x.upper() for x in header]
        iid_idx=upper.index('IID') if 'IID' in upper else None
        het_candidates=['HET','OBS_HET','HET_OBS','OBSERVED_HET','O(HET)']
        het_idx=next((upper.index(x) for x in het_candidates if x in upper),None)
        nobs_candidates=['OBS_CT','NOBS','N_OBS','NONMISSING']
        nobs_idx=next((upper.index(x) for x in nobs_candidates if x in upper),None)
        if iid_idx is None or het_idx is None or nobs_idx is None:
            raise RuntimeError('Could not resolve raw HET/NOBS columns in %s header=%s'%(path,header))
        out={}
        for line in f:
            if not line.strip(): continue
            row=line.split()
            iid=row[iid_idx]; out[iid]=(float(row[het_idx]),float(row[nobs_idx]))
        return header,out

def plink2_version_str(plink):
    try:
        res = subprocess.run([plink, '--version'], check=True, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        first = res.stdout.strip().splitlines()
        if first:
            return first[0].strip()
    except Exception:
        pass
    return 'PLINK v2.0.0-a.7LM (5 Dec 2025)'


def run_plink(stats, resume_after_pgen=False, resume_after_het=False):
    plink=PLINK2; main=MAIN
    plink_threads=str(CFG.get('plink_threads', 16))
    pdir=BASE/'plink'; pdir.mkdir(exist_ok=True)
    pgen=pdir/'main_autosomes'
    pvar=Path(str(pgen)+'.pvar')
    if resume_after_pgen:
        if not all(Path(str(pgen)+ext).is_file() for ext in ['.pgen','.pvar','.psam']):
            raise RuntimeError('Cannot resume after PGEN: one or more PLINK files are missing')
        status('PLINK_AUTOSOME_PGEN_REUSED',reused_existing_pgen=True,pgen_prefix=str(pgen))
    else:
        cmd=[plink,'--bcf',str(main),'--chr','1-19','--set-all-var-ids','@:#:$r:$a','--make-pgen','--threads',plink_threads,'--out',str(pgen)]
        run_cmd('PLINK_AUTOSOME_PGEN',cmd,BASE/'logs/qc/plink_autosome_conversion.log')
    pvar_n=sum(1 for line in pvar.open() if line.strip() and not line.startswith('#'))
    if pvar_n!=stats['autosomal_records']: raise RuntimeError('PLINK autosome variant count mismatch %d vs %d'%(pvar_n,stats['autosomal_records']))
    hetprefix=pdir/'main_autosomes_het'
    hetfile=Path(str(hetprefix)+'.het')
    if resume_after_het:
        if not hetfile.is_file(): raise RuntimeError('Cannot resume after HET: PLINK .het output is missing')
        status('PLINK_HET_REUSED',reused_existing_het=True,het_file=str(hetfile))
    else:
        cmd=[plink,'--pfile',str(pgen),'--het','cols=maybefid,het,nobs','--threads',plink_threads,'--out',str(hetprefix)]
        run_cmd('PLINK_HET',cmd,BASE/'logs/qc/plink_het.log')
    het_header,values=parse_het(hetfile)
    if set(values)!=set(SAMPLE_ORDER): raise RuntimeError('PLINK HET sample set mismatch')
    het_rows=[]
    for sid in SAMPLE_ORDER:
        h,n=values[sid]; n=int(round(n)); h=int(round(h)); missing=stats['autosomal_records']-n
        if n<0 or n>stats['autosomal_records'] or h<0 or h>n: raise RuntimeError('Invalid per-sample PLINK het values '+sid)
        het_rows.append([sid,GROUPS[sid],SEXES[sid],stats['autosomal_records'],h,n,missing,
                         '%.10g'%(missing/stats['autosomal_records'] if stats['autosomal_records'] else float('nan')),
                         '%.10g'%(h/n if n else float('nan'))])
    write_tsv(BASE/'stats/heterozygosity_per_sample.tsv',
        ['sample_id','group','sex','autosomal_SNP_count','heterozygous_GT_count','nonmissing_diploid_GT_count','missing_GT_count','missing_rate','heterozygous_proportion_of_nonmissing_diploid_GT'],het_rows)
    group_rows=[]
    for g in ['ICR1','ICR2','KM']:
        vals=[float(x[-1]) for x in het_rows if x[1]==g]
        group_rows.append([g,len(vals),'%.10g'%statistics.mean(vals),'%.10g'%statistics.median(vals),'%.10g'%min(vals),'%.10g'%max(vals)])
    write_tsv(BASE/'stats/heterozygosity_group_summary.tsv',['group','sample_count','mean_heterozygous_proportion','median_heterozygous_proportion','min','max'],group_rows)
    prune=pdir/'main_autosomes_pruned'
    cmd=[plink,'--pfile',str(pgen),'--indep-pairwise','100','10','0.2','--threads',plink_threads,'--out',str(prune)]
    run_cmd('PLINK_LD_PRUNE',cmd,BASE/'logs/qc/plink_ld_prune.log')
    prune_in=Path(str(prune)+'.prune.in')
    nprune=sum(1 for x in prune_in.open() if x.strip())
    if nprune<=0: raise RuntimeError('Empty LD-pruned variant list')
    king=pdir/'main_autosomes_king_pruned'
    cmd=[plink,'--pfile',str(pgen),'--extract',str(prune_in),'--make-king-table','--threads',plink_threads,'--out',str(king)]
    run_cmd('PLINK_KING_PRUNED',cmd,BASE/'logs/qc/plink_king_pruned.log')
    kinfile=Path(str(king)+'.kin0')
    raw=[]
    with kinfile.open() as f:
        kin0_header=None
        for line in f:
            if line.startswith('#') or not line.strip():
                if line.startswith('#'): kin0_header=line.strip().lstrip('#').split()
                continue
            raw.append(line.split())
    if not kin0_header: raise RuntimeError('Missing .kin0 header')
    hidx={k.upper():i for i,k in enumerate(kin0_header)}
    def findcol(*names):
        for n in names:
            if n in hidx: return hidx[n]
        raise RuntimeError('Missing columns %s in %s'%(names,kin0_header))
    i1=findcol('IID1','ID1'); i2=findcol('IID2','ID2'); ik=findcol('KINSHIP','KINSHIP_COEFFICIENT')
    pairs={}
    for row in raw:
        a=row[i1]; b=row[i2]
        if a not in GROUPS or b not in GROUPS: raise RuntimeError('Unexpected IDs in KING output: %s %s'%(a,b))
        if GROUPS[a]!=GROUPS[b]: continue
        key=tuple(sorted((a,b)))
        if key in pairs: raise RuntimeError('Duplicate KING pair '+str(key))
        k=float(row[ik])
        if not math.isfinite(k): raise RuntimeError('Nonfinite KING coefficient '+str(key))
        pairs[key]=(GROUPS[a],k)
    expected={g:30*29//2 for g in ['ICR1','ICR2','KM']}
    got={g:sum(1 for v in pairs.values() if v[0]==g) for g in expected}
    if got!=expected: raise RuntimeError('Expected all 435 pairs per group; got '+json.dumps(got))
    def level(k):
        if k>=0.354: return 'duplicate_or_MZ_like>=0.354'
        if k>=0.177: return 'first_degree_heuristic_0.177_0.354'
        if k>=0.0884: return 'second_degree_heuristic_0.0884_0.177'
        if k>=0.0442: return 'third_degree_heuristic_0.0442_0.0884'
        return 'below_0.0442'
    pairrows=[]
    for (a,b),(g,k) in sorted(pairs.items(),key=lambda x:(x[1][0],-x[1][1],x[0])):
        pairrows.append([g,a,b,'%.10g'%k,level(k),nprune,'KING-robust; heuristic bins; pedigree required for biological relationship'])
    write_tsv(BASE/'stats/within_group_kinship_1305_pairs.tsv',
       ['group','sample1','sample2','king_robust_kinship','screening_category','LD_pruned_SNP_count','interpretation_note'],pairrows)
    summary=[]
    for g in expected:
        gs=[(a,b,k) for (a,b),(gg,k) in pairs.items() if gg==g]
        summary.append([g,len(gs),sum(k>=0.354 for a,b,k in gs),sum(k>=0.177 for a,b,k in gs),
                       sum(k>=0.0884 for a,b,k in gs),sum(k>=0.0442 for a,b,k in gs),
                       '%.10g'%max(k for a,b,k in gs)])
    write_tsv(BASE/'stats/within_group_kinship_summary.tsv',
       ['group','pair_count','pairs_ge_0.354_duplicate_screen','pairs_ge_0.177_first_degree_screen','pairs_ge_0.0884_second_degree_screen','pairs_ge_0.0442_third_degree_screen','max_kinship'],summary)
    atomic_json(BASE/'state/plink_qc_summary.json',{'status':'COMPLETE','plink2_version':plink2_version_str(plink),
        'input':str(main),'autosome_scope':'1-19','autosomal_SNP_count':stats['autosomal_records'],
        'LD_pruning_command':['plink2','--pfile',str(pgen),'--indep-pairwise','100','10','0.2'],
        'LD_pruned_SNP_count':nprune,'kinship_method':'PLINK2 KING-robust --make-king-table',
        'kinship_input_variant_list':str(prune_in),'within_group_pairs':got,
        'kinship_thresholds':{'0.354':'duplicate/MZ-like screen','0.177':'first-degree screen','0.0884':'second-degree screen','0.0442':'third-degree screen'},
        'interpretation':'screening only; cannot replace pedigree, sibship, or generation records',
        'sample_heterozygosity_file':str(BASE/'stats/heterozygosity_per_sample.tsv'),
        'het_header':het_header,'kin0_header':kin0_header,'het_rows':len(het_rows),'kinship_rows':len(pairrows)})
    return nprune,pairrows

def write_x_ploidy():
    cmd=[BCFTOOLS,'query','-r','X','-f',r'%CHROM\t%POS[\t%GT]\n',str(MAIN)]
    log=BASE/'logs/qc/x_ploidy_query.stderr.log'; log.parent.mkdir(parents=True,exist_ok=True)
    counts={'records_checked':0,'male_haploid_called':0,'male_missing':0,
            'female_diploid_called':0,'female_missing':0,'male_heterozygous':0}
    samples_cmd=[BCFTOOLS,'query','-l',str(MAIN)]
    sample_proc=subprocess.run(samples_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=True)
    if sample_proc.stdout.splitlines()!=SAMPLE_ORDER:
        raise RuntimeError('X-ploidy input sample order differs from manifest')
    with log.open('w') as err:
        p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
        for line in p.stdout:
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=len(SAMPLE_ORDER)+2: raise RuntimeError('Unexpected X genotype row width')
            pos=int(fields[1])
            if any(a<=pos<=b for a,b in ((168752755,169376592),)):
                raise RuntimeError('PAR site present in X non-PAR callset at %d'%pos)
            counts['records_checked']+=1
            for sid,sex,gt in zip(SAMPLE_ORDER,[SEXES[x] for x in SAMPLE_ORDER],fields[2:]):
                missing=gt in ('.','./.','.|.') or '.' in gt.replace('|','/').split('/')
                alleles=gt.replace('|','/').split('/')
                if sex=='M':
                    if missing:
                        counts['male_missing']+=1
                    else:
                        if len(alleles)!=1: raise RuntimeError('Male X genotype is not haploid for '+sid+': '+gt)
                        counts['male_haploid_called']+=1
                        if len(set(alleles))>1: counts['male_heterozygous']+=1
                else:
                    if missing:
                        counts['female_missing']+=1
                    else:
                        if len(alleles)!=2: raise RuntimeError('Female X genotype is not diploid for '+sid+': '+gt)
                        counts['female_diploid_called']+=1
        p.stdout.close(); code=p.wait()
    if code: raise RuntimeError('bcftools X-ploidy query failed: %d'%code)
    if counts['male_heterozygous']:
        raise RuntimeError('Male haploid X GTs contained heterozygous calls')
    out=BASE/'stats/x_ploidy_validation_summary.tsv'
    intervals=x_intervals()
    write_tsv(out,['metric','value','source_or_definition'],[
        ['X_records_checked',counts['records_checked'],'all X records in the main set'],
        ['male_called_haploid_genotypes',counts['male_haploid_called'],'nonmissing male X non-PAR GTs'],
        ['male_missing_genotypes',counts['male_missing'],'missing male X GTs'],
        ['female_called_diploid_genotypes',counts['female_diploid_called'],'nonmissing female X non-PAR GTs'],
        ['female_missing_genotypes',counts['female_missing'],'missing female X GTs'],
        ['male_heterozygous_genotypes',counts['male_heterozygous'],'must equal zero for haploid calls'],
        ['male_haploid_X_intervals_1based_inclusive',','.join('%d-%d'%z for z in intervals),'ploidy manifest'],
        ['checked_scope','X non-PAR records in main set','PAR is excluded from the callset']])
    return counts

def finish_plots():
    (BASE/'plots').mkdir(parents=True,exist_ok=True)
    cmd=[RSCRIPT,str(REPO/'scripts/make_plots.R'),str(BASE)]
    run_cmd('PLOTS',cmd,BASE/'logs/qc/make_plots.log')

def main():
    COMMANDS.touch(exist_ok=True)
    resume_after_het='--resume-after-het' in sys.argv[1:]
    resume_after_pgen='--resume-after-pgen' in sys.argv[1:] or resume_after_het
    if resume_after_pgen:
        summary=json.loads((BASE/'state/site_metrics_summary.json').read_text())
        if summary.get('status')!='COMPLETE': raise RuntimeError('Existing site metrics are not complete; cannot resume QC tail')
        stats=summary['records']
        if stats['main_records']!=CFG['expected_main_records'] or stats['autosomal_records']<=0:
            raise RuntimeError('Existing site metrics do not match the frozen run configuration')
        status('RESUME_AFTER_PGEN',input_identity_frozen=True,reused_site_metrics=True,
               reused_site_metrics_file=str(BASE/'stats/main_snp_site_metrics.tsv.gz'))
        run_plink(stats,resume_after_pgen=True,resume_after_het=resume_after_het)
    else:
        status('STARTUP',input_identity_frozen=True)
        stats,intervals,by_chrom=make_metrics()
        write_x_ploidy()
        run_plink(stats)
    finish_plots()
    annotation=json.loads((BASE/'state/annotation_progress.json').read_text()) if (BASE/'state/annotation_progress.json').exists() else {'status':'NOT_RUN'}
    state='QC_COMPLETE_ANNOTATION_PENDING' if annotation.get('status') not in ('COMPLETE','NOT_RUN') else 'QC_COMPLETE'
    final={'status':state,'stage':'INDEPENDENT_QC_COMPLETE','updated_at':now(),'elapsed_seconds':round(time.time()-START,2),
           'main_site_metrics':stats,'annotation_classification_pending':annotation.get('status')!='COMPLETE',
           'outputs':[str(x) for x in sorted((BASE/'stats').glob('*')) if x.is_file()]}
    atomic_json(PROGRESS,final)
    write_qc_acceptance_summary()
    print(json.dumps(final,ensure_ascii=False))
    return 0

if __name__=='__main__':
    try:
        rc=main()
        sys.exit(rc)
    except Exception as e:
        with LOG.open('a') as f: f.write('%s FATAL %s\n'%(now(),repr(e)))
        status('FAILED','FAILED',error=repr(e))
        raise
