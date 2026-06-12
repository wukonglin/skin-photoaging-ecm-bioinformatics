#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20_public_data_extension_analysis.py  --  ROUND 3 public-data audit: NEW replicate-aware
in-vivo confirmation from an independent public mouse UVB photoaging dataset (GSE284483)
========================================================================================
Why: the current public-only Figure 1 (Panel B) re-analyzes GSE173385, which has ONE biological
sample per condition (exploratory, cell-level only). To answer the reviewers' obvious question --
"is the photoaging-driven ECM/collagen suppression reproducible with biological replicates?" --
this script re-analyzes an INDEPENDENT, REPLICATE-AWARE public dataset:

  GSE284483  (Mus musculus, bulk RNA-seq, in vivo dorsal skin; Public 2025; PMID 41233777)
  Groups: NoUVB (n=3)  |  UVB-Saline (n=3)  |  UVB-hCOL3A1 mRNA (n=3)
  We use ONLY the public normal-vs-photoaged contrast: NoUVB (control) vs UVB-Saline (photoaged).
  The UVB-hCOL3A1 *treatment* arm is EXCLUDED (out of scope; like the protective third arm of GSE173385).

Method (transparent, replicate-level; no cell-level p-value inflation):
  - raw counts -> CPM -> log2(CPM+1)
  - ECM/collagen module = the SAME mouse gene set used for GSE173385 (script 12), for consistency
  - per-sample module score = mean of z-scored log2CPM across module genes (n=3 vs n=3)
  - replicate-level Welch t-test (NoUVB vs UVB-Saline) + Cohen's d
  - per-gene log2FC (UVB-Saline / NoUVB) for the structural collagen / matrisome core

  conda run -n skin-ecm-bioinformatics python 04_code/20_public_data_extension_analysis.py

Inputs (public; downloaded into the project, originals untouched):
  05_data_raw/public_extension_GSE284483/GSE284483_raw_counts_All_Samples.txt.gz
Outputs:
  07_results/GSE284483_UVB_ECM_module.tsv         (per-sample module score + group test)
  07_results/GSE284483_UVB_collagen_genes.tsv     (per-gene log2FC for the matrisome core)
"""
from __future__ import annotations
import os, gzip
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT=Path(os.environ.get("PROJECT_ROOT",Path(__file__).resolve().parent.parent)).resolve()
RAW=ROOT/"05_data_raw"/"public_extension_GSE284483"/"GSE284483_raw_counts_All_Samples.txt.gz"
OUTR=ROOT/"07_results"

# SAME mouse ECM/collagen module as GSE173385 (04_code/12_scRNA_UVB_GSE173385.py)
MODULE=["Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col6a1","Fn1","Dcn","Lum",
        "Eln","Fbn1","Postn","Tnc","Lox","Mmp2","Mmp3","Mmp9","Timp1"]
# structural matrix-production core (for the per-gene panel; excludes the MMP/TIMP regulators)
CORE=["Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col6a1","Fn1","Dcn","Lum","Eln","Fbn1","Postn"]

CTRL=["NoUVB.1","NoUVB.2","NoUVB.3"]
PHOTO=["UVB-Saline.1","UVB-Saline.2","UVB-Saline.3"]   # photoaged; treatment arm excluded

# ---------- load ----------
df=pd.read_csv(RAW,sep="\t")
df["sym"]=df["Gene"].astype(str).str.split("_").str[-1]
counts=df[[ "Gene","sym"]+CTRL+PHOTO].copy()
mat=counts[CTRL+PHOTO].apply(pd.to_numeric,errors="coerce").fillna(0.0)
# collapse duplicate symbols (sum counts)
mat.index=counts["sym"].values
mat=mat.groupby(level=0).sum()
# CPM + log2
cpm=mat/ mat.sum(axis=0) * 1e6
logcpm=np.log2(cpm+1.0)

# ---------- module score (z-score per gene across the 6 samples, then average) ----------
present=[g for g in MODULE if g in logcpm.index]
missing=[g for g in MODULE if g not in logcpm.index]
sub=logcpm.loc[present]
z=sub.sub(sub.mean(axis=1),axis=0).div(sub.std(axis=1,ddof=1).replace(0,np.nan),axis=0).dropna(how="any")
modscore=z.mean(axis=0)   # per-sample module score
g_ctrl=modscore[CTRL].values; g_photo=modscore[PHOTO].values
t,p=stats.ttest_ind(g_ctrl,g_photo,equal_var=False)
# Cohen's d (pooled)
sp=np.sqrt(((len(g_ctrl)-1)*g_ctrl.var(ddof=1)+(len(g_photo)-1)*g_photo.var(ddof=1))/(len(g_ctrl)+len(g_photo)-2))
d=(g_photo.mean()-g_ctrl.mean())/sp if sp>0 else np.nan

rows=[{"sample":s,"group":"NoUVB(control)","module_score":round(modscore[s],4)} for s in CTRL]
rows+=[{"sample":s,"group":"UVB-Saline(photoaged)","module_score":round(modscore[s],4)} for s in PHOTO]
rows.append({"sample":"__GROUP_TEST__","group":"NoUVB vs UVB-Saline (Welch t, n=3 vs 3)",
             "module_score":f"mean_ctrl={g_ctrl.mean():.3f}; mean_photo={g_photo.mean():.3f}; "
                            f"delta={g_photo.mean()-g_ctrl.mean():+.3f}; t={t:.3f}; p={p:.4f}; cohens_d={d:.2f}; "
                            f"module_genes_used={len(z)}/{len(MODULE)} (missing: {','.join(missing) or 'none'})"})
pd.DataFrame(rows).to_csv(OUTR/"GSE284483_UVB_ECM_module.tsv",sep="\t",index=False)

# ---------- per-gene log2FC for the structural core ----------
def l2fc(g):
    if g not in cpm.index: return np.nan,np.nan,np.nan
    c=cpm.loc[g,CTRL].mean(); pn=cpm.loc[g,PHOTO].mean()
    return np.log2((pn+1)/(c+1)), c, pn
gene_rows=[]
for g in CORE:
    fc,c,pn=l2fc(g)
    gene_rows.append({"gene":g,"mean_CPM_NoUVB":round(c,2) if not np.isnan(c) else "NA",
                      "mean_CPM_UVB_Saline":round(pn,2) if not np.isnan(pn) else "NA",
                      "log2FC_UVBvsNoUVB":round(fc,3) if not np.isnan(fc) else "NA"})
pd.DataFrame(gene_rows).to_csv(OUTR/"GSE284483_UVB_collagen_genes.tsv",sep="\t",index=False)

print("=== GSE284483 (mouse UVB, in vivo, n=3/group) ECM/collagen module ===")
print(f"module genes used: {len(z)}/{len(MODULE)}  (missing: {missing or 'none'})")
print(f"NoUVB module mean   = {g_ctrl.mean():+.3f}  ({np.round(g_ctrl,3)})")
print(f"UVB-Saline mean     = {g_photo.mean():+.3f}  ({np.round(g_photo,3)})")
print(f"delta (photo-ctrl)  = {g_photo.mean()-g_ctrl.mean():+.3f}")
print(f"Welch t-test p      = {p:.4f}   Cohen's d = {d:.2f}")
print("per-gene log2FC (UVB-Saline / NoUVB):")
for r in gene_rows: print("  ",r["gene"],r["log2FC_UVBvsNoUVB"])
