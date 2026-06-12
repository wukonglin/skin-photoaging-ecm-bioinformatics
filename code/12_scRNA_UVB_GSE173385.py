#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12_scRNA_UVB_GSE173385.py  --  Panel B (UVB photoaging single-cell evidence)
============================================================================
Mouse UVB-irradiated skin scRNA-seq (GSE173385: Control / UV / UV+VitD).
Tissue-level test of whether UVB dysregulates the ECM / adhesion / collagen
program specifically in DERMAL FIBROBLASTS — the model of interest.

  conda run -n skin-ecm-bioinformatics python 04_code/12_scRNA_UVB_GSE173385.py

Input:  05_data_raw/GSE173385/GSM5266942_C5_matrix.tsv.gz (+ UV, VD)  [genes x cells]
Outputs:
  06_data_processed/GSE173385_annotated.h5ad
  07_results/GSE173385_fibroblast_module_scores.tsv
  08_figures/PanelB_scRNA_UVB.pdf/.png
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse, stats
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parent.parent)).resolve()
RAW = ROOT/"05_data_raw"/"GSE173385"
PROC= ROOT/"06_data_processed"; RES=ROOT/"07_results"; FIG=ROOT/"08_figures"
for d in (PROC,RES,FIG): d.mkdir(parents=True, exist_ok=True)

FILES = {"Control":"GSM5266942_C5_matrix.tsv.gz","UV":"GSM5266943_UV_matrix.tsv.gz","UV_VitD":"GSM5266944_VD_matrix.tsv.gz"}

# mouse axis gene sets
GS = {
 "ECM_collagen": ["Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col6a1","Fn1","Dcn","Lum","Eln","Fbn1","Postn","Tnc","Lox","Mmp2","Mmp3","Mmp9","Timp1"],
 "Adhesion_focal_adhesion": ["Itga1","Itga2","Itga5","Itgav","Itgb1","Ptk2","Src","Vcl","Pxn","Tln1","Actn1","Flna","Cav1","Ilk","Fermt2","Lims1","Lims2","Parva"],
 "Senescence_SASP_ROS": ["Cdkn1a","Cdkn2a","Trp53","Il6","Cxcl1","Serpine1","Gdf15","Apoe","Jun","Fos","Mapk14","Sod1","Sod2","Gpx1","Cat"],
}
FIB = ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"]
MARK = {"Fibroblast":FIB,"Keratinocyte":["Krt14","Krt5","Krt10"],"Immune":["Ptprc","Cd52","Lyz2"],
        "Endothelial":["Pecam1","Cdh5"],"Melanocyte":["Mlana","Pmel"],"Muscle":["Acta2","Myl9"]}

def load_one(cond, fn):
    df = pd.read_csv(RAW/fn, sep="\t", index_col=0)
    df = df[~df.index.duplicated(keep="first")]
    X = sparse.csr_matrix(df.values.T.astype(np.float32))     # cells x genes
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=df.columns.astype(str)),
                    var=pd.DataFrame(index=df.index.astype(str)))
    ad.obs["condition"] = cond
    ad.obs_names = [f"{cond}_{b}" for b in ad.obs_names]
    print(f"  {cond}: {ad.n_obs} cells x {ad.n_vars} genes")
    return ad

def main():
    for fn in FILES.values():
        if not (RAW/fn).exists(): sys.exit(f"[MISSING] {RAW/fn}")
    print("Loading 3 samples...")
    ad = sc.concat([load_one(c,f) for c,f in FILES.items()], join="outer", index_unique=None)
    ad.X = sparse.csr_matrix(np.nan_to_num(ad.X.toarray())) if not sparse.issparse(ad.X) else ad.X
    print(f"Combined: {ad.n_obs} cells x {ad.n_vars} genes")

    # QC
    ad.var["mt"] = ad.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(ad, min_genes=200); sc.pp.filter_genes(ad, min_cells=3)
    ad = ad[(ad.obs.n_genes_by_counts < 6000) & (ad.obs.pct_counts_mt < 25)].copy()
    print(f"After QC: {ad.n_obs} cells")

    ad.layers["counts"]=ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    ad.raw=ad
    sc.pp.highly_variable_genes(ad, n_top_genes=2000)
    adh=ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(adh, max_value=10); sc.tl.pca(adh, n_comps=30)
    sc.pp.neighbors(adh, n_neighbors=15, n_pcs=30); sc.tl.leiden(adh, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adh)
    ad.obs["leiden"]=adh.obs["leiden"]; ad.obsm["X_umap"]=adh.obsm["X_umap"]

    # cell-type score per cluster -> assign label by top marker module
    for ct,gs in MARK.items():
        sc.tl.score_genes(ad, [g for g in gs if g in ad.var_names], score_name=f"ct_{ct}")
    cl_scores = ad.obs.groupby("leiden")[[f"ct_{c}" for c in MARK]].mean()
    cl_label = {cl: cl_scores.loc[cl].idxmax().replace("ct_","") for cl in cl_scores.index}
    ad.obs["cell_type"]=ad.obs["leiden"].map(cl_label).astype("category")
    print("Cluster -> cell type:", cl_label)

    # axis module scores
    for name,gs in GS.items():
        sc.tl.score_genes(ad, [g for g in gs if g in ad.var_names], score_name=name)

    # fibroblast-only: module score by condition
    fib = ad[ad.obs.cell_type=="Fibroblast"]
    print(f"Fibroblasts: {fib.n_obs} cells  (by condition: "
          f"{dict(fib.obs.condition.value_counts())})")
    rows=[]
    for name in GS:
        for cond in ["Control","UV","UV_VitD"]:
            v=fib.obs.loc[fib.obs.condition==cond, name]
            rows.append({"module":name,"condition":cond,"mean":v.mean(),"n":len(v)})
        # UV vs Control test
        a=fib.obs.loc[fib.obs.condition=="Control",name]; b=fib.obs.loc[fib.obs.condition=="UV",name]
        u,p=stats.mannwhitneyu(a,b,alternative="two-sided")
        rows.append({"module":name,"condition":"UV_vs_Control_MWU_p","mean":p,"n":""})
    md=pd.DataFrame(rows); md.to_csv(RES/"GSE173385_fibroblast_module_scores.tsv",sep="\t",index=False)
    print("\nFibroblast module scores (mean by condition):")
    print(md[md.condition.isin(["Control","UV","UV_VitD"])].pivot(index="module",columns="condition",values="mean").to_string())

    # ---- Panel B figure ----
    fig=plt.figure(figsize=(15,8))
    ax1=fig.add_subplot(2,3,1); ax2=fig.add_subplot(2,3,2)
    sc.pl.umap(ad,color="cell_type",ax=ax1,show=False,legend_loc="on data",legend_fontsize=7,title="B1  cell types",frameon=False)
    sc.pl.umap(ad,color="condition",ax=ax2,show=False,title="B2  condition",frameon=False)
    ax3=fig.add_subplot(2,3,3)
    sc.pl.umap(ad,color="ECM_collagen",ax=ax3,show=False,title="B3  ECM/collagen score",frameon=False,cmap="RdBu_r")
    # violins: fibroblast module by condition
    order=["Control","UV","UV_VitD"]; colors={"Control":"#4daf4a","UV":"#e41a1c","UV_VitD":"#377eb8"}
    for i,name in enumerate(GS):
        ax=fig.add_subplot(2,3,4+i)
        data=[fib.obs.loc[fib.obs.condition==c,name].values for c in order]
        parts=ax.violinplot(data,showmeans=True,showextrema=False)
        for pc,c in zip(parts["bodies"],order): pc.set_facecolor(colors[c]); pc.set_alpha(.6)
        ax.set_xticks([1,2,3]); ax.set_xticklabels(order,fontsize=8)
        a=fib.obs.loc[fib.obs.condition=="Control",name]; b=fib.obs.loc[fib.obs.condition=="UV",name]
        _,p=stats.mannwhitneyu(a,b);
        ax.set_title(f"{name}\nfibroblasts; UV vs Ctrl p={p:.1e}",fontsize=8.5)
        ax.set_ylabel("score")
    fig.suptitle("Panel B — UVB photoaging mouse skin scRNA (GSE173385): fibroblast ECM/adhesion/senescence programs",fontsize=11,y=1.0)
    fig.tight_layout()
    fig.savefig(FIG/"PanelB_scRNA_UVB.pdf",bbox_inches="tight"); fig.savefig(FIG/"PanelB_scRNA_UVB.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    ad.write(PROC/"GSE173385_annotated.h5ad")
    print(f"\n[written] {FIG/'PanelB_scRNA_UVB.pdf'}")
    print(f"[written] {PROC/'GSE173385_annotated.h5ad'}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
