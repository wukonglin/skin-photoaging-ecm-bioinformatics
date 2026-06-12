#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13_panelC_enrichment.py  --  Panel C (enrichment of UV-altered fibroblast genes)
================================================================================
Take the genes DOWN-regulated in fibroblasts under UV (from the GSE173385 scRNA
analysis) and test which pathways they enrich for. Headline: ECM / collagen /
adhesion programs are over-represented among the UV-suppressed fibroblast genes.

  conda run -n skin-ecm-bioinformatics python 04_code/13_panelC_enrichment.py

Input:  06_data_processed/GSE173385_annotated.h5ad  (from step 12)
Outputs:
  07_results/GSE173385_fibroblast_UVdown_DE.tsv
  07_results/PanelC_enrichment_terms.tsv
  08_figures/PanelC_enrichment.pdf/.png   (overwrites the bulk placeholder)
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); sc.settings.verbosity=0

ROOT=Path(os.environ.get("PROJECT_ROOT",Path(__file__).resolve().parent.parent)).resolve()
H5=ROOT/"06_data_processed"/"GSE173385_annotated.h5ad"
RES=ROOT/"07_results"; FIG=ROOT/"08_figures"
RES.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)

# offline curated gene-set ORA fallback sets (mouse symbols)
CURATED={
 "ECM organization / collagen":["Col1a1","Col1a2","Col3a1","Col4a1","Col5a1","Col5a2","Col6a1","Col6a2","Col6a3","Col14a1","Fn1","Lox","Loxl1","Postn","Sparc","Bgn","Dcn","Lum","Mmp2","Timp1","Lama1","Lamb1","Lamc1","Nid1","Hspg2","Eln","Fbn1"],
 "Focal adhesion / integrin":["Itga1","Itga3","Itga5","Itga8","Itgav","Itgb1","Ilk","Fermt2","Lims1","Lims2","Parva","Vcl","Tln1","Pxn","Actn1","Flna","Cav1","Zyx","Ptk2"],
 "Cytoskeleton / actin":["Acta2","Actb","Actg1","Myl9","Myh9","Tpm1","Tpm2","Cnn1","Tagln","Vim","Flna","Actn1"],
 "Senescence / SASP / ROS":["Cdkn1a","Cdkn2a","Trp53","Il6","Cxcl1","Serpine1","Gdf15","Apoe","Jun","Fos","Sod1","Sod2","Gpx1","Cat"],
}

def main():
    if not H5.exists(): sys.exit(f"[MISSING] {H5} (run 12_scRNA_UVB_GSE173385.py first)")
    ad=sc.read_h5ad(H5)
    fib=ad[ad.obs.cell_type=="Fibroblast"].copy()
    fib=fib[fib.obs.condition.isin(["Control","UV"])].copy()
    sc.tl.rank_genes_groups(fib,"condition",groups=["UV"],reference="Control",method="wilcoxon")
    de=sc.get.rank_genes_groups_df(fib,group="UV")
    de.to_csv(RES/"GSE173385_fibroblast_UVdown_DE.tsv",sep="\t",index=False)
    down=de[(de.logfoldchanges<-0.25)&(de.pvals_adj<0.05)]["names"].tolist()
    up  =de[(de.logfoldchanges> 0.25)&(de.pvals_adj<0.05)]["names"].tolist()
    print(f"Fibroblast UV vs Control: {len(down)} down, {len(up)} up (|lfc|>0.25, padj<0.05)")
    bg=set(fib.var_names)

    terms=None
    # try Enrichr (online); fall back to curated ORA (offline, hypergeometric)
    try:
        import gseapy as gp
        enr=gp.enrichr(gene_list=down, organism="mouse",
                       gene_sets=["GO_Biological_Process_2021","Reactome_2022"],
                       outdir=str(RES/"enrichr_tmp"), no_plot=True)
        r=enr.results.copy()
        r["Adjusted P-value"]=pd.to_numeric(r["Adjusted P-value"],errors="coerce")
        kw=r["Term"].str.contains("extracellular matrix|collagen|integrin|adhesion|laminin|ECM|cytoskel|actin|elastic|basement",case=False,na=False)
        terms=r[kw].sort_values("Adjusted P-value").head(14)[["Term","Adjusted P-value","Overlap","Gene_set"]]
        terms.to_csv(RES/"PanelC_enrichment_terms.tsv",sep="\t",index=False)
        terms["neglog10"]=-np.log10(terms["Adjusted P-value"].clip(lower=1e-300))
        labels=terms["Term"].str.replace(r"\(GO:\d+\)","",regex=True).str.slice(0,52)
        vals=terms["neglog10"]; src="Enrichr (GO_BP + Reactome), mouse"
        print("[enrichr OK]")
    except Exception as e:
        print(f"[enrichr unavailable -> offline curated ORA] {e}")
        from scipy.stats import hypergeom
        N=len(bg); downset=set(down); rows=[]
        for name,genes in CURATED.items():
            gs=set(g for g in genes if g in bg); K=len(gs)
            k=len(downset & gs); n=len(downset)
            if K>=3:
                p=hypergeom.sf(k-1,N,K,n)
                rows.append({"Term":name,"Adjusted P-value":p,"Overlap":f"{k}/{K}","Gene_set":"curated"})
        terms=pd.DataFrame(rows).sort_values("Adjusted P-value")
        terms.to_csv(RES/"PanelC_enrichment_terms.tsv",sep="\t",index=False)
        terms["neglog10"]=-np.log10(terms["Adjusted P-value"].clip(lower=1e-300))
        labels=terms["Term"]; vals=terms["neglog10"]; src="curated ORA (hypergeometric)"

    # figure
    fig,ax=plt.subplots(figsize=(7.5,0.5*len(terms)+1.4))
    y=np.arange(len(terms))[::-1]
    ax.barh(y,vals,color="#2166ac",edgecolor="black",lw=.4)
    ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=8)
    ax.axvline(-np.log10(0.05),color="red",ls="--",lw=.8,label="FDR 0.05")
    ax.set_xlabel("-log10(adjusted p)"); ax.legend(fontsize=7,frameon=False)
    ax.set_title(f"Panel C — pathways enriched among UV-suppressed fibroblast genes\n(GSE173385; {src})",fontsize=9.5)
    fig.tight_layout(); fig.savefig(FIG/"PanelC_enrichment.pdf",bbox_inches="tight"); fig.savefig(FIG/"PanelC_enrichment.png",dpi=300,bbox_inches="tight"); plt.close(fig)
    print(terms[["Term","Adjusted P-value","Overlap"]].to_string(index=False))
    print(f"\n[written] {FIG/'PanelC_enrichment.pdf'}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
