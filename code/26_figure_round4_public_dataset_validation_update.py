#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_figure_round4_public_dataset_validation_update.py  --  ROUND 4 public-data supplement update
===============================================================================================
Updated supplementary public-data validation figure. It ADDS the Round-4 reanalyses to the
public evidence base and, crucially, fills the gap the Round-3 supplement could not:

  * GSE110978 (mouse dermal-fibroblast aging microarray, n=4 young vs n=4 old) is the FIRST
    in-vivo, FIBROBLAST-RESOLVED, REPLICATE-AWARE ECM read-out in this package. The single-sample
    mouse UVB scRNA (GSE173385) and the diluted whole-skin bulk (GSE284483) could not provide
    replicate-level fibroblast statistics; cultured human fibroblasts (GSE113957) lose the in-vivo
    ECM decline. GSE110978 closes that gap with a clean, replicate-aware ECM/collagen decline.

  * GSE275491 (human in-vivo photoaging scRNA, single individual) is shown HONESTLY as an
    exploratory cross-check: per-fibroblast ECM module is essentially comparable between
    photoaged arm and normal buttock in this one individual (and is reported as exploratory),
    while photoaged skin captured ~40% fewer fibroblasts. It is NOT used to claim ECM decline.

Panels:
  A  Cross-source concordance map (updated): each public evidence line x ECM readouts, with an
     explicit "replicate-aware design?" column.
  B  GSE110978 per-gene ECM/collagen log2FC (old vs young), replicate-aware (Welch per gene).  [NEW headline]
  C  GSE110978 ECM module score: young (n=4) vs old (n=4), per-sample points + Welch test.      [NEW]
  D  GSE275491 exploratory single-individual cross-check (honestly bounded).                    [NEW]

PUBLIC-DATA ONLY. No experimental data from this study appears here. All numbers are read from the
Round-4 reanalysis tables (07_results/GSE110978_*.tsv, GSE275491_*.tsv); nothing is fabricated.

  conda run -n skin-ecm-bioinformatics python 04_code/26_figure_round4_public_dataset_validation_update.py

Inputs:
  07_results/GSE110978_mouse_fibroblast_aging_gene_trends.tsv
  07_results/GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv
  07_results/GSE275491_public_photoaging_reanalysis.tsv
Outputs:
  08_figures/revision_round_4_public_reanalysis_figure_polish_20260612/FigureS_public_dataset_validation_round4.pdf / .png
"""
from __future__ import annotations
import os, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
    "axes.linewidth":0.8,"axes.edgecolor":"#222","xtick.labelsize":7.5,"ytick.labelsize":7.5,
    "axes.titlesize":9.5,"axes.labelsize":8.5,"figure.dpi":120,
    "axes.spines.right":False,"axes.spines.top":False,
})
DOWN_C, FLAT_C, UP_C, NA_C = "#3b6fb0", "#9aa3ab", "#d73027", "#eef1f4"
YOUNG_C, OLD_C = "#3b6fb0", "#b2182b"

ROOT=Path(os.environ.get("PROJECT_ROOT",Path(__file__).resolve().parent.parent)).resolve()
OUT=ROOT/"08_figures"/"revision_round_4_public_reanalysis_figure_polish_20260612"; OUT.mkdir(parents=True,exist_ok=True)
RES=ROOT/"07_results"

# ==================== load Round-4 reanalysis tables ====================
gt=pd.read_csv(RES/"GSE110978_mouse_fibroblast_aging_gene_trends.tsv",sep="\t")
gt=gt[gt["direction"]!="not_on_array"].copy()
gt["log2FC_old_vs_young"]=pd.to_numeric(gt["log2FC_old_vs_young"],errors="coerce")
gt["p_value"]=pd.to_numeric(gt["p_value"],errors="coerce")

mod=pd.read_csv(RES/"GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv",sep="\t")
per=mod[~mod["sample"].astype(str).str.startswith("__")].copy()
per["module_score"]=pd.to_numeric(per["module_score"],errors="coerce")
y_scores=per.loc[per["group"]=="young","module_score"].values
o_scores=per.loc[per["group"]=="old","module_score"].values
t110,p110=stats.ttest_ind(o_scores,y_scores,equal_var=False)

g275=pd.read_csv(RES/"GSE275491_public_photoaging_reanalysis.tsv",sep="\t")
g275m={r["metric"]:r["value"] for _,r in g275.iterrows()}
def _f(k):
    try: return float(g275m[k])
    except: return np.nan
norm275=_f("fibroblast_ECM_mean_Normal_buttock"); photo275=_f("fibroblast_ECM_mean_Photoaged_arm")
nnorm=int(_f("fibroblast_n_Normal_buttock")); nphoto=int(_f("fibroblast_n_Photoaged_arm"))

# ==================== scaffold (2x2) ====================
fig=plt.figure(figsize=(12.0,9.8))
gs=fig.add_gridspec(100,100,left=0.075,right=0.985,top=0.905,bottom=0.075,hspace=0.0,wspace=0.0)
def L(x,y,t): fig.text(x,y,t,fontsize=14,fontweight="bold",ha="left",va="top")

# ---------------- Panel A: cross-source concordance (updated) ----------------
axA=fig.add_subplot(gs[0:48,0:46]); axA.set_xlim(0,3.2); axA.set_ylim(0,3.4); axA.axis("off"); L(0.020,0.900,"A")
axA.set_title("Independent public evidence converges on ECM/collagen decline in vivo;\n"
              "GSE110978 adds the missing fibroblast-resolved, replicate-aware read-out",
              fontsize=8.6,loc="left",pad=4)
rows=[("Mouse aged dermal fibroblasts,\nin vivo, REPLICATED n=4/4\n(GSE110978, this round)","r0"),
      ("Mouse UVB skin, in vivo scRNA,\nfibroblast-resolved, 1 sample\n(GSE173385)","r1"),
      ("Human photoaging skin, in vivo\nscRNA, 1 individual\n(GSE275491, this round)","r2"),
      ("Mouse UVB skin, in vivo BULK,\nreplicated n=3 (GSE284483)","r3"),
      ("Human aging/photoaging in vivo\nproteomic signature (129 genes)","r4"),
      ("Human dermal fibroblast,\nin vitro bulk aging (GSE113957)","r5")]
cols=["Fibroblast\nECM/collagen","Structural\nmatrisome","Replicate-\naware design"]
cells={
 "r0":[("down",f"down\nΔ={float(mod.loc[mod['sample']=='__delta_old_minus_young__','module_score'].iloc[0]):+.2f}"),
       ("down","17/19\ngenes down"),("down","YES\nn=4 vs 4")],
 "r1":[("down","down\n(exploratory)"),("na","—"),("flat","no\n1 sample")],
 "r2":[("flat","~ flat\n(exploratory)"),("na","—"),("flat","no\n1 individual")],
 "r3":[("flat","~ bulk\n(diluted)"),("down","BM coll, PG\ndown; Eln↑"),("down","YES\nn=3 vs 3")],
 "r4":[("down","↓ net"),("down","coll/BM/PG\ndown"),("down","cohort\nsignature")],
 "r5":[("flat","~ flat\n(in vitro)"),("flat","no in-vivo\ndecline"),("down","YES\nn=133")],
}
cc={"down":DOWN_C,"up":UP_C,"flat":FLAT_C,"na":NA_C}
cw,ch=0.62,0.46; x0,y0=1.30,2.55
for j,cl in enumerate(cols):
    axA.text(x0+cw/2+j*cw, y0+ch+0.06, cl, ha="center", va="bottom", fontsize=6.2, fontweight="bold", linespacing=0.95)
for i,(rl,key) in enumerate(rows):
    yy=y0-i*ch
    axA.text(x0-0.05, yy+ch/2, rl, ha="right", va="center", fontsize=5.7, linespacing=0.95,
             fontweight="bold" if key=="r0" else "normal")
    for j,(dirn,txt) in enumerate(cells[key]):
        xx=x0+j*cw
        axA.add_patch(FancyBboxPatch((xx+0.03,yy+0.03),cw-0.06,ch-0.06,
            boxstyle="round,pad=0.004,rounding_size=0.02",fc=cc[dirn],ec="#5b6670",lw=0.8,
            alpha=0.92 if dirn!="na" else 0.5))
        tcol="white" if dirn=="down" else "#3a3f44"
        axA.text(xx+cw/2, yy+ch/2, txt, ha="center", va="center", fontsize=5.2,
                 color=tcol, fontweight="bold", linespacing=0.92)
axA.legend(handles=[Patch(fc=DOWN_C,label="suppressed / down with aging"),
                    Patch(fc=FLAT_C,label="flat / diluted / exploratory"),
                    Patch(fc=NA_C,label="not assessed in this source")],
           fontsize=5.9,loc="lower center",bbox_to_anchor=(0.50,0.0),ncol=1,frameon=False)

# ---------------- Panel B: GSE110978 per-gene log2FC (NEW headline) ----------------
axB=fig.add_subplot(gs[6:44,58:99]); L(0.545,0.900,"B")
gb=gt.sort_values("log2FC_old_vs_young").copy()
ybar=np.arange(len(gb))
barc=[OLD_C if v>0 else DOWN_C for v in gb["log2FC_old_vs_young"]]
axB.barh(ybar,gb["log2FC_old_vs_young"].values,color=barc,edgecolor="#222",lw=.45,height=.72)
axB.axvline(0,color="#222",lw=.9)
axB.set_yticks(ybar); axB.set_yticklabels(gb["gene"].values,fontsize=7.0,fontstyle="italic")
axB.set_xlabel("log2 fold-change (old / young), per-gene",fontsize=8)
axB.set_title("GSE110978 — replicate-aware aged dermal fibroblasts (n=4 vs 4):\n"
              "17/19 ECM/collagen genes decline with aging",fontsize=8.4)
for yv,v,p in zip(ybar,gb["log2FC_old_vs_young"],gb["p_value"]):
    star="**" if p<0.01 else ("*" if p<0.05 else "")
    if star:
        axB.text(v+(0.03 if v<0 else 0.03)*np.sign(v) - (0.05 if v<0 else -0.02), yv, star,
                 va="center", ha="right" if v<0 else "left", fontsize=7.5, color="#333")
axB.annotate("* p<0.05  ** p<0.01 (Welch, per gene)  ·  blue = down with aging",
             (0.0,-0.13),xycoords="axes fraction",fontsize=6.0,style="italic",color="#555")

# ---------------- Panel C: GSE110978 module score young vs old (NEW) ----------------
axC=fig.add_subplot(gs[58:100,0:42]); L(0.020,0.470,"C")
rng=np.random.RandomState(0)
xs_y=1+ (rng.rand(len(y_scores))-0.5)*0.18
xs_o=2+ (rng.rand(len(o_scores))-0.5)*0.18
axC.scatter(xs_y,y_scores,s=58,color=YOUNG_C,edgecolor="#222",lw=.6,zorder=3,label="young (n=4)")
axC.scatter(xs_o,o_scores,s=58,color=OLD_C,edgecolor="#222",lw=.6,zorder=3,label="old (n=4)")
for xc,vals,col in [(1,y_scores,YOUNG_C),(2,o_scores,OLD_C)]:
    axC.hlines(vals.mean(),xc-0.22,xc+0.22,color="#222",lw=2.0,zorder=4)
axC.set_xticks([1,2]); axC.set_xticklabels(["young","old"],fontsize=8.5)
axC.set_xlim(0.5,2.5); axC.set_ylabel("ECM/collagen module score (per sample)",fontsize=8)
axC.set_title("GSE110978 — replicate-aware ECM module:\nclean young vs old separation",fontsize=8.6)
ymax=max(y_scores.max(),o_scores.max()); ymin=min(y_scores.min(),o_scores.min()); pad=0.18*(ymax-ymin)
axC.set_ylim(ymin-pad,ymax+pad*2.4)
axC.annotate(f"Welch t = {t110:.1f}, p = {p110:.1e}\n(old − young Δ = {o_scores.mean()-y_scores.mean():+.2f})",
             (0.5,0.985),xycoords="axes fraction",fontsize=7.2,fontweight="bold",ha="center",va="top",color="#333")
axC.legend(fontsize=6.6,loc="lower left",frameon=False)

# ---------------- Panel D: GSE275491 exploratory single-individual cross-check (NEW) ----------------
axD=fig.add_subplot(gs[58:100,55:99]); L(0.545,0.470,"D")
vals=[norm275,photo275]; labs=["Normal\n(buttock)","Photoaged\n(arm)"]; bc=[FLAT_C,FLAT_C]
xb=[1,2]
axD.bar(xb,vals,width=0.55,color=bc,edgecolor="#222",lw=.7)
for x,v,n in zip(xb,vals,[nnorm,nphoto]):
    axD.text(x,v+0.02,f"{v:.2f}\n(n={n} fib)",ha="center",va="bottom",fontsize=6.8,color="#333")
axD.set_xticks(xb); axD.set_xticklabels(labs,fontsize=8.0)
axD.set_xlim(0.4,2.6); axD.set_ylim(0,max(vals)*1.32)
axD.set_ylabel("fibroblast ECM/collagen module score",fontsize=8)
axD.set_title("GSE275491 — human in-vivo photoaging scRNA (1 individual):\n"
              "per-fibroblast ECM module ~ comparable (EXPLORATORY)",fontsize=8.2)
axD.annotate("Single individual, paired arm vs buttock → no biological replication; body site confounds\n"
             "sun exposure. Per-fibroblast ECM module is not lower in photoaged skin (exploratory), but\n"
             f"photoaged arm captured ~{100*(1-nphoto/nnorm):.0f}% fewer fibroblasts ({nphoto} vs {nnorm}). NOT used to claim ECM decline.",
             (0.0,-0.15),xycoords="axes fraction",fontsize=5.8,style="italic",color="#555",linespacing=1.05)

fig.suptitle("Supplementary public-data validation (Round 4): a replicate-aware, fibroblast-resolved "
             "ECM/collagen decline, honestly bounded",fontsize=11.4,fontweight="bold",y=0.962)
fig.text(0.075,0.014,"Public-data only (no experimental data from this study). GSE110978: replicate-aware in-vivo "
         "dermal-fibroblast aging (n=4 vs 4, Welch). GSE275491: single individual (exploratory, cell-level). "
         "Module gene set identical to Figure 1 Panel B.",fontsize=6.2,style="italic",color="#555")

fig.savefig(OUT/"FigureS_public_dataset_validation_round4.pdf",bbox_inches="tight")
fig.savefig(OUT/"FigureS_public_dataset_validation_round4.png",dpi=300,bbox_inches="tight")
print("[written]",OUT/"FigureS_public_dataset_validation_round4.pdf")
print("GSE110978 module: young",round(y_scores.mean(),3),"old",round(o_scores.mean(),3),
      "delta",round(o_scores.mean()-y_scores.mean(),3),"Welch p",f"{p110:.2e}")
print("GSE110978 per-gene down/up:",int((gt.direction=='down_with_aging').sum()),int((gt.direction=='up_with_aging').sum()))
print("GSE275491 fib module normal/photoaged:",round(norm275,3),round(photo275,3),"n",nnorm,nphoto)
