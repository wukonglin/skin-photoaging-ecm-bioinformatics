#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_public_supplement_validation_figure.py  --  ROUND 3 public-data audit: supplementary
public-data validation figure (PUBLIC DATA ONLY; no experimental data from this study)
========================================================================================
Decision (see PublicData_FigureDecision_EN.md): the Round-3 public-only Figure 1 is KEPT as the
recommended front figure; new 2022-2026 public datasets are NOT crammed into Figure 1. This
SUPPLEMENTARY figure strengthens the public rationale by (i) triangulating the ECM/collagen-decline
argument across independent in-vivo sources, (ii) adding a NEW replicate-aware in-vivo reanalysis
(GSE284483, mouse UVB, n=3/group) that the single-sample GSE173385 could not provide, and
(iii) HONESTLY bounding the claim -- cultured fibroblasts and whole-skin bulk dilute the
fibroblast-specific signal, which is exactly why the rationale rests on fibroblast-resolved /
ECM-resolved in-vivo data.

Panels:
  A  Cross-source concordance map: independent public evidence lines vs ECM/collagen readouts
  B  GSE113957 (human dermal fibroblast aging, n=133): in-vitro caveat -- senescence up, ECM flat
  C  GSE284483 (mouse UVB, in vivo, n=3/group) NEW replicate-aware reanalysis: per-gene log2FC --
     basement-membrane collagens + proteoglycans down, elastin up (solar elastosis) = matrisome
     remodeling (module-level n.s. because whole-skin bulk dilutes the fibroblast-specific signal)
  D  Screened 2018-2026 public dataset landscape (from public_dataset_candidates_2022_2026.tsv)

PUBLIC-DATA ONLY. No experimental data from this study appears here.

  conda run -n junye_ecm python 04_code/22_public_supplement_validation_figure.py

Inputs (all on disk; GSE284483 was downloaded as public data, originals untouched):
  07_results/GSE113957_module_vs_age.tsv
  07_results/GSE173385_fibroblast_module_scores.tsv
  07_results/GSE284483_UVB_ECM_module.tsv
  07_results/GSE284483_UVB_collagen_genes.tsv
  09_outputs/revision_round_3_public_data_audit_20260612/public_dataset_candidates_2022_2026.tsv
Outputs:
  08_figures/revision_round_3_public_data_audit_20260612/FigureS_public_dataset_validation.pdf / .png
"""
from __future__ import annotations
import os, re, warnings
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
SIG_C, NS_C = "#2c5f8a", "#b9c2cb"

ROOT=Path(os.environ.get("PROJECT_ROOT",Path(__file__).resolve().parent.parent)).resolve()
OUT=ROOT/"08_figures"/"revision_round_3_public_data_audit_20260612"; OUT.mkdir(parents=True,exist_ok=True)
AUDIT=ROOT/"09_outputs"/"revision_round_3_public_data_audit_20260612"
RES=ROOT/"07_results"

# ==================== load local, already-computed public results ====================
mv=pd.read_csv(RES/"GSE113957_module_vs_age.tsv",sep="\t")
mv["spearman_r"]=pd.to_numeric(mv["spearman_r"],errors="coerce"); mv["p"]=pd.to_numeric(mv["p"],errors="coerce")
rmap={r["module"]:(r["spearman_r"],r["p"]) for _,r in mv.iterrows()}

ms=pd.read_csv(RES/"GSE173385_fibroblast_module_scores.tsv",sep="\t")
ms["mean"]=pd.to_numeric(ms["mean"],errors="coerce")
def _ms(mod,cond):
    s=ms[(ms["module"]==mod)&(ms["condition"]==cond)]["mean"]; return float(s.iloc[0]) if len(s) else np.nan
ecm_ctrl,ecm_uv=_ms("ECM_collagen","Control"),_ms("ECM_collagen","UV")

# GSE284483 replicate-aware reanalysis
mod284=pd.read_csv(RES/"GSE284483_UVB_ECM_module.tsv",sep="\t")
per=mod284[mod284["sample"]!="__GROUP_TEST__"].copy()
per["module_score"]=pd.to_numeric(per["module_score"],errors="coerce")
gc=per.loc[per["group"].str.startswith("NoUVB"),"module_score"].values
gp=per.loc[per["group"].str.startswith("UVB"),"module_score"].values
t284,p284=stats.ttest_ind(gc,gp,equal_var=False)
gen284=pd.read_csv(RES/"GSE284483_UVB_collagen_genes.tsv",sep="\t")
gen284["lfc"]=pd.to_numeric(gen284["log2FC_UVBvsNoUVB"],errors="coerce")

# ==================== figure scaffold (2x2) ====================
fig=plt.figure(figsize=(12.0,9.6))
gs=fig.add_gridspec(100,100,left=0.075,right=0.985,top=0.905,bottom=0.075,hspace=0.0,wspace=0.0)
def L(x,y,t): fig.text(x,y,t,fontsize=14,fontweight="bold",ha="left",va="top")

# ---------------- Panel A: cross-source concordance map ----------------
axA=fig.add_subplot(gs[0:46,0:45]); axA.set_xlim(0,3); axA.set_ylim(0,3.2); axA.axis("off"); L(0.020,0.900,"A")
axA.set_title("Independent public evidence converges on ECM/collagen decline\nin vivo — honestly bounded in vitro and in bulk",
              fontsize=9.0,loc="left",pad=4)
rows=[("Mouse UVB skin, in vivo\nscRNA, fibroblast-resolved\n(GSE173385)","m1"),
      ("Human skin aging/photoaging,\nin vivo proteomic signature\n(129 genes; PXD015982/018430)","m2"),
      ("Mouse UVB skin, in vivo\nBULK, replicated n=3\n(GSE284483, this audit)","m3"),
      ("Human dermal fibroblast,\nin vitro bulk aging\n(GSE113957, n=133)","m4")]
cols=["Fibroblast\nECM/collagen\nprogram","Structural\nmatrisome\n(coll/BM/PG/elastin)","ECM/collagen\npathway\nenrichment"]
cells={
 "m1":[("down",f"↓ Δ={ecm_uv-ecm_ctrl:+.2f}"),("na","—"),("down","↓ ECM org\nGO p=6e-14")],
 "m2":[("down","↓ net"),("down","↓ 10/1 coll;\n11/0 BM; 9/1 PG"),("na","—")],
 "m3":[("flat","~ bulk\n(diluted)"),("down","↓ BM Col4, PG\nDcn/Lum; Eln↑"),("na","—")],
 "m4":[("flat",f"~ r={rmap['ECM_collagen'][0]:+.2f}\n(ns)"),("flat","~ down-sig\nr=+0.25 (not↓)"),("na","—")],
}
cc={"down":DOWN_C,"up":UP_C,"flat":FLAT_C,"na":NA_C}
cw,ch=0.60,0.50; x0,y0=1.14,2.10
for j,cl in enumerate(cols):
    axA.text(x0+cw/2+j*cw, y0+ch+0.10, cl, ha="center", va="bottom", fontsize=6.3, fontweight="bold", linespacing=0.95)
for i,(rl,key) in enumerate(rows):
    yy=y0-i*ch
    axA.text(x0-0.05, yy+ch/2, rl, ha="right", va="center", fontsize=5.9, linespacing=0.95)
    for j,(dirn,txt) in enumerate(cells[key]):
        xx=x0+j*cw
        axA.add_patch(FancyBboxPatch((xx+0.03,yy+0.03),cw-0.06,ch-0.06,
            boxstyle="round,pad=0.004,rounding_size=0.025",fc=cc[dirn],ec="#5b6670",lw=0.8,
            alpha=0.92 if dirn!="na" else 0.5))
        tcol="white" if dirn=="down" else "#3a3f44"
        axA.text(xx+cw/2, yy+ch/2, txt, ha="center", va="center", fontsize=5.3,
                 color=tcol, fontweight="bold", linespacing=0.92)
axA.legend(handles=[Patch(fc=DOWN_C,label="suppressed / down with aging"),
                    Patch(fc=FLAT_C,label="flat / diluted / discordant"),
                    Patch(fc=NA_C,label="not assessed in this source")],
           fontsize=5.9,loc="lower center",bbox_to_anchor=(0.52,0.005),ncol=1,frameon=False)

# ---------------- Panel B: GSE113957 module-vs-age Spearman r ----------------
axB=fig.add_subplot(gs[6:42,58:99]); L(0.545,0.900,"B")
modB=[("Public aging-DOWN signature","Public_ECM_aging_DOWN"),
      ("Public aging-UP signature","Public_ECM_aging_UP"),
      ("Senescence / SASP / ROS","Senescence_ROS_SASP"),
      ("Focal adhesion (curated)","Adhesion_focal_adhesion"),
      ("ECM / collagen (curated)","ECM_collagen")]
labs=[m[0] for m in modB]; rs=[rmap[m[1]][0] for m in modB]; ps=[rmap[m[1]][1] for m in modB]
yb=np.arange(len(modB))[::-1]
axB.barh(yb,rs,color=[SIG_C if p<0.05 else NS_C for p in ps],edgecolor="#222",lw=.5,height=.66)
axB.axvline(0,color="#222",lw=.9); axB.set_yticks(yb); axB.set_yticklabels(labs,fontsize=7.0)
axB.set_xlim(-0.30,0.42); axB.set_xlabel("Spearman r (module score vs donor age)",fontsize=8)
axB.set_title("In-vitro caveat — cultured human dermal fibroblasts (GSE113957, n=133):\n"
              "senescence rises with age, ECM/collagen does not decline",fontsize=8.6)
for y,r,p in zip(yb,rs,ps):
    star="**" if p<0.01 else ("*" if p<0.05 else "ns")
    axB.text(r+(0.012 if r>=0 else -0.012), y, f"r={r:+.2f} {star}", va="center",
             ha="left" if r>=0 else "right", fontsize=6.2, color="#333")
axB.legend(handles=[Patch(fc=SIG_C,label="p < 0.05"),Patch(fc=NS_C,label="ns")],
           fontsize=6.2,loc="lower right",frameon=False)
axB.annotate("Literature 'down-with-aging' genes do NOT decline in cultured fibroblasts —\n"
             "in-vitro fibroblast bulk is not used to claim ECM decline.",
             (0.0,-0.255),xycoords="axes fraction",fontsize=5.8,style="italic",color="#555",linespacing=1.0)

# ---------------- Panel C: GSE284483 NEW replicate-aware reanalysis ----------------
axC=fig.add_subplot(gs[58:100,0:42]); L(0.020,0.470,"C")
CAT={"Col4a1":"BM collagen","Col4a2":"BM collagen",
     "Col1a1":"Fibrillar collagen","Col1a2":"Fibrillar collagen","Col3a1":"Fibrillar collagen",
     "Col5a1":"Fibrillar collagen","Col6a1":"Fibrillar collagen",
     "Dcn":"Proteoglycan","Lum":"Proteoglycan","Fn1":"Glycoprotein","Postn":"Glycoprotein",
     "Eln":"Elastic fiber","Fbn1":"Elastic fiber"}
CAT_C={"BM collagen":"#1f6f8b","Fibrillar collagen":"#5aa1c0","Proteoglycan":"#3b6fb0",
       "Glycoprotein":"#9aa3ab","Elastic fiber":"#c98a2b"}
g=gen284.dropna(subset=["lfc"]).copy(); g["cat"]=g["gene"].map(CAT)
g=g.sort_values("lfc")
yy=np.arange(len(g))
axC.barh(yy,g["lfc"].values,color=[CAT_C.get(c,"#888") for c in g["cat"]],edgecolor="#222",lw=.45,height=.7)
axC.axvline(0,color="#222",lw=.9); axC.set_yticks(yy); axC.set_yticklabels(g["gene"].values,fontsize=7.0,fontstyle="italic")
axC.set_xlabel("log2 fold-change (UVB-saline / no-UVB)",fontsize=8)
axC.set_title("NEW replicate-aware in-vivo reanalysis — mouse UVB skin (GSE284483, n=3/group):\n"
              "basement-membrane collagens + proteoglycans down, elastin up (solar elastosis)",fontsize=8.4)
seen=[]
for c in ["BM collagen","Fibrillar collagen","Proteoglycan","Glycoprotein","Elastic fiber"]:
    if (g["cat"]==c).any(): seen.append(c)
axC.legend(handles=[Patch(fc=CAT_C[c],label=c) for c in seen],fontsize=5.8,loc="lower right",frameon=False,
           title="matrisome category",title_fontsize=6.0)
axC.annotate(f"whole-skin module-level n.s. (Welch t p={p284:.2f}, n=3 vs 3): the fibroblast-specific\n"
             "suppression is diluted in bulk tissue — consistent with the Panel A bound.",
             (0.0,-0.235),xycoords="axes fraction",fontsize=5.7,style="italic",color="#555",linespacing=1.0)

# ---------------- Panel D: screened public dataset landscape ----------------
axD=fig.add_subplot(gs[58:100,55:99]); L(0.545,0.470,"D")
cand_path=AUDIT/"public_dataset_candidates_2022_2026.tsv"
dec_color={"use_in_Figure1_v2":"#1a9850","use_in_supplement":"#1f6f8b","cite_only":"#7a5195",
           "needs_manual_download_or_permission":"#e6731c","screened_out":"#b9c2cb"}
if cand_path.exists():
    cand=pd.read_csv(cand_path,sep="\t",dtype=str).fillna("")
    def _yr(s):
        m=re.search(r"(20\d{2})",str(s)); return int(m.group(1)) if m else np.nan
    cand["yr"]=cand["release_or_publication_year"].map(_yr)
    cand["rk"]=pd.to_numeric(cand["candidate_rank"],errors="coerce")
    top=cand.sort_values("rk").head(12).reset_index(drop=True)
    yy=np.arange(len(top))[::-1]
    for y,(_,row) in zip(yy,top.iterrows()):
        dec=row.get("decision","").strip(); col=dec_color.get(dec,"#555")
        xr=row["yr"] if not pd.isna(row["yr"]) else 2020
        axD.hlines(y,2013.5,xr,color=col,lw=1.0,alpha=0.45)
        axD.scatter([xr],[y],s=62,color=col,edgecolor="#222",lw=.6,zorder=3)
        sp=row.get("species","")[:5]
        axD.text(2013.2,y,row.get("accession_or_id",""),ha="right",va="center",fontsize=6.0,fontweight="bold")
        axD.text(xr+0.18,y,f"{sp} · {row.get('assay_type','')}",ha="left",va="center",fontsize=5.5,color="#444")
    axD.set_xlim(2009,2028.5); axD.set_ylim(-0.8,len(top)-0.2)
    axD.set_yticks([]); axD.set_xlabel("dataset release / publication year",fontsize=8)
    axD.set_title("Screened public datasets (top-ranked) and recommended use",fontsize=8.6,loc="left")
    seen=[d for d in dec_color if (top["decision"]==d).any()]
    axD.legend(handles=[Patch(fc=dec_color[d],label=d.replace("_"," ")) for d in seen],
               fontsize=5.6,loc="lower right",frameon=False,ncol=1)
    for s in ("left","top","right"): axD.spines[s].set_visible(False)
    axD.annotate("Full screening of 33 datasets in public_dataset_candidates_2022_2026.tsv",
                 (0.0,-0.235),xycoords="axes fraction",fontsize=5.7,style="italic",color="#555")
else:
    axD.axis("off"); axD.text(0.5,0.5,"candidates TSV not found",ha="center",va="center",color="#999")

fig.suptitle("Supplementary public-data validation: cross-source convergence on ECM/collagen "
             "decline, replicate-aware and honestly bounded",fontsize=11.6,fontweight="bold",y=0.962)
fig.text(0.075,0.016,"Public-data only (no experimental data from this study). GSE113957: replicate-aware "
         "(Spearman vs donor age, n=133). GSE284483: replicate-aware (n=3 vs 3, Welch t). GSE173385 module "
         "difference is exploratory (cell-level; one biological sample per condition).",
         fontsize=6.2,style="italic",color="#555")

fig.savefig(OUT/"FigureS_public_dataset_validation.pdf",bbox_inches="tight")
fig.savefig(OUT/"FigureS_public_dataset_validation.png",dpi=300,bbox_inches="tight")
print("[written]",OUT/"FigureS_public_dataset_validation.pdf")
print("GSE284483 module: NoUVB mean",round(gc.mean(),3),"UVB-saline mean",round(gp.mean(),3),
      "delta",round(gp.mean()-gc.mean(),3),"Welch p",round(p284,4))
print("GSE284483 per-gene log2FC:",list(zip(gen284['gene'],gen284['lfc'])))
print("GSE113957 r:",[(m[0],round(rmap[m[1]][0],3)) for m in modB])
