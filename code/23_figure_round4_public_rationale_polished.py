#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
23_figure_round4_public_rationale_polished.py  --  ROUND 4 polish of the public-data Figure 1
=============================================================================================
Polish-only revision of 17_figure_round3_public_rationale.py. Same public-data-only content and
the same exploratory caveats; only readability/typography are changed, per Round-4 feedback:

  1. SHORTER in-figure title (full nuance stays in the caption).
  2. Panel B right (violin) -- delta-mean and the exploratory caveat are moved OUT of the violin
     body: the delta sits in clear headroom above the violins, the caveat is a footnote below the
     axis. Violin density is never covered.
  3. Panel B condition labels -> "Normal (Control)" / "Photoaged (UVB)" (the public dataset's third
     combination arm stays absent).
  4. Panel E -> a clean 3-step schematic with larger fonts and three short component chips.

Still PUBLIC-DATA ONLY. No experimental data from this study appears in this figure. Terms referring
to experimental groups specific to this study, the public dataset's third (combination) condition arm,
or post-hoc verbs (rescue/restore/reverse) are kept out of the figure. ("remodeling" is the biological
ECM-remodeling term.)

  conda run -n skin-ecm-bioinformatics python 04_code/23_figure_round4_public_rationale_polished.py

Inputs (already on disk; nothing re-downloaded, nothing fabricated):
  06_data_processed/GSE173385_annotated.h5ad
  07_results/PanelC_enrichment_terms.tsv
  06_data_processed/public_ecm_aging_signature.tsv
Outputs:
  08_figures/revision_round_4_public_reanalysis_figure_polish_20260612/Figure1_public_ECM_rationale_polished.pdf / .png
"""
from __future__ import annotations
import os, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
warnings.filterwarnings("ignore"); sc.settings.verbosity=0

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
    "axes.linewidth":0.8,"axes.edgecolor":"#222","xtick.labelsize":7.5,"ytick.labelsize":7.5,
    "axes.titlesize":9.5,"axes.labelsize":8.5,"figure.dpi":120,
    "axes.spines.right":False,"axes.spines.top":False,
})
NORMAL_C, PHOTO_C = "#3b6fb0", "#e6731c"
PUB_DOWN, PUB_UP  = "#4575b4", "#d73027"
GO_BLUE, REAC_BLUE= "#4575b4", "#74add1"
BOX_PUB, BOX_AXIS, BOX_MAT, BOX_VAL = "#dbe7f1", "#c6dbef", "#f7e2d3", "#dceede"

ROOT=Path(os.environ.get("PROJECT_ROOT",Path(__file__).resolve().parent.parent)).resolve()
OUT=ROOT/"08_figures"/"revision_round_4_public_reanalysis_figure_polish_20260612"; OUT.mkdir(parents=True,exist_ok=True)

def _pre(s, prefixes): return any(s.startswith(p) for p in prefixes)
def matrisome_cat(sym:str)->str:
    s=str(sym).upper()
    if _pre(s,("LAMA","LAMB","LAMC")) or s in {"NID1","NID2","HSPG2","COL4A1","COL4A2","COL7A1","COL15A1","COL17A1","COL18A1"}:
        return "Basement membrane / laminins"
    if s.startswith("COL"): return "Collagens"
    if _pre(s,("ELN","FBN","FBLN","MFAP","LTBP","EMILIN")) or s in {"MGP","THSD4"}: return "Elastic fiber"
    if _pre(s,("DCN","LUM","FMOD","BGN","OGN","VCAN","ASPN","PRELP","PODN","ACAN","PRG4","HAPLN","GPC")) or s=="CILP":
        return "Proteoglycans"
    if _pre(s,("MMP","TIMP","LOX","CTS","HTRA","TGM","SERPIN","PCOLCE","PLOD")) or s in {"AMBP"} or s.startswith("PLA2"):
        return "ECM regulators"
    if _pre(s,("FN1","POSTN","TNC","TNX","THBS","VTN","DPT","MATN","CTHRC1","SPARC","NPNT","EFEMP")):
        return "ECM glycoproteins"
    return "Other ECM-assoc."
def net_dir(s):
    t=[x.split(":")[-1].strip().lower() for x in str(s).split(";") if ":" in x]
    return "down" if t.count("down")>t.count("up") else ("up" if t.count("up")>t.count("down") else "mixed")

# ==================== load data (public only) ====================
ad=sc.read_h5ad(ROOT/"06_data_processed"/"GSE173385_annotated.h5ad")
adB=ad[ad.obs.condition.isin(["Control","UV"])].copy()
fibB=adB[adB.obs.cell_type=="Fibroblast"].copy()

enr=pd.read_csv(ROOT/"07_results"/"PanelC_enrichment_terms.tsv",sep="\t")
enr["padj"]=pd.to_numeric(enr["Adjusted P-value"],errors="coerce")

sig=pd.read_csv(ROOT/"06_data_processed"/"public_ecm_aging_signature.tsv",sep="\t",dtype=str)
sig=sig[~sig["gene_symbol"].str.contains(r"\s|GLOBAL|ALL ",case=False,na=False)].copy()
sig["sym"]=sig["gene_symbol"].str.upper().str.replace(r"/.*","",regex=True)
sig["cat"]=sig["sym"].map(matrisome_cat); sig["aging"]=sig["directions"].map(net_dir)

# ==================== figure scaffold ====================
fig=plt.figure(figsize=(11.0,13.0))
gs=fig.add_gridspec(100,100,left=0.055,right=0.985,top=0.945,bottom=0.045)
def L(x,y,t): fig.text(x,y,t,fontsize=14,fontweight="bold",ha="left",va="top")

# ---------- Panel A: public-data rationale workflow ----------
axA=fig.add_subplot(gs[0:13,2:100]); axA.set_xlim(0,1); axA.set_ylim(0,1); axA.axis("off")
def box(ax,x,y,w,h,t,fc,fs=7.6,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.008,rounding_size=0.03",fc=fc,ec="#444",lw=1.0))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,fontweight="bold" if bold else "normal")
def arr(ax,x1,x2,y=0.5):
    ax.add_patch(FancyArrowPatch((x1,y),(x2,y),arrowstyle="-|>",mutation_scale=12,lw=1.3,color="#666"))
xs=[0.005,0.205,0.405,0.605,0.805]; w=0.185
labs=[("Public skin aging /\nphotoaging datasets","#dbe7f1"),
      ("Fibroblast ECM /\ncollagen program\nsuppression","#dbe7f1"),
      ("Structural\nmatrisome\ndecline","#dbe7f1"),
      ("ECM-collagen\nremodeling axis\nnominated",BOX_AXIS),
      ("SFMP-HA-SF\nECM-mimetic adhesive\ndesign rationale",BOX_MAT)]
for (lx,(t,fc)) in zip(xs,labs): box(axA,lx,0.27,w,0.52,t,fc)
for i in range(4): arr(axA,xs[i]+w+0.002,xs[i+1]-0.002,0.53)
axA.text(0.5,0.99,"Public databases nominate the target  (pre-experiment rationale)",
         ha="center",va="top",fontsize=8.2,style="italic",color="#555")
L(0.018,0.945,"A")

# ---------- Panel B (3 subpanels): normal vs photoaged only ----------
um=adB.obsm["X_umap"]
axB1=fig.add_subplot(gs[17:40,3:32]); L(0.018,0.812,"B")
order_ct=["Keratinocyte","Immune","Melanocyte","Fibroblast"]
cmap_ct={"Keratinocyte":"#cfd8dc","Immune":"#b0bec5","Melanocyte":"#90a4ae","Fibroblast":"#c0392b"}
for ct in order_ct:
    m=(adB.obs.cell_type==ct).values
    axB1.scatter(um[m,0],um[m,1],s=2 if ct!="Fibroblast" else 4,c=cmap_ct[ct],
                 label=ct,rasterized=True,linewidths=0,alpha=.8 if ct=="Fibroblast" else .5)
axB1.set_title("Skin cell classes (public scRNA, GSE173385)",fontsize=9); axB1.set_xticks([]); axB1.set_yticks([])
axB1.set_xlabel("UMAP1",fontsize=7.5); axB1.set_ylabel("UMAP2",fontsize=7.5)
for sp in ("left","bottom"): axB1.spines[sp].set_visible(True)
axB1.legend(fontsize=6.2,loc="lower left",frameon=False,markerscale=2,handletextpad=0.2)
fx,fy=um[(adB.obs.cell_type=="Fibroblast").values].mean(0)
axB1.annotate("Fibroblasts",(fx,fy),fontsize=7.5,fontweight="bold",color="#c0392b",ha="center")

axB2=fig.add_subplot(gs[17:40,37:64])
sca=axB2.scatter(um[:,0],um[:,1],s=2.5,c=adB.obs["ECM_collagen"],cmap="RdBu_r",
                 vmin=-1,vmax=2,rasterized=True,linewidths=0)
axB2.set_title("ECM/collagen module score",fontsize=9); axB2.set_xticks([]); axB2.set_yticks([])
axB2.set_xlabel("UMAP1",fontsize=7.5)
for sp in ("left","bottom"): axB2.spines[sp].set_visible(True)
cb=fig.colorbar(sca,ax=axB2,fraction=0.045,pad=0.02); cb.ax.tick_params(labelsize=6.5)

# B-iii fibroblast ECM module: normal vs photoaged -- annotations OUTSIDE the violin body
axB3=fig.add_subplot(gs[17:40,71:99])
order=["Control","UV"]; disp=["Normal\n(Control)","Photoaged\n(UVB)"]
cols_b={"Control":NORMAL_C,"UV":PHOTO_C}
data=[fibB.obs.loc[fibB.obs.condition==c,"ECM_collagen"].values for c in order]
parts=axB3.violinplot(data,showmeans=False,showextrema=False,widths=0.82)
for pc,c in zip(parts["bodies"],order): pc.set_facecolor(cols_b[c]); pc.set_alpha(.75); pc.set_edgecolor("#333"); pc.set_linewidth(.6)
for i,c in enumerate(order):
    v=fibB.obs.loc[fibB.obs.condition==c,"ECM_collagen"]; axB3.hlines(np.percentile(v,50),i+0.80,i+1.20,color="black",lw=1.6)
axB3.set_xticks([1,2]); axB3.set_xticklabels(disp,fontsize=7.8)
axB3.set_xlim(0.4,2.6)
axB3.set_ylabel("fibroblast ECM/collagen module score",fontsize=8)
axB3.set_title("Photoaging suppresses the fibroblast\nECM/collagen program",fontsize=9,pad=18)
mctrl=fibB.obs.loc[fibB.obs.condition=="Control","ECM_collagen"].mean()
muv  =fibB.obs.loc[fibB.obs.condition=="UV","ECM_collagen"].mean()
pcell=stats.mannwhitneyu(fibB.obs.loc[fibB.obs.condition=="Control","ECM_collagen"],
                         fibB.obs.loc[fibB.obs.condition=="UV","ECM_collagen"],alternative="two-sided")[1]
# add clear headroom above the violins, then put the delta-mean in that empty band (no density covered)
y0,y1=axB3.get_ylim(); axB3.set_ylim(y0, y1+0.42*(y1-y0))
axB3.annotate(f"$\\Delta$mean = {muv-mctrl:+.2f}   ({muv:.2f} vs {mctrl:.2f})",
              (0.5,0.985),xycoords="axes fraction",fontsize=7.4,fontweight="bold",color="#333",
              ha="center",va="top")
# exploratory caveat as a footnote BELOW the axis (off the data entirely)
axB3.annotate(f"exploratory: cell-level only, 1 sample / condition (MWU p={pcell:.1e})",
              (0.5,-0.185),xycoords="axes fraction",fontsize=6.0,style="italic",color="#666",
              ha="center",va="top",linespacing=0.95)

# ---------- Panel C ----------
axC=fig.add_subplot(gs[48:67,4:47]); L(0.018,0.530,"C")
CUR_C=[
    ("extracellular matrix organization (GO:0030198)","ECM organization (GO)","GO"),
    ("Extracellular Matrix Organization R-HSA-1474244","ECM organization (Reactome)","Reactome"),
    ("Collagen Formation R-HSA-1474290","Collagen formation","Reactome"),
    ("collagen fibril organization (GO:0030199)","Collagen fibril organization","GO"),
    ("Collagen Biosynthesis And Modifying Enzymes R-HSA-1650814","Collagen biosynthesis & modification","Reactome"),
    ("ECM Proteoglycans R-HSA-3000178","ECM proteoglycans","Reactome"),
    ("Elastic Fibre Formation R-HSA-1566948","Elastic fibre formation","Reactome"),
    ("Laminin Interactions R-HSA-3000157","Laminin interactions","Reactome"),
]
padj_map=dict(zip(enr["Term"],enr["padj"]))
cc=[(lab,src,-np.log10(max(float(padj_map[term]),1e-300)))
    for term,lab,src in CUR_C if term in padj_map and pd.notna(padj_map[term])]
yc=np.arange(len(cc))[::-1]
ccolors=[GO_BLUE if s=="GO" else REAC_BLUE for _,s,_ in cc]
axC.barh(yc,[v for *_,v in cc],color=ccolors,edgecolor="#222",lw=.4,height=.72)
axC.set_yticks(yc); axC.set_yticklabels([lab for lab,_,_ in cc],fontsize=7.2)
axC.axvline(-np.log10(0.05),color="#b2182b",ls="--",lw=.9)
axC.set_xlabel("$-\\log_{10}$ adjusted p",fontsize=8)
axC.set_title("Photoaging-suppressed fibroblast genes enrich\nfor ECM / collagen organization",fontsize=9)
axC.legend(handles=[Patch(color=GO_BLUE,label="GO-BP 2021"),Patch(color=REAC_BLUE,label="Reactome 2022")],
           fontsize=6.0,loc="lower right",frameon=False,title="source",title_fontsize=6.2)
axC.annotate("mouse genes via Enrichr; Reactome IDs are human reference R-HSA; red dashed = adj p 0.05",
             (0.98,-0.27),xycoords="axes fraction",fontsize=5.7,ha="right",style="italic",color="#555")

# ---------- Panel D ----------
axD=fig.add_subplot(gs[48:67,56:99]); L(0.535,0.530,"D")
cats=["Collagens","Basement membrane / laminins","ECM glycoproteins","Proteoglycans","Elastic fiber","ECM regulators"]
dn=[int(((sig.cat==c)&(sig.aging=="down")).sum()) for c in cats]
upc=[int(((sig.cat==c)&(sig.aging=="up")).sum()) for c in cats]
yy=np.arange(len(cats))[::-1]
axD.barh(yy,[-x for x in dn],color=PUB_DOWN,edgecolor="#222",lw=.4,height=.7,label="down with aging")
axD.barh(yy,upc,color=PUB_UP,edgecolor="#222",lw=.4,height=.7,label="up with aging")
axD.set_yticks(yy); axD.set_yticklabels(cats,fontsize=7.4); axD.axvline(0,color="#222",lw=.8)
for i,(d,u) in zip(yy,zip(dn,upc)):
    if d: axD.text(-d-0.3,i,str(d),va="center",ha="right",fontsize=6.5,color=PUB_DOWN)
    if u: axD.text(u+0.3,i,str(u),va="center",ha="left",fontsize=6.5,color=PUB_UP)
axD.set_xlim(-15.5,6.5)
axD.set_xlabel("# signature genes ($\\leftarrow$down with aging | up$\\rightarrow$)",fontsize=8)
axD.set_title("Human public skin aging/photoaging:\nbroad structural matrisome decline",fontsize=9)
axD.legend(fontsize=6.6,loc="upper right",frameon=True,facecolor="white",edgecolor="#cccccc",framealpha=0.92)
axD.annotate("129-gene signature curated from published human skin aging / photoaging studies",
             (0.5,-0.27),xycoords="axes fraction",fontsize=5.9,ha="center",style="italic",color="#555")

# ---------- Panel E: clean 3-step design-rationale schematic ----------
axE=fig.add_subplot(gs[74:97,2:100]); axE.set_xlim(0,1); axE.set_ylim(0,1); axE.axis("off"); L(0.018,0.290,"E")
axE.text(0.5,0.99,"Design rationale  (how the public-nominated axis motivates the material -- not a Figure 1 result)",
         ha="center",va="top",fontsize=8.4,style="italic",color="#555")
# three larger step boxes
box(axE,0.020,0.50,0.255,0.34,"Public-nominated\nECM-collagen\nremodeling axis",BOX_AXIS,fs=9.2,bold=True)
arr(axE,0.280,0.360,0.67)
box(axE,0.370,0.50,0.255,0.34,"SFMP-HA-SF\nECM-mimetic\nadhesive hydrogel",BOX_MAT,fs=9.2,bold=True)
arr(axE,0.630,0.730,0.67)
box(axE,0.735,0.50,0.250,0.34,"Validate in this study:\ncollagen I / III / IV,\nfibronectin, vinculin,\nhistology",BOX_VAL,fs=8.2)
# three simplified component chips under the material step (bigger text, short labels)
chips=[("HA","hydration / filling"),
       ("Silk fibroin","ECM-like matrix"),
       ("SF microspheres","adhesive microinterfaces")]
cw=0.150; cx=[0.300,0.470,0.640]
axE.text(0.4975,0.455,"three components:",ha="center",va="top",fontsize=7.0,style="italic",color="#7a4a22")
for (h,t),x in zip(chips,cx):
    axE.add_patch(FancyBboxPatch((x,0.07),cw,0.30,boxstyle="round,pad=0.006,rounding_size=0.025",
                                 fc="#fbf1e7",ec="#caa57f",lw=0.9))
    axE.text(x+cw/2,0.315,h,ha="center",va="top",fontsize=7.8,fontweight="bold",color="#9c5a22")
    axE.text(x+cw/2,0.20,t,ha="center",va="center",fontsize=7.0,color="#7a4a22",linespacing=0.95)

# ---------- shortened title + honest footer ----------
fig.suptitle("Public skin photoaging data nominate ECM-collagen remodeling for hydrogel design",
             fontsize=13.5,fontweight="bold",y=0.972)
fig.text(0.055,0.012,"Public-data rationale only (no experimental data from this study on this figure). "
         "scRNA module-score comparisons are exploratory (cell-level; one biological sample per condition). "
         "Panel E is a design rationale, not a result.",
         fontsize=6.6,style="italic",color="#555")

fig.savefig(OUT/"Figure1_public_ECM_rationale_polished.pdf",bbox_inches="tight")
fig.savefig(OUT/"Figure1_public_ECM_rationale_polished.png",dpi=300,bbox_inches="tight")
print("[written]",OUT/"Figure1_public_ECM_rationale_polished.pdf")
print("Panel B normal vs photoaged fibroblast ECM/collagen mean:",round(mctrl,3),round(muv,3),
      "| delta",round(muv-mctrl,3),"| cell-level MWU p",f"{pcell:.3e}")
print("Panel C terms:",[(lab,round(v,2)) for lab,_,v in cc])
print("Panel D categories (down/up):",list(zip(cats,dn,upc)))
