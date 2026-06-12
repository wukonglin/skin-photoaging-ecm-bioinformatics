#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
25_public_dataset_reanalysis_round4.py  --  ROUND 4 public-data reanalysis
==========================================================================
Two NEW public-data reanalyses that strengthen the public ECM/collagen-decline rationale with
data that the existing sources could not provide:

  (1) GSE110978  -- mouse dermal-fibroblast aging microarray (Affymetrix Mouse 430 2.0 PM, GPL11180).
      FIBROBLAST-RESOLVED and REPLICATE-AWARE (n=4 young vs n=4 old). This is the dataset class the
      Round-3 supplement lacked: an in-vivo, fibroblast-resolved, replicated ECM read-out (the
      single-sample mouse UVB scRNA and the diluted whole-skin bulk could not give replicate-level
      fibroblast statistics). Salzer et al., Cell 2018 (PMID 30415840).

  (2) GSE275491  -- human in-vivo photoaging scRNA-seq (10x), paired sun-exposed forearm (photoaged)
      vs unexposed buttock (normal) within ONE individual. Human in-vivo photoaging at single-cell
      resolution, directly on-axis (fibroblast + ECM). Because it is a single individual it carries the
      SAME exploratory caveat class as GSE173385 (cell-level only; no biological replication). Wu et al.
      (PMID 39540047).

The ECM/collagen gene module is identical to the one used for Figure 1 Panel B (script 12):
  Col1a1 Col1a2 Col3a1 Col4a1 Col4a2 Col5a1 Col6a1 Fn1 Dcn Lum Eln Fbn1 Postn Tnc Lox Mmp2 Mmp3 Mmp9 Timp1
(human reanalysis uses the upper-case human orthologs).

PUBLIC DATA ONLY. No experimental data generated in this study enters this script.
No values are fabricated: every number is computed from the downloaded public matrices.

  conda run -n skin-ecm-bioinformatics python 04_code/25_public_dataset_reanalysis_round4.py

Inputs (downloaded public data; originals untouched):
  05_data_raw/round4_public/GSE110978/GSE110978_series_matrix.txt.gz
  05_data_raw/round4_public/GSE110978/GPL11180.annot.gz
  05_data_raw/round4_public/GSE275491/GSE275491_{barcodes,features,matrix}.{tsv,mtx}.gz
Outputs:
  07_results/GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv
  07_results/GSE110978_mouse_fibroblast_aging_gene_trends.tsv
  07_results/GSE275491_public_photoaging_reanalysis.tsv
"""
from __future__ import annotations
import os, gzip, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse, stats
import scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity = 0

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parent.parent)).resolve()
RAW  = ROOT/"05_data_raw"/"round4_public"
RES  = ROOT/"07_results"; RES.mkdir(parents=True, exist_ok=True)

# Same ECM/collagen module as Figure 1 Panel B (mouse symbols); human orthologs are the upper-case forms.
ECM_MOUSE = ["Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col6a1","Fn1","Dcn","Lum",
             "Eln","Fbn1","Postn","Tnc","Lox","Mmp2","Mmp3","Mmp9","Timp1"]
ECM_HUMAN = [g.upper() for g in ECM_MOUSE]


# ======================================================================================
# (1) GSE110978  --  mouse dermal fibroblast aging microarray (replicate-aware)
# ======================================================================================
def reanalyze_GSE110978():
    print("\n" + "="*78 + "\n[GSE110978] mouse dermal fibroblast aging microarray (n=4 young vs n=4 old)\n" + "="*78)
    smf = RAW/"GSE110978"/"GSE110978_series_matrix.txt.gz"
    annf= RAW/"GSE110978"/"GPL11180.annot.gz"

    # --- sample group labels from series-matrix metadata (young=GSM..651-654, old=..655-658) ---
    gsm, grp = [], []
    with gzip.open(smf, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsm = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "group:" in line:
                grp = [x.strip().strip('"').replace("group:","").strip() for x in line.rstrip("\n").split("\t")[1:]]
    grp = ["young" if "young" in g else "old" for g in grp]
    sample_group = dict(zip(gsm, grp))
    print("  samples:", sample_group)

    # --- expression table (RMA log2, probes x samples) ---
    expr = pd.read_csv(smf, sep="\t", comment="!", index_col=0)
    expr.columns = [c.strip().strip('"') for c in expr.columns]
    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    expr.index = [str(i).strip().strip('"') for i in expr.index]
    print(f"  expression table: {expr.shape[0]} probes x {expr.shape[1]} samples; "
          f"value range {np.nanmin(expr.values):.2f}..{np.nanmax(expr.values):.2f} (RMA log2)")

    # --- probe -> gene symbol (GPL11180 GEO annotation) ---
    p2s = {}
    with gzip.open(annf, "rt") as fh:
        intab = False; hdr = None
        for line in fh:
            if line.startswith("!platform_table_begin"): intab = True; continue
            if line.startswith("!platform_table_end"): break
            if intab:
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = f; iid = hdr.index("ID"); isym = hdr.index("Gene symbol"); continue
                if len(f) > isym and f[isym]:
                    p2s[f[iid]] = f[isym].split("///")[0].strip()
    print(f"  probe->symbol annotation: {len(p2s)} probes mapped")

    young = [g for g in gsm if sample_group[g]=="young"]
    old   = [g for g in gsm if sample_group[g]=="old"]

    # --- collapse to one probe per ECM gene (max mean expression = most reliable probe) ---
    gene_rows = []
    gene_sample = {}  # gene -> per-sample log2 vector (Series indexed by gsm)
    for gene in ECM_MOUSE:
        probes = [p for p, s in p2s.items() if s == gene and p in expr.index]
        if not probes:
            gene_rows.append({"gene": gene, "n_probes": 0, "probe_used": "NA",
                              "young_mean_log2": np.nan, "old_mean_log2": np.nan,
                              "log2FC_old_vs_young": np.nan, "welch_t": np.nan, "p_value": np.nan,
                              "direction": "not_on_array"})
            continue
        sub = expr.loc[probes]
        best = sub.mean(axis=1).idxmax()           # highest-expressed probe for this gene
        vec = expr.loc[best]
        gene_sample[gene] = vec
        yv = vec[young].astype(float).values; ov = vec[old].astype(float).values
        t, p = stats.ttest_ind(ov, yv, equal_var=False)
        lfc = ov.mean() - yv.mean()                # values are log2 -> difference is log2FC
        gene_rows.append({"gene": gene, "n_probes": len(probes), "probe_used": best,
                          "young_mean_log2": round(yv.mean(),4), "old_mean_log2": round(ov.mean(),4),
                          "log2FC_old_vs_young": round(lfc,4), "welch_t": round(t,3),
                          "p_value": round(p,5),
                          "direction": "down_with_aging" if lfc < 0 else "up_with_aging"})
    gtrend = pd.DataFrame(gene_rows)
    gtrend.to_csv(RES/"GSE110978_mouse_fibroblast_aging_gene_trends.tsv", sep="\t", index=False)

    # --- replicate-aware module score: per-gene z across the 8 samples, averaged per sample ---
    gmat = pd.DataFrame(gene_sample).T            # gene x sample (only detected genes)
    gmat = gmat.reindex(columns=gsm)
    z = gmat.sub(gmat.mean(axis=1), axis=0).div(gmat.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    module = z.mean(axis=0)                        # per-sample module score
    yv = module[young].values; ov = module[old].values
    tt, pp = stats.ttest_ind(ov, yv, equal_var=False)

    mrows = [{"sample": g, "group": sample_group[g], "module_score": round(float(module[g]),4)} for g in gsm]
    mrows.append({"sample": "__young_mean__", "group": "young", "module_score": round(float(yv.mean()),4)})
    mrows.append({"sample": "__old_mean__",   "group": "old",   "module_score": round(float(ov.mean()),4)})
    mrows.append({"sample": "__delta_old_minus_young__", "group": "", "module_score": round(float(ov.mean()-yv.mean()),4)})
    mrows.append({"sample": "__welch_t__", "group": "", "module_score": round(float(tt),3)})
    mrows.append({"sample": "__welch_p_old_vs_young__", "group": "", "module_score": round(float(pp),5)})
    mrows.append({"sample": "__n_genes_in_module__", "group": "", "module_score": int(gmat.shape[0])})
    pd.DataFrame(mrows).to_csv(RES/"GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv", sep="\t", index=False)

    ndown = int((gtrend["direction"]=="down_with_aging").sum())
    nup   = int((gtrend["direction"]=="up_with_aging").sum())
    nsig  = int((gtrend["p_value"]<0.05).sum())
    print(f"  ECM module (replicate-aware): young {yv.mean():+.3f} vs old {ov.mean():+.3f}  "
          f"delta={ov.mean()-yv.mean():+.3f}  Welch t={tt:.2f}  p={pp:.4g}  (n=4 vs 4)")
    print(f"  per-gene: {ndown} down / {nup} up with aging; {nsig} genes p<0.05; "
          f"{int((gtrend['direction']=='not_on_array').sum())} not on array")
    print("  notable down genes:",
          list(gtrend.loc[(gtrend.direction=='down_with_aging')&(gtrend.p_value<0.10)].sort_values('p_value')['gene']))
    print("  [written]", RES/"GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv")
    print("  [written]", RES/"GSE110978_mouse_fibroblast_aging_gene_trends.tsv")
    return dict(delta=float(ov.mean()-yv.mean()), p=float(pp), ndown=ndown, nup=nup, nsig=nsig,
                young=float(yv.mean()), old=float(ov.mean()))


# ======================================================================================
# (2) GSE275491  --  human in-vivo photoaging scRNA (single individual; exploratory)
# ======================================================================================
def reanalyze_GSE275491():
    print("\n" + "="*78 + "\n[GSE275491] human photoaging scRNA: photoaged arm vs normal buttock (1 individual)\n" + "="*78)
    d = RAW/"GSE275491"
    # condition mapping comes from GEO sample order (GSM8478429 arm/photoaged first -> aggregation
    # suffix '-1'; GSM8478430 buttock/normal second -> suffix '-2'). Stated as an assumption below.
    SUFFIX_COND = {"1": "Photoaged_arm", "2": "Normal_buttock"}

    feats = pd.read_csv(d/"GSE275491_features.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(d/"GSE275491_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).values
    M = sc.read_mtx(d/"GSE275491_matrix.mtx.gz").T          # cells x genes
    ad = sc.AnnData(X=M.X.tocsr(),
                    obs=pd.DataFrame(index=barcodes),
                    var=pd.DataFrame(index=feats[1].astype(str).values))   # use gene symbols
    ad.var_names_make_unique()
    ad.obs["condition"] = [SUFFIX_COND.get(b.split("-")[-1], "unknown") for b in barcodes]
    print("  raw cells per condition:", dict(ad.obs.condition.value_counts()))

    # QC (same thresholds as the GSE173385 pipeline)
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(ad, min_genes=200); sc.pp.filter_genes(ad, min_cells=3)
    ad = ad[(ad.obs.n_genes_by_counts < 7000) & (ad.obs.pct_counts_mt < 25)].copy()
    print(f"  after QC: {ad.n_obs} cells x {ad.n_vars} genes; per condition {dict(ad.obs.condition.value_counts())}")

    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad); ad.raw = ad
    sc.pp.highly_variable_genes(ad, n_top_genes=2000)
    adh = ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(adh, max_value=10); sc.tl.pca(adh, n_comps=30)
    sc.pp.neighbors(adh, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adh, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adh)
    ad.obs["leiden"] = adh.obs["leiden"].values; ad.obsm["X_umap"] = adh.obsm["X_umap"]

    MARK = {"Fibroblast":["COL1A1","COL1A2","PDGFRA","DCN","LUM","COL3A1"],
            "Keratinocyte":["KRT14","KRT5","KRT10","KRT1"],
            "Immune":["PTPRC","CD52","LYZ","CD3D"],
            "Endothelial":["PECAM1","CDH5","VWF"],
            "Melanocyte":["MLANA","PMEL","TYR"],
            "Muscle":["ACTA2","MYL9","TAGLN"]}
    for ct, gs in MARK.items():
        sc.tl.score_genes(ad, [g for g in gs if g in ad.var_names], score_name=f"ct_{ct}")
    cl = ad.obs.groupby("leiden")[[f"ct_{c}" for c in MARK]].mean()
    lab = {c: cl.loc[c].idxmax().replace("ct_","") for c in cl.index}
    ad.obs["cell_type"] = ad.obs["leiden"].map(lab).astype("category")
    print("  cluster->cell type:", lab)

    sc.tl.score_genes(ad, [g for g in ECM_HUMAN if g in ad.var_names], score_name="ECM_collagen")
    used = [g for g in ECM_HUMAN if g in ad.var_names]
    print(f"  ECM module genes present: {len(used)}/{len(ECM_HUMAN)} -> {used}")

    fib = ad[ad.obs.cell_type=="Fibroblast"]
    a = fib.obs.loc[fib.obs.condition=="Normal_buttock","ECM_collagen"]
    b = fib.obs.loc[fib.obs.condition=="Photoaged_arm","ECM_collagen"]
    mw_u, mw_p = (stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan))
    rows = [
        {"metric":"fibroblast_n_total","value":int(fib.n_obs)},
        {"metric":"fibroblast_n_Normal_buttock","value":int(len(a))},
        {"metric":"fibroblast_n_Photoaged_arm","value":int(len(b))},
        {"metric":"fibroblast_ECM_mean_Normal_buttock","value":round(float(a.mean()),4)},
        {"metric":"fibroblast_ECM_mean_Photoaged_arm","value":round(float(b.mean()),4)},
        {"metric":"delta_Photoaged_minus_Normal","value":round(float(b.mean()-a.mean()),4)},
        {"metric":"cell_level_MWU_p_exploratory","value":float(f"{mw_p:.3e}") if mw_p==mw_p else np.nan},
        {"metric":"n_module_genes_present","value":len(used)},
        {"metric":"total_cells_after_QC","value":int(ad.n_obs)},
        {"metric":"biological_replication","value":"NONE (single individual; cell-level comparison is exploratory)"},
        {"metric":"condition_mapping_assumption","value":"suffix -1=GSM8478429 photoaged arm; -2=GSM8478430 normal buttock (GEO sample order)"},
    ]
    pd.DataFrame(rows).to_csv(RES/"GSE275491_public_photoaging_reanalysis.tsv", sep="\t", index=False)
    print(f"  fibroblast ECM module: Normal buttock {a.mean():+.3f} (n={len(a)}) vs "
          f"Photoaged arm {b.mean():+.3f} (n={len(b)})  delta={b.mean()-a.mean():+.3f}  "
          f"cell-level MWU p={mw_p:.2e} (EXPLORATORY, single individual)")
    print("  [written]", RES/"GSE275491_public_photoaging_reanalysis.tsv")
    return dict(normal=float(a.mean()), photoaged=float(b.mean()), delta=float(b.mean()-a.mean()),
                p=float(mw_p), n_fib=int(fib.n_obs))


if __name__ == "__main__":
    r1 = reanalyze_GSE110978()
    r2 = reanalyze_GSE275491()
    print("\n" + "#"*78)
    print("ROUND 4 REANALYSIS SUMMARY")
    print(f"  GSE110978 (replicate-aware, n=4v4): ECM module delta(old-young)={r1['delta']:+.3f} "
          f"Welch p={r1['p']:.4g}; per-gene {r1['ndown']} down/{r1['nup']} up, {r1['nsig']} p<0.05")
    print(f"  GSE275491 (exploratory, 1 indiv):   fibroblast ECM delta(photoaged-normal)={r2['delta']:+.3f} "
          f"cell-level MWU p={r2['p']:.2e}; n_fib={r2['n_fib']}")
    print("#"*78)
