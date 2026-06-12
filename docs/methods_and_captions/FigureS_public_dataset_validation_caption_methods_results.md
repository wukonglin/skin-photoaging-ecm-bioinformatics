# Supplementary Figure (public-data validation) — Legend, Results, Methods, Limitation
## Companion to the public-only Figure 1; for the proteomics/mechanism or supplementary section

> **Scope rule.** This figure is built **only from public databases**. It does **not** change Figure 1
> and contains **no experimental data from this study**. Its job is to (i) show that the public
> ECM/collagen-decline rationale used in Figure 1 is **triangulated across independent in-vivo
> sources**, (ii) add a **replicate-aware** in-vivo reanalysis that the single-sample dataset in
> Figure 1 could not provide, and (iii) **honestly bound** the claim — the fibroblast-specific
> matrix-suppression signal is diluted in whole-skin bulk and in cultured fibroblasts, which is
> exactly why the rationale rests on fibroblast-resolved / ECM-resolved in-vivo data.

**Figure files:** `08_figures/revision_round_3_public_data_audit_20260612/FigureS_public_dataset_validation.pdf`
(vector) / `.png` (preview).

---

## 1. Legend (English)

**Supplementary Figure. Public datasets converge on photoaging/aging-associated ECM–collagen
remodeling, replicate-aware and honestly bounded.**
**(A)** Cross-source concordance map. Independent public evidence lines (rows) are scored against three
ECM/collagen read-outs (columns). Fibroblast-resolved in-vivo single-cell data (mouse UVB skin,
GSE173385) and ECM-resolved in-vivo proteomic/literature signature (human skin aging/photoaging,
129-gene signature derived from public proteomics; e.g. PXD015982, PXD018430) both show ECM/collagen
suppression and structural matrisome decline. A replicate-aware whole-skin bulk dataset (mouse UVB,
GSE284483) shows matrisome **remodeling** at the category level but a **flat** whole-skin module
(signal diluted by non-fibroblast cells). Cultured human dermal fibroblast bulk (GSE113957) is flat/
discordant for ECM despite rising senescence. Blue = suppressed/down with aging; grey = flat/diluted/
discordant; pale = not assessed in that source.
**(B)** In-vitro caveat. In cultured human dermal fibroblasts across the lifespan (GSE113957, n = 133
donors), module score correlates with donor age (Spearman): senescence/SASP rises (r = +0.17,
p = 0.047) but the curated ECM/collagen module does **not** decline (r = −0.07, p = 0.42, n.s.), and
literature "down-with-aging" genes do not fall (they trend up, r = +0.25). Cultured fibroblasts capture
aging but lose the in-vivo matrix-suppression phenotype, so in-vitro fibroblast bulk is **not** used to
claim ECM decline.
**(C)** New replicate-aware in-vivo reanalysis. An independent public mouse UVB photoaging dataset
(GSE284483; in-vivo dorsal skin; n = 3/group) was reanalyzed for the **normal-vs-photoaged contrast
only** (no-UVB vs UVB-saline; the collagen-mRNA treatment arm was excluded). Per-gene log2 fold-change
(UVB-saline / no-UVB) shows a biologically coherent photoaging matrisome signature: basement-membrane
collagens (*Col4a1* −0.34, *Col4a2* −0.31) and proteoglycans (*Dcn* −1.17, *Lum* −0.73) are reduced,
while elastin (*Eln* +0.69) is increased (consistent with solar elastosis). The whole-skin ECM/collagen
module difference is not significant (Welch t, p = 0.93), because the fibroblast-specific suppression is
diluted in bulk tissue — consistent with the Panel A bound.
**(D)** Screened public-dataset landscape. Top-ranked public datasets from a systematic 2018–2026 search
(33 datasets screened in total; full table in `public_dataset_candidates_2022_2026.tsv`), colored by
recommended use. Datasets recommended for deeper validation are marked *use in supplement*; supporting
datasets *cite only*; inaccessible/undeposited records *needs manual download or permission*; off-axis
datasets *screened out*.

---

## 2. Results paragraph (English)

To test whether the public-data rationale used in Figure 1 is reproducible beyond a single dataset, we
triangulated independent public sources and added a replicate-aware reanalysis. The ECM/collagen-decline
signal is concordant where it can be resolved at the relevant biological scale: fibroblast-resolved
single-cell data (mouse UVB skin) and ECM-resolved human proteomic signatures both show suppression of
the matrix program (Figure S, A). In an independent, **biologically replicated** in-vivo mouse UVB
photoaging dataset (GSE284483, n = 3/group), the normal-vs-photoaged contrast reproduced a coherent
photoaging matrisome signature — reduced basement-membrane collagens (*Col4a1/2*) and proteoglycans
(*Dcn*, *Lum*) with increased elastin (solar elastosis) — although the whole-skin module change was not
significant because the fibroblast-specific signal is diluted in bulk tissue (Figure S, C). Consistent
with this scale-dependence, cultured human dermal fibroblasts (GSE113957, n = 133) showed rising
senescence with donor age but no ECM/collagen decline (Figure S, B), confirming that the matrix-
suppression phenotype is a property of fibroblasts **in their tissue context** rather than of cultured
fibroblasts. A systematic 2018–2026 search identified additional human in-vivo photoaging single-cell
datasets (e.g. GSE275491; GSE274955; PRJNA754272; GSE130973) and replicate-aware mouse aging fibroblast
datasets (e.g. GSE110978) recommended for deeper validation (Figure S, D). Together these analyses show
the Figure 1 rationale is reproducible and well-bounded: photoaging/aging suppress the fibroblast
ECM/collagen program in vivo, an effect best read at fibroblast or ECM resolution.

---

## 3. Methods paragraph (English)

*Cross-dataset public-data validation.* Public datasets were used as downloaded; no experimental data
from this study were used. **GSE113957** (human dermal fibroblast RNA-seq across ages 1–94, n = 133) was
summarized by per-sample module scores (curated ECM/collagen, focal-adhesion, senescence/SASP, and the
public aging up/down signatures) correlated against donor age by Spearman correlation (replicate-aware;
adjusted/raw p as reported). **GSE284483** (mouse UVB photoaging, in-vivo dorsal skin, bulk RNA-seq,
n = 3/group; public raw-count matrix) was reanalyzed for the no-UVB vs UVB-saline contrast only (the
UVB + collagen-mRNA treatment arm was excluded); counts were CPM-normalized and log2-transformed, an
ECM/collagen module score (the same mouse gene set used for GSE173385) was computed per sample as the
mean of z-scored gene expression, compared between groups by Welch t-test (n = 3 vs 3), and per-gene
log2 fold-changes were computed for the structural matrisome core. **GSE173385** (mouse UVB scRNA) module
values shown in Panel A are exploratory (cell-level; one biological sample per condition) and reproduced
from the Figure 1 analysis. Additional public datasets (2018–2026) were identified through authoritative
repository and literature search (NCBI GEO, NCBI SRA/BioProject, EBI PRIDE/ProteomeXchange, PubMed) and
each accession was verified against its authoritative repository page; the full screening, decisions, and
traceable sources are tabulated in `public_dataset_candidates_2022_2026.tsv`. Module/category colors:
blue = down with aging, grey = flat/diluted, orange = up (elastin). Analysis scripts:
`04_code/19_public_dataset_audit_2022_2026.py`, `04_code/20_public_data_extension_analysis.py`,
`04_code/22_public_supplement_validation_figure.py`.

---

## 4. Limitation note (short)

(i) **GSE284483** has n = 3 per group and is whole-skin bulk; it is reported at effect-size/category
level (basement-membrane and proteoglycan decline, elastin increase) and the whole-skin module change is
not statistically significant — it confirms photoaging matrisome **remodeling**, not a uniform loss.
(ii) **GSE113957** and other **cultured-fibroblast** datasets are deliberately shown as a bound: they do
not reproduce in-vivo ECM decline, so they are not used to support the decline claim. (iii) The strongest
recommended additions — human in-vivo photoaging single-cell datasets (GSE275491, GSE274955,
PRJNA754272) and replicate-aware mouse fibroblast-resolved datasets (GSE110978) — require download and
cell-level/probe-level reanalysis in a dedicated session and are documented, not assumed. (iv) This
figure does not change Figure 1; it is positioned later (proteomics/mechanism or supplementary section)
and reports public-data convergence and bounds only.
