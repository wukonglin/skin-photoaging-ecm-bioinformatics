# Figure 1 (public-data only) — Legend, Results, Methods, Limitation
## Public-data rationale figure (for the start of the paper)

> **Scope rule (Round-3 design feedback):** Figure 1 is a *pre-experiment rationale* figure built
> **only from public databases/literature**. It compares **normal vs photoaged** skin and nominates the
> ECM–collagen remodeling axis as the design target. It contains **no treated-sample proteomics result**
> and **no protective third arm**. Downstream experimental validation of this same axis is reported
> separately in the main study.
>
> Wording deliberately uses *nominate / disease-relevant / targetable / design rationale / to be validated*
> and avoids any post-experiment verb (no "shows", "re-engages", "rescues", "restored").

**Figure files:** `08_figures/revision_round_3_PI_feedback_20260612/Figure1_public_ECM_rationale.pdf` (vector,
for submission) / `.png` (preview).

---

## 1. Figure 1 legend (English)

**Figure 1. Public skin aging/photoaging datasets nominate the ECM–collagen remodeling axis as a
targetable repair program for hydrogel design.**
**(A)** Public-data rationale workflow. Public skin aging/photoaging datasets reveal fibroblast ECM/collagen
program suppression and structural matrisome decline, nominating the ECM–collagen remodeling axis as a
disease-relevant target that motivates the SFMP–HA–SF ECM-mimetic adhesive hydrogel design.
**(B)** Single-cell RNA-seq of UVB-exposed mouse skin (GSE173385; normal/Control vs photoaged/UVB).
Left, UMAP of major skin cell classes with dermal fibroblasts highlighted (annotated by
Col1a1/Col1a2/Col3a1/Dcn/Lum/Pdgfra/Pdgfrb). Middle, ECM/collagen module score across cells. Right,
fibroblast ECM/collagen module score in normal versus photoaged skin; photoaging lowers the fibroblast
ECM/collagen program (Δmean = −0.37; 0.92 photoaged vs 1.29 normal). Module-score comparisons are
exploratory (cell-level; one biological sample per condition; cell-level Mann–Whitney p = 3.7×10⁻⁵, which
is not a biological-replicate test and is reported only as an exploratory descriptor).
**(C)** Pathway enrichment of genes suppressed in photoaged fibroblasts (Enrichr; GO Biological Process 2021
and Reactome 2022; mouse gene symbols, Reactome shown as human reference R-HSA identifiers). Representative
non-redundant terms are shown (complete output in `07_results/PanelC_enrichment_terms.tsv`); the two "ECM
organization" rows are the independent GO and Reactome terms, both reaching adjusted p ≈ 6–8×10⁻¹⁴.
Suppressed fibroblast genes are over-represented for extracellular-matrix organization, collagen formation,
collagen-fibril organization, collagen biosynthesis & modification, ECM proteoglycans, elastic-fibre
formation and laminin interactions (red dashed line = adjusted p 0.05).
**(D)** A 129-gene, literature-curated human skin aging/photoaging ECM signature shows predominant
down-regulation across structural matrisome categories: collagens (10 down / 1 up), basement-membrane /
laminins (11 / 0), ECM glycoproteins (3 / 1), proteoglycans (9 / 1), elastic fibre (4 / 4) and ECM
regulators (14 / 5). Counts are shown at the category level; not every signature gene is individually
evaluated in this study.
**(E)** Design rationale (schematic, not a result). The public-nominated ECM–collagen remodeling axis
motivates an SFMP–HA–SF hydrogel designed to rebuild an ECM-mimetic adhesive interface — hyaluronic acid for
hydration and dermal filling, silk fibroin as a bioactive ECM-like polymer matrix, and silk-fibroin
microspheres to provide microscale adhesive interfaces for cells — to be validated in this study by collagen
I/III/IV, fibronectin, vinculin and histological collagen remodeling.

---

## 2. Results paragraph (English — for the start of the paper)

To define a disease-relevant tissue program that an injectable regenerative material could target, and
before generating any experimental data of our own, we first re-analyzed publicly available skin
aging/photoaging datasets. In UVB-exposed mouse skin single-cell RNA-seq (GSE173385), dermal fibroblasts
in photoaged skin showed a reduced ECM/collagen module score relative to normal skin (Δmean = −0.37;
exploratory, cell-level), indicating suppression of the fibroblast matrix-production program. Pathway
enrichment of the genes suppressed in photoaged fibroblasts identified extracellular-matrix organization as
the top term (GO adjusted p = 6.27×10⁻¹⁴; Reactome adjusted p = 8.43×10⁻¹⁴), together with collagen
formation, collagen-fibril organization, collagen biosynthesis and modification, ECM proteoglycans,
elastic-fibre formation and laminin interactions. Independently, a 129-gene literature-curated human skin
aging/photoaging ECM signature showed predominant down-regulation across structural matrisome categories,
including collagens, basement-membrane/laminins, proteoglycans and ECM regulators. Together, these public
skin aging/photoaging datasets nominate ECM/collagen remodeling as a disease-relevant and targetable repair
program in photoaged skin, providing a rationale for designing the SFMP–HA–SF hydrogel as an ECM-mimetic
adhesive interface.

---

## 3. Methods paragraph (English)

*Public-data re-analysis for target nomination.* Publicly available skin aging/photoaging datasets were used
as downloaded. UVB-exposed mouse skin single-cell RNA-seq (GSE173385) was processed in Scanpy (quality
control: >200 genes per cell, <25% mitochondrial reads; total-count normalization to 1×10⁴, log1p, 2,000
highly variable genes, PCA, Leiden clustering, UMAP). Major cell classes were annotated from canonical
markers, and dermal fibroblasts confirmed by Col1a1/Col1a2/Col3a1/Dcn/Lum/Pdgfra/Pdgfrb. An ECM/collagen
program module score was computed per cell with `score_genes`; for the rationale figure only normal (Control)
and photoaged (UVB) cells were retained, and condition comparisons within fibroblasts are reported as
exploratory (cell-level Mann–Whitney; one biological sample per condition; cell-level p-values are
descriptive, not biological-replicate significance). Genes suppressed in photoaged fibroblasts were tested
for pathway enrichment with Enrichr (GO Biological Process 2021 and Reactome 2022; mouse gene symbols;
Reactome terms shown with human reference R-HSA identifiers); a representative non-redundant set of terms is
plotted, with the complete Enrichr output tabulated in `07_results/PanelC_enrichment_terms.tsv`. A human
skin aging/photoaging ECM signature (129 genes) was curated from published studies and classified into
matrisome categories (Naba matrisome framework); per-category up/down counts reflect the net reported
direction with aging. Color code: normal blue, photoaged orange; public aging-down blue, aging-up red. No
experimental data generated in this study are shown in this figure.

---

## 4. Limitation note (short)

This figure is a public-data rationale, not proof of mechanism. (i) The mouse single-cell dataset has one
biological sample per condition, so fibroblast module-score comparisons are exploratory and reported at the
cell level only. (ii) The human ECM signature is curated from the published literature at the
category/pathway level; individual signature genes are not each evaluated experimentally here, and the figure
therefore reports structural-category trends rather than a validated gene list. (iii) Figure 1 establishes
*why* the ECM–collagen remodeling axis is a disease-relevant design target; evidence that the SFMP–HA–SF
material actually re-engages this axis is reported separately in the main study. Keeping the front figure to public data avoids implying that the design target
was defined post hoc from the study's own results.
