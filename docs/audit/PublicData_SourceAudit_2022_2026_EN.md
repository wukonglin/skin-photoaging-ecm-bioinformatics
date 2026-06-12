# Public-Data Source Audit and 2022–2026 Dataset Expansion
### Round-3 deep audit of the public-only Figure 1 evidence base — 2026-06-12

This report (1) audits every public dataset, statistic, and figure claim behind the current
public-only Figure 1; (2) describes a fresh authoritative 2022–2026 dataset search; (3) summarizes the
screening; (4) lists recommended and rejected datasets with reasons; (5) gives the figure decision; and
(6) specifies safe vs. unsafe claim wording. **No manuscript or unpublished experimental content
was sent to any external service; only public accessions and generic scientific terms were queried.**

Companion files in this folder: `public_dataset_candidates_2022_2026.tsv` (full screening, 33 datasets,
traceable sources + decisions), `PublicData_FigureDecision_EN.md` (decision detail),
`FigureS_public_dataset_validation_caption_methods_results.md`, `PublicDataExpansion_validation_report.md`.

---

## 1. Audit of the current public evidence

### 1.1 GSE173385 — mouse UVB skin scRNA (Figure 1, Panel B)
| Field | Finding |
|---|---|
| Accession / identity | **GSE173385**, verified on the authoritative NCBI GEO record. Title: *"Single-cell RNA-seq of UVB-radiated Skin Reveals Landscape of Photoaging-related Inflammation and Protection by Vitamin D."* |
| Species / tissue | *Mus musculus*; whole skin (in vivo). |
| Condition labels | Control / UV(UVB) / protective third arm — **3 libraries, one sample per condition.** |
| Platform | Illumina NovaSeq 6000 (GPL24247); scRNA-seq. PMID 35577040; public May 2022. |
| Data state / access | Processed `RAW.tar` supplementary + SRA (SRP316525). Sufficient for our re-analysis (already processed locally to `06_data_processed/GSE173385_annotated.h5ad`). |
| **"Photoaging-like" label** | **Accurate and appropriately cautious.** UVB irradiation is the standard acute photoaging model; the figure says "UVB-exposed/photoaged," not "human disease," which is correct. |
| **Why the protective third arm is excluded** | Per the design brief: Figure 1 must compare normal vs photoaged only; this study uses no protective-cofactor arm, so showing it would create an out-of-scope reviewer question. Correctly dropped. |
| Fibroblast annotation | Dermal fibroblasts confirmed by canonical markers (Col1a1/Col1a2/Col3a1/Dcn/Lum/Pdgfra/Pdgfrb). Module score is biologically reasonable (fibroblasts carry the ECM/collagen program). |
| **Statistical caveat** | **Correctly stated.** One biological sample per condition → cell-level Mann–Whitney p (3.7×10⁻⁵) is exploratory only, **not** biological-replicate inference. Figure and caption both flag this. |
| Effect | Fibroblast ECM/collagen module Control 1.286 → UV 0.921, **Δ = −0.37** (reproduced from `GSE173385_fibroblast_module_scores.tsv`). |

**Verdict:** identity, labeling, annotation, and the exploratory caveat are all sound. The single weakness
is the absence of biological replication — addressed in this audit by the replicate-aware reanalyses below.

### 1.2 Enrichment of photoaging-suppressed fibroblast genes (Figure 1, Panel C)
| Field | Finding |
|---|---|
| Source gene list | Genes **suppressed in photoaged (UV) fibroblasts** from the scanpy differential-expression table (`07_results/GSE173385_fibroblast_UVdown_DE.tsv`, 21,594 ranked genes; the up-in-control / down-in-UV side). Correct directional input. |
| Tool / gene sets | Enrichr; **GO Biological Process 2021** and **Reactome 2022** (mouse symbols; Reactome shown as human reference R-HSA IDs). |
| GO vs Reactome | **Not conflated.** The two "ECM organization" rows are the independent GO term (GO:0030198, adj p 6.27×10⁻¹⁴, 31/300) and Reactome term (R-HSA-1474244, adj p 8.43×10⁻¹⁴, 30/291). Verified in `PanelC_enrichment_terms.tsv`. |
| Adjusted p | Benjamini–Hochberg adjusted p-values, interpreted as over-representation (ORA). |
| **Gene-universe caveat** | Enrichr uses its **own default background**, not a skin/fibroblast-matched universe. This is standard but should be acknowledged: ORA p-values are relative to Enrichr's background, so they index strong over-representation rather than a custom-background test. |
| Cross-species caveat | Mouse symbols mapped to human-reference Reactome IDs; conventional but a minor mapping assumption. |

**Verdict:** the enrichment is correctly built and labeled; the only additions are the Enrichr-background and
mouse→human-Reactome caveats (both already non-decisive and worth a one-line methods note).

### 1.3 Human 129-gene skin aging/photoaging ECM signature (Figure 1, Panel D)
| Field | Finding |
|---|---|
| Source basis | Literature-curated from **published human skin aging/photoaging studies**: McCabe et al. 2020 (Matrix Biology Plus; **PXD015982**), Li et al. 2021 (Front Cell Dev Biol; dermal ECM atlas), Ma et al. 2020 (Aging; **PXD018430**), Solé-Boldo et al. 2020 (Commun Biol), Zou et al. 2021 (Dev Cell), Tsitsipatis et al. (per-gene `directions` column in `public_ecm_aging_signature.tsv`). |
| Human / skin / aging relevance | **Yes** on all three — every source is human skin and aging/photoaging. |
| Matrisome categorization | Naba matrisome framework, reproduced deterministically in `04_code/17` (`matrisome_cat`). |
| **Reproducibility of counts** | **Reproduced from the local table** (128 usable genes after dropping the one "all-collagen (global)" row): collagens **10↓/1↑**, basement-membrane/laminins **11↓/0↑**, ECM glycoproteins **3↓/1↑**, proteoglycans **9↓/1↑**, elastic fiber **4↓/4↑**, ECM regulators **14↓/5↑**. Matches the figure exactly. |
| Caveat | Category-level, literature-curated (net reported direction); individual signature genes are not each re-tested here — correctly stated in the caption ("not every signature gene is individually evaluated"). |

**Verdict:** the signature is human, skin-relevant, reproducible, and honestly labeled at category level.

### 1.4 GSE113957 — already analyzed locally; an honest in-vitro bound (not currently on Figure 1)
A gold-standard human dermal-fibroblast aging dataset (n = 133, ages 1–94) was already processed locally.
It shows, by Spearman correlation vs donor age (`07_results/GSE113957_module_vs_age.tsv`): senescence/SASP
**rises** (r = +0.17, p = 0.047) but the ECM/collagen module is **flat** (r = −0.07, p = 0.42, n.s.) and
the literature "down-with-aging" genes do **not** fall (r = +0.25). GSEA agrees (Public_ECM_aging_DOWN
NES = 1.38, FDR = 0.24, n.s.). **This is the key reason GSE113957 is correctly kept off Figure 1**:
cultured fibroblasts lose the in-vivo matrix-suppression phenotype. It is valuable as an honest bound and
is shown in the supplementary figure (Panel B).

### 1.5 Figure-claim audit (each public-data claim graded)
| # | Figure 1 claim | Grade | Note / wording action |
|---|---|---|---|
| 1 | Photoaging lowers the fibroblast ECM/collagen program (Δ = −0.37) | **Supported with caveat** | Keep "exploratory, cell-level, one sample/condition" — already present. |
| 2 | Photoaging-suppressed genes enrich for ECM/collagen organization (GO 6.27×10⁻¹⁴; Reactome 8.43×10⁻¹⁴) | **Fully supported** | Add one-line Enrichr-background note in methods. |
| 3 | Human signature shows structural matrisome decline (category counts) | **Fully supported** | Keep category-level wording; do not promote to gene-level. |
| 4 | Public data **nominate** ECM/collagen remodeling as a targetable axis | **Fully supported** | "nominate / design rationale" wording is correct. |
| 5 | The axis **motivates** the SFMP–HA–SF design (Panel E schematic) | **Fully supported (schematic, not result)** | Keep "design rationale, not a result." |
| — | Any claim that the material **restores/rescues/reverses** the public genes | **Not present (correctly)** | Must remain absent from Figure 1. |

**Overall:** every current Figure-1 claim is either fully supported or supported-with-caveat, and all caveats
are already on the figure. No claim needs softening; two minor methods-line caveats (Enrichr background,
mouse→human Reactome) are recommended.

---

## 2. New 2022–2026 dataset search strategy

A parallel authoritative search was run across complementary angles (human scRNA photoaging; human bulk/
microarray photoaging; mouse UVB/aging; skin-aging proteomics; spatial/senescence), using **only generic
scientific terms** and querying primary repositories and literature: **NCBI GEO, NCBI SRA/BioProject, EBI
PRIDE/ProteomeXchange, PubMed/PMC, ArrayExpress/BioStudies, GSA/NGDC.** Every candidate accession was then
**independently verified against its authoritative repository page**; candidates that could not be confirmed
(404, "not in database," or a PubMed ID with no deposited dataset) were flagged as unverified and excluded
from numeric use. **36 candidates were surfaced; 28 verified; 8 unverified.** The full result with per-
accession metadata, traceable URLs, and decisions is in `public_dataset_candidates_2022_2026.tsv`
(33 curated rows after consolidation).

---

## 3. Dataset screening summary

| Decision | Count | Meaning |
|---|---|---|
| `use_in_supplement` | 5 | On-axis in-vivo datasets worth reanalysis for the supplement (one already reanalyzed here). |
| `cite_only` | 16 | Supportive references (incl. the proteomic signature sources and the in-vitro bound). |
| `needs_manual_download_or_permission` | 4 | Cited but not currently retrievable from the repository, or access-gated/not-yet-public. |
| `screened_out` | 8 | Off-axis (vascular, epidermal, immune/psoriasis, wound/KO, in-vitro perturbation) or mislabeled. |

**Important audit correction:** three proteomics accessions that an earlier internal shortlist listed as
"verified via source paper" are **not currently retrievable from PRIDE** — **PXD016440** (dermal ECM atlas,
Li 2021) and **PXD021194** (peptide-fingerprint photoaging) return "not in database"/404, and **PXD045887**
does not exist in PRIDE; **PXD050746** is still under review. The **published tables** of these studies
remain valid literature support, but the **raw datasets cannot be used for new numeric claims** until the
deposits are accessible. This is reflected in their decisions.

---

## 4. Recommended datasets (and why)

1. **GSE284483** — mouse UVB photoaging, in-vivo dorsal skin, **bulk RNA-seq, n = 3/group** (2025; PMID
   41233777). *Reanalyzed in this audit* (no-UVB vs UVB-saline; treatment arm excluded). Provides the
   **replicate-aware in-vivo** evidence the single-sample GSE173385 lacks: basement-membrane collagens
   (*Col4a1/2*) and proteoglycans (*Dcn*, *Lum*) down, elastin up (solar elastosis); whole-skin module
   change n.s. (fibroblast signal diluted in bulk). → **Supplement (done).**
2. **GSE275491** — human sun-exposed vs unexposed skin scRNA (2024; PMID 39540047). Collagen-degradation +
   ECM-receptor interaction enriched in sun-exposed fibroblasts; all fibroblast subclusters decreased.
   Strongest **human in-vivo photoaging** single-cell match found; single individual (no replication).
   → **Supplement (download + fibroblast-module reanalysis recommended).**
3. **PRJNA754272** — human forearm skin scRNA, young (n = 3) vs old (n = 4), ~50,000 fibroblasts (2022;
   PMID 35069694). **Replicate-aware human in-vivo** aging with large fibroblast capture (intrinsic aging).
   → **Supplement (recommended reanalysis).**
4. **GSE130973 / GSE274955** — human in-vivo aging (intrinsic; PMID linked) and human photoexposed-vs-
   photoprotected scRNA (prior-verified). Strong human in-vivo references; GSE274955 needs GEO re-verify.
5. **GSE110978** — mouse aged dermal fibroblasts, **in-vivo, n = 4/group**, with an explicit "old fibroblasts
   reduce ECM-formation genes" finding (Salzer et al., Cell 2018; PMID 30415840). The cleanest replicate-
   aware **fibroblast-resolved** confirmation; intrinsic aging, microarray. → recommended next reanalysis.
6. **PXD015982 / PXD018430** — human in-vivo skin aging/photoaging proteomics; **the proteomic basis of the
   129-gene matrisome signature already in Figure 1D**. → cite.
7. **GSE117763** — mouse fibroblasts fresh vs cultured (PMID 30415840): directly demonstrates that **culture
   distorts the age-related ECM signature**, corroborating the in-vitro bound. → cite (supports caveat).

---

## 5. Datasets rejected (and why)

- **GSE196395** (adventitial **aorta** fibroblasts) — vascular, not skin. *screened_out.*
- **GSE283557** (TSP2-knockout dermal fibroblasts, wound healing) — genotype/wound, no normal-vs-aged
  baseline. *screened_out.*
- **GSE137176** (epidermal stem cells) — wrong cell type. *screened_out.*
- **GSE250390** (skin-digestion method study) / **GSE151177** (psoriasis emigrating cells) — immune/method
  focus, not stromal ECM. *screened_out.*
- **GSE311942** (cultured HDF + tRF-34 under UVA1) — in-vitro perturbation, not a tissue normal-vs-photoaged
  contrast. *screened_out.*
- **GSE175011** — accession **mismatch** (an ENCODE vascular sample; the claimed skin-senescence paper
  actually reused other public datasets). Flagged as hallucination-adjacent and excluded. *screened_out.*
- **PXD045887** (does not exist in PRIDE), **PXD016440 / PXD021194** (cited but not retrievable),
  **PXD050746** (under review), **PMID:39199288 / 39370688 / 41293818** (paper IDs, no deposited dataset) —
  not usable for numeric claims; *needs_manual_download_or_permission* or *cite_only*.

---

## 6. Figure decision (summary; full rationale in `PublicData_FigureDecision_EN.md`)

**Keep the current Round-3 public-only Figure 1 as the recommended front figure; do NOT build a Figure 1 v2.**
New datasets do not de-crowd Figure 1 and the strongest human in-vivo additions require download/cell-level
reanalysis (not faked here). Instead, the new evidence is delivered as a **supplementary public-data
validation figure** (`FigureS_public_dataset_validation`) that triangulates the rationale, adds the
**replicate-aware GSE284483 reanalysis**, and honestly bounds the claim. This keeps Figure 1 a clean
public-data pre-rationale (the design brief's requirement) while materially raising reviewer confidence.

---

## 7. Safe claim wording (use these)

- "Public skin aging/photoaging datasets **nominate** ECM/collagen remodeling as a disease-relevant,
  **targetable** repair program."
- "Photoaging/aging **suppress the fibroblast ECM/collagen program** in vivo."
- "Public data show **structural matrisome decline** (collagen, basement-membrane/laminin, proteoglycan)."
- "An independent **replicate-aware** in-vivo dataset shows photoaging-associated **matrisome remodeling**
  (basement-membrane and proteoglycan decline, elastin increase)."
- "These public data **motivate the SFMP–HA–SF ECM-mimetic adhesive hydrogel design** (design rationale)."
- "to be **tested/validated in this study**."

## 8. Claim wording to avoid (or soften)

- "**proves / confirms the therapeutic target**" → use "nominate / support."
- "**restores / rescues / reverses aging**" (for any public-data gene) → not on a public-data figure.
- "the material **re-engages / corrects** the public-nominated genes" → belongs only to later validation,
  framed as pathway-level support, never gene-level rescue.
- "**every** signature gene declines / a **uniform** ECM loss" → use "predominant/structural decline" and,
  for the replicate-aware mouse data, "remodeling (BM/proteoglycan down, elastin up)," not "uniform loss."
- Citing **PXD016440 / PXD021194 / PXD045887 / PXD050746** as available raw data → cite the published
  tables only until the deposits are retrievable.
