# Figure Decision — Public-Data Audit (Round 3 extension)
### 2026-06-12

## Decision

**Keep the current Round-3 public-only Figure 1 as the recommended front (rationale) figure. Do NOT
create a Figure 1 v2. Deliver the new public evidence as a separate supplementary public-data validation
figure (`FigureS_public_dataset_validation`).**

## The six decision questions (from the brief)

| # | Question | Answer | Reason |
|---|---|---|---|
| 1 | Does new public evidence **materially strengthen** the current Figure 1? | Partly — but as *support*, not as new front-figure panels | The replicate-aware in-vivo reanalysis (GSE284483) and the cross-source triangulation raise confidence, but they answer a *reviewer-robustness* question, not the *pre-experiment nomination* that Figure 1 exists to make. |
| 2 | Can it be added **without crowding** Figure 1? | No | Figure 1 already has 5 panels (A–E). the project lead explicitly asked for a cleaner, less crowded front figure; adding a 6th/7th panel (a second UMAP, a new dataset bar) works against that. |
| 3 | Does the new evidence remain **public-data only**? | Yes | All additions are public datasets; the supplement is public-only. |
| 4 | Does it **avoid the protective third arm** (no VitD)? | Yes | GSE284483's contrast is no-UVB vs UVB-saline; its collagen-mRNA treatment arm is excluded, mirroring the Figure-1 rule. |
| 5 | Does it **avoid unvalidated long gene lists**? | Yes | The supplement reports category/pathway-level results and a small structural-matrisome core (Col4a1/2, Dcn, Lum, Eln) consistent with later wet-lab validation; no long unvalidated gene list reaches the front figure. |
| 6 | Does it **improve reviewer confidence more than it adds complexity**? | Yes — **in a supplement** | Triangulation + replicate-aware reanalysis + an honest in-vitro/bulk bound is exactly what a skeptical reviewer wants, but it belongs *behind* the clean rationale, not inside it. |

Because Q2 is "no" and Q6 is "yes only in a supplement," the rule in the brief points to: **keep Figure 1,
add a supplementary public-data validation figure.**

## Why a Figure 1 v2 would *weaken* clarity

- It would re-introduce visual density the PI asked to remove.
- The strongest human in-vivo additions (GSE275491, GSE274955, PRJNA754272) require fresh download and
  cell-level reanalysis; putting them on the front figure now would either fake analysis (disallowed) or
  add an under-analyzed panel.
- The most rigorous new result (GSE284483) is honestly *bounded* (whole-skin module n.s.; signal diluted
  in bulk). That nuance is appropriate for a supplement, not for the headline rationale figure, where it
  could muddy the clean "public data nominate the axis" message.

## What makes the supplement genuinely stronger (not just more panels)

- **Triangulation (Panel A):** the ECM/collagen-decline signal is concordant across fibroblast-resolved
  in-vivo single-cell data and ECM-resolved in-vivo proteomic signatures.
- **Replicate-aware in-vivo reanalysis (Panel C):** an independent 2025 mouse UVB dataset (n = 3/group)
  reproduces a coherent photoaging matrisome signature (BM collagen + proteoglycan down, elastin up).
- **Honest bound (Panels A–C):** cultured fibroblasts (GSE113957) and whole-skin bulk do not show the
  decline, which explains *why* the rationale rests on fibroblast-/ECM-resolved in-vivo data — a
  reviewer-disarming point.
- **Systematic search (Panel D):** a documented 2018–2026 screen (33 datasets) shows the evidence base was
  surveyed, not cherry-picked, and names the best future additions.

## Files produced by this decision

- `08_figures/revision_round_3_public_data_audit_20260612/FigureS_public_dataset_validation.{pdf,png}`
- `09_outputs/revision_round_3_public_data_audit_20260612/FigureS_public_dataset_validation_caption_methods_results.md`
- `…/public_dataset_candidates_2022_2026.tsv` (screening + decisions)
- `…/PublicData_SourceAudit_2022_2026_EN.md` (audit + search + claim wording)
- Scripts: `04_code/19_public_dataset_audit_2022_2026.py`, `…/20_public_data_extension_analysis.py`,
  `…/22_public_supplement_validation_figure.py`

**Figure 1 itself is unchanged** (`08_figures/revision_round_3_PI_feedback_20260612/Figure1_public_ECM_rationale.{pdf,png}`
and its caption file remain the recommended front figure).
