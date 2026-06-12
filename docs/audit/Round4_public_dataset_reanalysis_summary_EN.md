# Round 4 — public-dataset reanalysis summary

Public-data only. Four prioritized public datasets were verified, with reanalysis completed where
feasible and honest non-completion notes where not. No results are inferred from abstracts; every
number below is computed from the downloaded public matrices. The ECM/collagen module gene set is the
same one used for the public-data rationale figure (Figure 1, Panel B): `Col1a1 Col1a2 Col3a1 Col4a1
Col4a2 Col5a1 Col6a1 Fn1 Dcn Lum Eln Fbn1 Postn Tnc Lox Mmp2 Mmp3 Mmp9 Timp1` (human orthologs for the
human dataset).

| Dataset | Type | Status | Headline result |
|---|---|---|---|
| **GSE110978** | mouse dermal-fibroblast aging microarray (n=4 vs 4) | Reanalysis completed | ECM module **down** with aging (Δ = −1.18, Welch p = 1.2×10⁻⁵); **17/19 genes down**, 12 p<0.05 |
| **GSE275491** | human in-vivo photoaging scRNA (1 individual) | Reanalysis completed (exploratory) | per-fibroblast ECM module ~comparable (1.16 photoaged vs 1.02 normal); ~40% fewer fibroblasts in photoaged arm |
| **PRJNA754272** | human in-vivo aging scRNA (7 donors) | Accession-verified; deferred | SRA-only (~185 GB FASTQ, no processed matrix) → recommended future work |
| **PXD018430** | human skin-aging proteomics (PRIDE) | Verified; cite-only | raw-only deposit (~90 GB, no processed table) → cite the paper's published tables |

## GSE110978 — replicate-aware, fibroblast-resolved (the headline)
- Source: GEO `GSE110978`; Salzer et al., *Cell* 2018 (PMID 30415840). *Mus musculus*, dermal fibroblasts (in-vivo isolated), Affymetrix Mouse 430 2.0 PM (GPL11180). n = 4 young vs 4 old.
- Processing: GEO series matrix (RMA log2); probes mapped to symbols via the GPL11180 GEO annotation (one probe per gene, highest mean expression); per-gene Welch t-test (old vs young); replicate-aware module score = per-gene z-score across samples, averaged per sample.
- Result: **ECM/collagen module young +0.59 vs old −0.59 (Δ = −1.18, Welch t = −13.7, p = 1.18×10⁻⁵)**; 17/19 genes decline with aging (12 at p<0.05; e.g. Col5a1, Col4a2, Col6a1, Col3a1, Lox, Mmp3, Eln). This is the first in-vivo, fibroblast-resolved, replicate-aware ECM read-out in the package (intrinsic aging, not UVB).
- Tables: `results/tables/GSE110978_mouse_fibroblast_aging_ECM_reanalysis.tsv`, `results/tables/GSE110978_mouse_fibroblast_aging_gene_trends.tsv`.

## GSE275491 — human in-vivo photoaging scRNA (exploratory)
- Source: GEO `GSE275491` (PMID 39540047). *Homo sapiens*, photoaged arm vs normal buttock paired within **one individual**, 10x scRNA-seq. 9,871 cells after QC; 2,285 fibroblasts.
- Result (exploratory): per-fibroblast ECM/collagen module +1.16 photoaged arm (n=858) vs +1.02 normal buttock (n=1,427), Δ = +0.14 (cell-level only; single individual → not a biological-replicate test). Photoaged arm captured ~40% fewer fibroblasts.
- Honest bound: single individual; arm-vs-buttock confounds sun exposure with body site. The per-fibroblast ECM module is not lower in photoaged skin here, so this dataset is shown as an exploratory human in-vivo cross-check and is **not** used to claim ECM decline.
- Table: `results/tables/GSE275491_public_photoaging_reanalysis.tsv`.

## PRJNA754272 — accession-verified; reanalysis deferred
- Source: BioProject `PRJNA754272` (PMID 35069694). *Homo sapiens*, sun-protected forearm (intrinsic aging), 10x scRNA-seq, young (n=3) vs old (n=4) → 7 donors; ~50,000 fibroblasts.
- Access: 7 SRA runs (SRR15440580–586), ~185 GB FASTQ, no processed matrix located. Reanalysis would need a full Cell Ranger alignment of 7 donors — out of scope this round; recommended as future donor-level (replicate-aware) validation.
- Report: `results/validation_reports/PRJNA754272_accession_verification.md`.

## PXD018430 — verified; cite-only
- Source: PRIDE `PXD018430`; Ma et al., *Aging (Albany NY)* 2020, DOI 10.18632/aging.103461 (PMID 32602849). *Homo sapiens*, volar-forearm skin, young vs elderly.
- Files: 31 files = 30 raw `.raw` (~90 GB) + 1 search binary; no processed protein-quantification table deposited. New numbers would require re-running the full MS search on the raw spectra — out of scope and unnecessary (processed results are in the paper's supplementary tables).
- Decision: cite-only — an independent human in-vivo skin-aging proteomic support line for ECM/collagen decline. Note: `PARTIAL` submission. See `docs/audit/PXD018430_public_use_note.md`.

## Supplement decision
A Round-4 supplement update is warranted because GSE110978 is a completed, reliable, replicate-aware, fibroblast-resolved reanalysis. The updated supplement (`figures/supplement_public_data_validation/FigureS_public_dataset_validation_round4`) adds GSE110978 (headline) and the GSE275491 exploratory cross-check, and makes the "replicate-aware design?" dimension explicit. Under-analyzed datasets (PRJNA754272 raw-only; PXD018430 cite-only) are not plotted as data.
