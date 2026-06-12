# PXD018430 — public-data verification and use decision

**Decision: cite-only.** PXD018430 is verified as a real, public human skin-aging proteomics deposit,
but it cannot supply new numeric claims for this package and is therefore cited (not reanalyzed).

## Accession (verified)
- PRIDE: https://www.ebi.ac.uk/pride/archive/projects/PXD018430 ·
  ProteomeCentral: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD018430
- Title (PRIDE API): "Quantitative analysis of skin proteomics in Chinese young and old populations by mass spectrometry"
- Species: *Homo sapiens*; instrument: Q Exactive HF; published 2020-08-05; submission type: PARTIAL.
- Paper: Ma et al., *Aging (Albany NY)* 2020, DOI 10.18632/aging.103461 (PMID 32602849). Volar-forearm skin, young vs elderly (~20 individuals).

## File-list verification (PRIDE API)
- 31 files total: **30 instrument raw files** (`*.raw`, ~1.6–3.6 GB each, ≈ 90 GB) + **1 search-engine result binary** (884 MB).
- **No processed protein-quantification table** (no `.txt`/`.tsv`/`.csv`/`.xlsx` protein matrix) is deposited — consistent with the `PARTIAL` submission type.

## Why cite-only
Deriving new protein-level numbers would require downloading ~90 GB of raw spectra and re-running the
entire MS search/quantification pipeline — out of scope for this package and unnecessary, since the
processed quantitative results are published in the article's supplementary tables. No numbers are
inferred from the abstract.

## Recommended use
Cite as an independent human in-vivo skin-aging proteomic line of support for ECM/collagen decline
(volar-forearm, ~20 donors), alongside the other human proteomic source that underpins the curated
matrisome signature. Do not add PXD018430-derived numbers to any figure or table; do not download the
raw files. This matches the prior screening decision (`cite_only`) recorded in the candidate manifest.
