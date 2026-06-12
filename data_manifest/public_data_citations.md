# Public-data citations

All datasets below are public. Cite them by their repository accession; the authoritative
citation (including the associated publication) resolves on each accession's repository
record. Verified PMIDs are listed where available; where a PMID was not independently
verified, the GEO/repository accession record is given as the authoritative source rather
than a guessed identifier.

## Datasets used in this package

### GSE173385 — mouse UVB skin, single-cell RNA-seq
- Repository record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173385
- Dataset title (as deposited): *"Single-cell RNA sequencing of UVB-radiated skin reveals
  landscape of photoaging-related inflammation and protection by vitamin D."*
- Note on use: only the **dermal fibroblast ECM/collagen** axis is re-analyzed here; the
  original study's inflammation / vitamin-D focus is not part of this package. One sample
  per condition — treated as exploratory / cell-level.
- Associated publication: see the linked reference on the GEO record above.

### GSE113957 — human dermal fibroblasts, bulk RNA-seq (aging)
- Repository record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113957
- Dataset title (as deposited): *"Predicting age from the transcriptome of human dermal
  fibroblasts."*
- Note on use: cultured-fibroblast aging series (n = 133 donors) used as an honest
  in-vitro/bulk bound. Associated publication: see the linked reference on the GEO record.

### GSE284483 — mouse UVB skin, bulk RNA-seq (replicated)
- Repository record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284483
- Associated publication (verified): PMID 41233777 — https://pubmed.ncbi.nlm.nih.gov/41233777/
- Note on use: NoUVB vs UVB-saline contrast only (n = 3/group); the collagen-mRNA treatment
  arm is excluded. Replicate-aware reanalysis.

### GSE110978 — mouse dermal fibroblasts, aging microarray (Round 4)
- Repository record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110978
- Platform: Affymetrix Mouse Genome 430 2.0 (PM), GPL11180.
- Associated publication (verified): PMID 30415840 — https://pubmed.ncbi.nlm.nih.gov/30415840/
  (Salzer et al., *Cell* 2018).
- Note on use: replicate-aware reanalysis (n = 4 young vs n = 4 old dermal fibroblasts);
  ECM/collagen module declines with aging (Welch p = 1.2×10⁻⁵; 17/19 genes down). First
  in-vivo, fibroblast-resolved, replicate-aware ECM read-out in this package (intrinsic aging).

### GSE275491 — human skin photoaging, single-cell RNA-seq (Round 4)
- Repository record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275491
- Associated publication (verified): PMID 39540047 — https://pubmed.ncbi.nlm.nih.gov/39540047/
- Note on use: paired photoaged arm vs normal buttock within **one individual**, so the
  per-fibroblast ECM comparison is reported as **exploratory** only (not used to claim ECM decline).

### PRJNA754272 — human skin aging, single-cell RNA-seq (Round 4; accession-verified, deferred)
- Repository record: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA754272
- Associated publication (verified): PMID 35069694 — https://pubmed.ncbi.nlm.nih.gov/35069694/
- Note on use: SRA-only (~185 GB FASTQ; 7 donors); accession-verified, reanalysis deferred to
  future donor-level validation. See `results/validation_reports/PRJNA754272_accession_verification.md`.

### PXD018430 — human skin aging proteomics (Round 4; cite-only)
- Repository record: https://www.ebi.ac.uk/pride/archive/projects/PXD018430
- Associated publication (verified): PMID 32602849 — https://pubmed.ncbi.nlm.nih.gov/32602849/
  (Ma et al., *Aging (Albany NY)* 2020; DOI 10.18632/aging.103461).
- Note on use: raw-only PRIDE deposit (~90 GB; no processed protein table) → **cite-only**;
  no new numbers are derived here. See `docs/audit/PXD018430_public_use_note.md`.

## Literature-curated signature

### Human matrisome aging/photoaging signature (129 genes)
- File: `data_processed/public_ecm_aging_signature.tsv` (gene-level; each gene tagged with
  the public study/studies and direction supporting it).
- Compiled from published, public human skin aging/photoaging studies. Per-gene provenance
  is recorded in the file's `directions` column.

## Full systematic dataset screen

The documented 2018–2026 screen of public skin-aging datasets, with one traceable source
URL per accession (and verified PMIDs where available), is in:
`data_manifest/public_dataset_candidates_2022_2026.tsv`
(33 candidates across GEO / SRA / BioProject / PRIDE; see `docs/audit/` for the screening
rationale and decisions). Inaccessible or unverifiable accessions are flagged in that table
and are not used for any numeric claim.

## How to cite

- Cite each **public dataset** by its accession and repository (and the linked publication
  on the accession record).
- If this **repository** is later made public, also cite it via `CITATION.cff`.
