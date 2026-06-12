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
