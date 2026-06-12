# PRJNA754272 — accession verification (reanalysis deferred to future work)

**Date:** 2026-06-12  ·  **Status:** accession-verified; full reanalysis **not completed this round** (SRA-only, too heavy).

## Accession and source
- **BioProject:** PRJNA754272 — https://www.ncbi.nlm.nih.gov/bioproject/PRJNA754272
- **Publication:** PMID 35069694 (2022).
- **Species / tissue:** *Homo sapiens*; skin, outer forearm (sun-protected — **intrinsic chronological aging**, not photoaging).
- **Assay:** single-cell RNA-seq (10x).
- **Sample structure:** young (n=3) vs old (n=4) female donors → **7 donors**, replicate-aware human in-vivo aging; ~50,000 fibroblasts reported.

## Data-access verification (this session)
- ENA `filereport` confirms **7 sequencing runs** (`SRR15440580`–`SRR15440586`), one per donor, paired-end 10x FASTQ.
- Read counts: ~178M–360M reads per run.
- FASTQ volume: each large run ≈ 27 GB (R1 ~6.8 GB + R2 ~20 GB); three runs ~12–13 GB. **Total ≈ 185 GB of raw FASTQ.**
- BioProject→GEO link (`elink bioproject→gds`) returns **no associated GEO series**, and **no processed cell-by-gene matrix** could be located in a public repository this session. The deposit is **raw SRA FASTQ only**.

## Why reanalysis is deferred
Completing a fibroblast ECM/collagen module comparison here would require downloading ~185 GB of FASTQ and running a full 10x alignment/quantification (Cell Ranger) on 7 donors before any module scoring — out of scope for this session, and **no shortcut processed matrix exists**. Per the project rule, **no analysis is fabricated**: this dataset is verified but not reanalyzed.

## Recommended use
- **High-value future validation:** PRJNA754272 is the one human in-vivo aging scRNA candidate with **biological replication** (7 donors) and a very large fibroblast capture — it could provide donor-level (replicate-aware) confirmation that single-individual sources (GSE173385, GSE275491) cannot.
- **Caveat:** sun-protected forearm = **intrinsic aging**, not photoaging, so it speaks to the chronological-aging arm of the ECM-decline axis rather than the UV/photoaging arm.
- **Action:** schedule a dedicated alignment run (download FASTQ → Cell Ranger → fibroblast subsetting → ECM/collagen module young-vs-old, donor-level statistics). Until then, cite as a known replicate-aware human in-vivo aging resource.

## Local output
`07_results/PRJNA754272_accession_verification.md` (this file). No numeric reanalysis table was produced (intentionally — not fabricated).
