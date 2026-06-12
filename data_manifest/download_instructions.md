# Public-data download instructions

This package does **not** redistribute large public raw matrices. Download them from their
public accessions into a working tree, then run the analysis scripts. Every accession is
listed in `public_data_accessions.tsv`.

## Quick start

```bash
export PROJECT_ROOT=/path/to/working_tree   # where 05_data_raw/ etc. will live
bash data_manifest/download_public_data.sh
```

## What gets downloaded

| Accession | What | Where it lands | Used by |
|---|---|---|---|
| **GSE284483** | mouse UVB bulk RNA-seq raw counts (`GSE284483_raw_counts_All_Samples.txt.gz`) | `05_data_raw/public_extension_GSE284483/` | `code/20`, `code/22` |
| **GSE173385** | mouse UVB scRNA-seq per-sample 10x matrices (Control / UV / UV+VitD) | `05_data_raw/GSE173385/` | `code/12`, `code/13`, `code/17` |
| **GSE113957** | human dermal fibroblast aging bulk RNA-seq matrix | `05_data_raw/GSE113957/` | GSE113957 outputs (provided in `results/tables/`) |
| **GSE110978** (Round 4) | mouse dermal-fibroblast aging microarray: GEO series matrix + GPL11180 probe→symbol annotation | `05_data_raw/round4_public/GSE110978/` | `code/25`, `code/26` |
| **GSE275491** (Round 4) | human photoaging scRNA-seq 10x matrices (barcodes/features/matrix) | `05_data_raw/round4_public/GSE275491/` | `code/25`, `code/26` |

Notes:
- For **GSE173385**, the public rationale (Figure 1) uses only the **Control** and **UV
  (photoaged)** conditions. The `UV+VitD` arm is the protective third arm and is not used
  in the public rationale.
- For **GSE284483**, only the **NoUVB vs UVB-saline** contrast is used; the collagen-mRNA
  treatment arm is excluded.
- For **GSE110978** (Round 4), the replicate-aware reanalysis (`code/25`) collapses the
  GPL11180 probes to one probe per gene and compares **n=4 young vs n=4 old** dermal
  fibroblasts; outputs are in `results/tables/GSE110978_*`.
- For **GSE275491** (Round 4), the single combined 10x matrix is split by barcode suffix
  (`-1` = photoaged arm, `-2` = normal buttock; GEO sample order). It is **one individual**,
  so the per-fibroblast comparison is reported as **exploratory** only.
- **Not auto-downloaded:** **PRJNA754272** (SRA-only ~185 GB; accession-verified, reanalysis
  deferred — see `results/validation_reports/PRJNA754272_accession_verification.md`) and
  **PXD018430** (PRIDE raw-only ~90 GB; **cite-only** — see `docs/audit/PXD018430_public_use_note.md`).
- GEO supplementary file names occasionally change. If a direct URL 404s, open the
  accession's GEO page (or its `.../suppl/` FTP directory) and download the listed file.

## Working layout expected by the scripts

```text
$PROJECT_ROOT/
  05_data_raw/
    public_extension_GSE284483/GSE284483_raw_counts_All_Samples.txt.gz
    GSE173385/<per-sample 10x matrices>
    GSE113957/<processed matrix>
    round4_public/GSE110978/GSE110978_series_matrix.txt.gz, GPL11180.annot.gz
    round4_public/GSE275491/GSE275491_{barcodes,features,matrix}.{tsv,mtx}.gz
  06_data_processed/    # intermediate (e.g. GSE173385_annotated.h5ad); also mirrored as data_processed/
  07_results/           # script outputs; published copies are in results/tables/
  08_figures/           # figure outputs; published copies are in figures/
```

The repository's `data_processed/`, `results/tables/`, and `figures/` folders hold the
already-generated public-data outputs (re-organized for readability) so the figures can be
inspected without re-downloading anything.

## Checksums

Checksums are not distributed with this minimal package because GEO redeposits can change
byte-level packaging. Verify integrity against the GEO record's reported file sizes, or
generate local checksums after download:

```bash
find "$PROJECT_ROOT/05_data_raw" -type f -exec shasum -a 256 {} \; > 05_data_raw_checksums.sha256
```
