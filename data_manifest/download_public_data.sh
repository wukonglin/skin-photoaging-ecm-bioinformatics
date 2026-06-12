#!/usr/bin/env bash
# =============================================================================
# download_public_data.sh
# Download the PUBLIC raw inputs used by this package from their public
# accessions. No private/manuscript data is fetched. Files land under
# $PROJECT_ROOT/05_data_raw (set PROJECT_ROOT, or run from your working tree).
# =============================================================================
set -euo pipefail
ROOT="${PROJECT_ROOT:-$(pwd)}"
RAW="$ROOT/05_data_raw"
mkdir -p "$RAW"
echo "[download] target: $RAW"

# --- GSE284483 : mouse UVB skin, bulk RNA-seq, n=3/group (replicated) --------
# Verified supplementary file used by code/20_public_data_extension_analysis.py
mkdir -p "$RAW/public_extension_GSE284483"
curl -fL --retry 3 -o "$RAW/public_extension_GSE284483/GSE284483_raw_counts_All_Samples.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284483/suppl/GSE284483_raw_counts_All_Samples.txt.gz" \
  || echo "[warn] GSE284483: if this 404s, list https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284483/suppl/ and adjust the filename."

# --- GSE173385 : mouse UVB skin, scRNA-seq -----------------------------------
# Per-sample 10x matrices (Control / UV / UV+VitD) are listed on the GEO record.
# code/12_scRNA_UVB_GSE173385.py expects them under 05_data_raw/GSE173385/.
mkdir -p "$RAW/GSE173385"
echo "[manual] GSE173385: download the per-sample matrices from"
echo "         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173385"
echo "         into $RAW/GSE173385/  (the public rationale uses the Control and UV conditions)."

# --- GSE113957 : human dermal fibroblast aging, bulk RNA-seq -----------------
mkdir -p "$RAW/GSE113957"
echo "[manual] GSE113957: download the processed expression matrix from"
echo "         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113957"
echo "         into $RAW/GSE113957/."

# --- GSE110978 (Round 4) : mouse dermal-fibroblast aging microarray (n=4 vs 4) ----
# Direct-downloadable: GEO series matrix + GPL11180 probe->symbol annotation.
# code/25_public_dataset_reanalysis_round4.py expects these under round4_public/GSE110978/.
mkdir -p "$RAW/round4_public/GSE110978"
curl -fL --retry 3 -o "$RAW/round4_public/GSE110978/GSE110978_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110978/matrix/GSE110978_series_matrix.txt.gz" \
  || echo "[warn] GSE110978 series matrix: check https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110978/matrix/"
curl -fL --retry 3 -o "$RAW/round4_public/GSE110978/GPL11180.annot.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11180/annot/GPL11180.annot.gz" \
  || echo "[warn] GPL11180 annotation: check https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11180/annot/"

# --- GSE275491 (Round 4) : human photoaging scRNA-seq, 10x (1 individual) ---------
# Direct-downloadable combined 10x matrix (barcodes/features/matrix).
# code/25_public_dataset_reanalysis_round4.py expects these under round4_public/GSE275491/.
mkdir -p "$RAW/round4_public/GSE275491"
for trio in "barcodes:tsv" "features:tsv" "matrix:mtx"; do
  name="${trio%%:*}"; ext="${trio##*:}"
  curl -fL --retry 3 -o "$RAW/round4_public/GSE275491/GSE275491_${name}.${ext}.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275491/suppl/GSE275491_${name}.${ext}.gz" \
    || echo "[warn] GSE275491 ${name}: check https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275491/suppl/"
done

# --- PRJNA754272 / PXD018430 : NOT auto-downloaded ---------------------------
echo "[skip] PRJNA754272: SRA-only (~185 GB FASTQ), reanalysis deferred — see results/validation_reports/PRJNA754272_accession_verification.md"
echo "[skip] PXD018430: PRIDE raw-only (~90 GB), cite-only — see docs/audit/PXD018430_public_use_note.md"

echo "[done] See data_manifest/download_instructions.md for manual steps and checksums."
