# Skin Photoaging ECM/Collagen Bioinformatics Reproducibility Package

A minimal, public-data reproducibility package for bioinformatics analyses of
**extracellular-matrix (ECM) / collagen remodeling in skin photoaging and aging**.
It contains the analysis code, public-dataset manifests, small processed public-data
tables, and the generated figures needed to reproduce a public-data rationale for an
ECM/collagen-centered repair program, plus a documented systematic screen of public
skin-aging datasets.

This repository is intentionally **scoped to public data and figure-generation code**.
It is a generalized reproducibility package, not a manuscript.

---

## 1. Scope

- Public single-cell and bulk transcriptomic analyses of fibroblast ECM/collagen
  programs in UVB photoaging and chronological skin aging.
- A curated, literature-derived human matrisome aging signature (gene-level, public).
- A documented 2018–2026 systematic screen of public skin-aging datasets
  (GEO / SRA / PRIDE / BioProject), with per-accession verification.
- Figure-generation code and the resulting publication figures.

## 2. Scope and confidentiality statement

This repository contains **public-data analyses and figure-generation code only**. It
**intentionally excludes** unpublished manuscript files, private raw experimental inputs,
internal notes, and any private collaboration material; only public datasets (referenced by
accession), public-data-derived processed tables, the analysis code, and the generated
public-data figures are included. See
[`excluded_private_files_report.md`](excluded_private_files_report.md) for what was left out.

## 3. Repository contents

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .github/
│   ├── CODEOWNERS
│   └── pull_request_template.md
├── code/                         # analysis & figure-generation scripts (+ env files)
├── data_manifest/                # public-dataset accessions, citations, download steps
├── data_processed/               # small public-data-derived processed tables
├── results/
│   ├── tables/                   # public-data analysis outputs (TSV)
│   └── validation_reports/       # public-data validation report
├── figures/
│   ├── figure1_public_ecm_rationale/        # public-data rationale figure
│   └── supplement_public_data_validation/   # public-data validation supplement
├── docs/
│   ├── methods_and_captions/     # figure legends / methods / results (public data)
│   └── audit/                    # public-data source audit & figure-decision notes
└── excluded_private_files_report.md
```

## 4. Figure overview

1. **Public-data Figure 1 — ECM/collagen rationale**
   (`figures/figure1_public_ecm_rationale/`). Public skin photoaging/aging data nominate
   fibroblast ECM/collagen remodeling as a disease-relevant, targetable axis. Built from
   public data only.
2. **Supplementary public-data validation**
   (`figures/supplement_public_data_validation/`). Cross-source concordance, an honest
   in-vitro/bulk bound, a replicate-aware in-vivo mouse UVB reanalysis (GSE284483), and
   the landscape of the systematic dataset screen.
## 5. How to reproduce

### 5.1 Environment

```bash
# Conda (recommended):
conda env create -f code/environment.yml      # creates env "skin-ecm-bioinformatics"
conda activate skin-ecm-bioinformatics
# or pip into a fresh venv:
python -m venv .venv && source .venv/bin/activate && pip install -r code/requirements.txt
```

### 5.2 Public-data download

Public raw inputs are **not** redistributed here; download them from their public
accessions:

```bash
bash data_manifest/download_public_data.sh    # fetches GEO processed matrices
```

See [`data_manifest/download_instructions.md`](data_manifest/download_instructions.md)
for manual steps and [`data_manifest/public_data_accessions.tsv`](data_manifest/public_data_accessions.tsv)
for every accession used.

### 5.3 Working layout and run order

The scripts resolve paths from a project root (env var `PROJECT_ROOT`, defaulting to the
script's parent directory's parent) and expect a working tree with numbered folders
(`05_data_raw/`, `06_data_processed/`, `07_results/`, `08_figures/`). For readability the
**published** copies of those outputs are reorganized here as `data_processed/`,
`results/tables/`, and `figures/`. To re-run end-to-end, set `PROJECT_ROOT` to a working
directory laid out with those numbered folders (download script targets `05_data_raw/`),
then:

```bash
export PROJECT_ROOT=/path/to/working_tree
python code/12_scRNA_UVB_GSE173385.py            # GSE173385 fibroblast module scores
python code/13_panelC_enrichment.py              # enrichment of photoaging-down fibroblast genes
python code/17_figure_round3_public_rationale.py # Figure 1 (public-data rationale)
python code/20_public_data_extension_analysis.py # GSE284483 replicate-aware reanalysis
python code/19_public_dataset_audit_2022_2026.py # regenerate the dataset candidate table
python code/22_public_supplement_validation_figure.py  # public-data validation supplement
```

(GSE113957 module-vs-age outputs in `results/tables/` are public-data analysis products;
the script that generated them is not required to rebuild the figures.)

## 6. Data availability

- **Public data** are available from their original public accessions
  (GEO / SRA / BioProject / PRIDE); see `data_manifest/`. They are redistributed here only
  as small processed/derived tables, or not at all.
- **Unpublished manuscript files and private raw experimental inputs are excluded**
  from this repository during manuscript preparation.

## 7. Citation

If you use this code or the public-data outputs, please cite the original public datasets
by their accessions (see `data_manifest/public_data_citations.md`) and, if this repository
is later made public, this repository (see `CITATION.cff`).

## 8. Branch protection / collaboration policy

- `main` is the default branch. Changes to `main` should be made via **pull request**
  and approved by the repository owner (**@wukonglin**) before merge.
- Code ownership is declared in `.github/CODEOWNERS`.
- Branch protection is enabled on `main`: pull-request review with code-owner approval is
  required before merge, and force-pushes and branch deletion are disabled.
