
# BFFF pediatric viral transcriptomics: reproducibility package

This local staging repository supports the manuscript **Cross-context integration of pediatric blood and airway transcriptomes identifies a ciliary–inflammatory core and phase-dependent repair response in respiratory syncytial virus infection: a multi-cohort secondary analysis**.

## Release status

This is the public reproducibility archive at <https://github.com/echoqf-lang/bfff-pediatric-viral-transcriptomics>. The validated package was uploaded and made public on 2026-10-07; no GitHub Release or DOI has yet been created. The authors confirmed the institutional public-release check on 2026-10-07 for the original analysis code and author-generated derived data, figures, documentation, and candidate-gene lists included in this archive. The author order, affiliations, emails, repository contact, planned depositor, repository URL, and dual-license selection have been recorded. Original code is assigned MIT and author-generated derived data, figures, and documentation are assigned CC BY 4.0. ORCIDs may be added if available but are not represented when unconfirmed.

## Evidence boundary

The study is a secondary analysis of public pediatric transcriptomic datasets. No analyzed cohort included exposure to Bufei Fanggan Formula (BFFF). Database-derived candidate genes were used as a constrained search space; the analyses do not establish BFFF efficacy, target engagement, regulatory direction, clinical benefit, or a causal treatment mechanism.

## Public input data

The expression datasets are available from NCBI GEO under accession numbers GSE38900, GSE77087, GSE103842, GSE97742, GSE41374, and GSE155925. The study used processed Series Matrix data and metadata. Raw IDAT, CEL, FASTQ, and SRA archives were not analyzed and are not redistributed here. See `metadata/geo_accessions.tsv` and `metadata/data_provenance.md`.

## Repository contents

- `analysis/R/`: analysis scripts. Scripts 12–31 implement the final V6 extension, including the clean-header and larger-text submission figures; earlier scripts are retained only when they provide helpers sourced by the final workflow.
- `analysis/config/`: frozen analysis configuration.
- `analysis/tests/`: executable R contracts retained from the project.
- `preregistration/`: original plan and documented amendments.
- `metadata/`: accession information, sample manifests, checksums, source inventory, and provenance.
- `derived_data/`: frozen target/module definitions and model diagnostics.
- `results/main_tables/`: three final main tables.
- `results/supplementary_tables/`: eight final supplementary table files.
- `results/figures/`: five final vector PDF figures.
- `results/supporting_results/`: compact result tables underlying the V6 synthesis.
- `renv.lock` and `DESCRIPTION`: frozen R dependency declarations.

## Environment restoration

From the repository root, use a compatible R installation and run:

```r
install.packages("renv")
renv::restore()
```

Package availability, archived Bioconductor versions, and platform-specific system libraries may affect restoration. This local staging phase validates archive structure and traceability but does not claim a fresh end-to-end rerun from GEO.

## Analysis entry points

Run `Rscript run_analysis.R` without arguments to inspect the frozen order and input boundary. Execution requires the explicit `--execute` flag and stops before analysis if the documented processed inputs are absent.

The final analysis extension is represented by:

1. `analysis/R/12_gse77087_revised_main.R`
2. `analysis/R/13_freeze_single_gene_upgrade.R` through `analysis/R/20_finalize_upgrade_audit.R`
3. `analysis/R/21_acquire_airway_processed_geo.R` through `analysis/R/25_build_full_86_evidence_matrix.R`
4. `analysis/R/27_make_v6_bmc_figures_tables.R`
5. `analysis/R/28_validate_v6_biological_figure_inputs.R` through `analysis/R/30_validate_v6_biological_figure_outputs.R`
6. `analysis/R/28_make_v6_clean_header_figures.R` and `analysis/R/31_make_v6_larger_text_figures.R` for the final submission-ready Figure 2–4 typography and layout

Several scripts expect processed GEO inputs at paths documented in their source and in `metadata/data_provenance.md`. Downloaded expression payloads are deliberately excluded from this archive. The archived checksum and sample-manifest records identify the inputs used in the completed analysis.

## Validate this archive

Run:

```bash
python3 scripts/validate_repository.py
```

The validator checks the SHA-256 manifest, accession contract, table and figure inventories, maximum file size, excluded file types, local absolute paths, evidence-boundary language, and traceability of the 86-gene set, 16 bridge candidates, and weighted-Stouffer result.

## Licensing

See `LICENSE`, `LICENSES/`, and `LICENSE_STATUS.md`. MIT applies to original code; CC BY 4.0 applies to author-generated derived data, figures, and documentation. GEO data and third-party database records remain governed by their original terms. Full TCMSP and BATMAN-TCM exports and the institution's internal clinical document are not redistributed or relicensed.

## Citation

`CITATION.cff` records the confirmed author order, MIT software license, repository URL, and Fan Qiu as the repository contact. `.zenodo.json` is a local Zenodo metadata draft, and `metadata/deposit_contact.md` records Fan Qiu as the planned uploader. These repository roles do not alter the supplied manuscript, which currently identifies Xuehui Yang as the manuscript corresponding author. The Zenodo DOI should be added only after the final release has been inspected and authorized.
