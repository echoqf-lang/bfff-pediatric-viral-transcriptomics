# Amendment: GSE77087 replaces GSE105450 in the revised main analysis

**Status:** FROZEN before downloading or loading the GSE77087 expression matrix

**Decision source:** The investigator requested removal of GSE105450 and use of GSE77087 after the original GSE105450/GSE103842 confirmatory results were known.

## Evidence-timing consequence

- This amendment does not alter or erase the original preregistered conclusion for GSE105450 plus GSE103842.
- The revised GSE77087 plus GSE103842 analysis is the manuscript's investigator-selected revised main analysis, but it is post-outcome and must not be described as the original preregistered confirmation.
- The original GSE105450 plus GSE103842 analysis remains reportable as the prior confirmatory analysis and will be moved to sensitivity/supplementary reporting rather than deleted.
- No result from GSE77087 has been downloaded, loaded, plotted, or calculated before freezing this amendment.

## Cohort role and independence boundary

- Revised main pair: GSE77087 plus GSE103842.
- GSE77087 is pediatric peripheral whole blood on GPL10558 and publicly reports 23 healthy controls and 81 acute RSV cases, including 20 outpatients and 61 inpatients before project-specific QC.
- GSE77087 and GSE105450 were recruited at Nationwide Children's Hospital during overlapping respiratory seasons by overlapping investigators. They must not enter the same analysis as independent cohorts.
- GSE105450 may appear only as a replacement/sensitivity cohort relative to GSE77087.
- GSE188427 also remains replacement-only and must not be combined with GSE105450 or GSE77087 as an independent Nationwide Children's Hospital cohort unless individual-level non-overlap is demonstrated.

## Processed-only input rule

- Permitted new downloads: official quick SOFT metadata, official supplementary file list, and `GSE77087_series_matrix.txt.gz` only.
- Prohibited new downloads: `GSE77087_RAW.tar`, CEL files, FASTQ, SRA, and redundant platform text files.
- The existing checked GPL10558 annotation is reused.

## Sample and model rule

- Include unambiguous acute RSV Day-1/diagnosis samples and healthy controls from peripheral whole blood.
- Include both outpatient and hospitalized RSV cases in the revised main contrast; hospitalized-only remains a sensitivity analysis.
- Exclude technical replicates, non-primary timepoints, samples with ambiguous infection status, and samples failing outcome-blind processed-data QC.
- The contrast is RSV minus healthy.
- Fit limma models separately by cohort.
- For GSE77087, adjust for continuous age in months, sex, and technical batch when these fields are available, identifiable, estimable, and not collinear with case status. Use complete cases without outcome-driven imputation. Any unavailable covariate is documented before target-set testing.
- Apply the same frozen target mapping, probe-selection, detectability, CAMERA, ROAST, fgsea, Stouffer, gene-level REML/HKSJ, multiplicity, and no-rescue rules used previously.

## Processed-matrix operationalization frozen before target testing

- GSE77087 does not provide Detection P values in the permitted series-matrix input. No additional non-normalized or raw archive will be downloaded to manufacture equivalence with the earlier cohorts.
- Outcome-blind sample QC uses missingness, median expression, expression IQR, median inter-sample Spearman correlation, five-PC distance, and intercept-only limma array weight. A sample is excluded only when at least two metrics exceed the previously frozen two-sided 5 scaled-MAD rule.
- After current Entrez mapping and outcome-blind highest-mean probe collapse, processed-matrix detectability is defined as positive across-sample IQR and all-sample mean at or above the cohort-wide 10th percentile of gene means. The same rule is invariant to sample order and never uses case status, effect size, or P value.
- This processed-matrix detectability rule is a documented platform-input limitation and is not claimed to reproduce Illumina Detection P filtering.

## Revised decision language

- The GSE77087 plus GSE103842 result is reported as a revised post-outcome main analysis.
- A positive result may support reproducibility within the revised pair but cannot retroactively convert the original failed preregistered H1 into a supported confirmatory result.
- A negative result is reported without further cohort replacement.
