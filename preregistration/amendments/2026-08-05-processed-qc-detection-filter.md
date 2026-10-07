# Pre-registration operational amendment: processed-only QC and detectable-gene filter

**Approved and recorded:** 2026-08-05, during blinded technical QC and before any case label, case-control effect, P value, target-set enrichment result, or confirmatory decision was accessed

## Reason for this operational amendment

The frozen plan required a processed-only fallback when source-level microarray processing was unavailable, but it did not define an exact Detection-P frequency threshold for calling a gene detectable. The user subsequently restricted acquisition and analysis to directly analyzable GEO processed data and prohibited further RAW/CEL/FASTQ/SRA acquisition. This amendment fixes the previously unspecified processed-data implementation before outcome analysis.

## Frozen implementation

- For GSE105450 and GSE103842, the GEO series matrices are the expression source. Their submitter-provided processing statements consistently report Illumina BeadStudio/GenomeStudio background subtraction followed by scaling of average signal intensity (average normalization).
- Because the series-matrix intensities are already average-normalized but remain on an intensity scale, the only additional expression transformation is `log2(pmax(intensity, 1))`.
- No second quantile normalization or other between-array normalization is applied.
- The submitter-provided supplementary non-normalized processed tables are used only for their alternating `Detection Pval` columns. Their intensity columns are not used, and no source-level background correction or re-normalization is attempted.
- Probe annotation uses the frozen official NCBI GEO GPL10558 annotation and current offline `org.Hs.eg.db` Entrez mapping. Empty, multiple-Entrez, and obsolete/unresolved mappings are excluded.
- For each current unique Entrez gene, the probe with the highest mean expression across all retained samples is selected without access to phenotype labels, effects, or P values.
- A gene is called detectable when its final blindly selected probe has `Detection Pval < 0.05` in at least 10% of the QC-retained samples in that cohort.

## Blinding and threshold-selection boundary

The 10% frequency rule was the only detectability threshold selected. No alternative threshold was run or compared. Before freezing this amendment, the analyst had seen only label-free technical distributions, blind QC flags, mapping counts, and aggregate target-coverage counts. These quantities were not used to search for or optimize the threshold. No case status, case-control contrast, gene effect, model P value, target-set enrichment statistic, or confirmatory decision had been loaded or calculated.

After this amendment is frozen, the transformation and detectability rules cannot change for the confirmatory analysis on the basis of downstream results. The processed-only boundary also means that the analysis will not claim CEL image quality, bead-level quality, RNA degradation, raw-background reprocessing, or source-level normalization diagnostics.
