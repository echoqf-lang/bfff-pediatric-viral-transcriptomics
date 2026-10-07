# Single-gene upgrade deviations

- Initialization: no deviations.

## 2026-08-05: non-interactive renv activation workaround

- Symptom: `Rscript` invoked from the project root silently exited after `.Rprofile` sourced `renv/activate.R`; the requested script body did not run.
- Resolution: all phase-1 commands used `Rscript --vanilla` with `R_LIBS_USER` explicitly set to the existing project-local R 4.6 renv library.
- Scope: no package was installed or downloaded, and package versions remained those already present in the project renv library.

## 2026-08-05: post-result diagnostic for small-k HKSJ interval contraction

- Trigger: the first frozen three-cohort REML/HKSJ run produced five BH-significant genes, all with tau2=0 and unusually narrow intervals relative to the cohort-specific standard errors.
- Root cause: for highly concordant effects, unmodified Knapp-Hartung allowed the residual scale factor to fall below one; this is a known small-k behavior, not a mapping or data error.
- Resolution: retain the frozen unmodified HKSJ result as the primary exploratory analysis and add metafor `test="adhoc"` as a clearly post-result conservative diagnostic that prevents the scale factor from falling below one.
- Evidence boundary: the modified result cannot be relabeled as pre-specified and neither result is confirmatory because all three cohort outcomes had been viewed previously.

## 2026-08-05: specificity classification wording clarification

- The design memo used `broad_viral_host_response` for genes that repeat in RSV-versus-healthy comparisons but are not significant in RSV-versus-other-virus analysis.
- A non-significant contrast does not establish equivalence or a shared pan-viral effect because no equivalence margin was frozen.
- The implemented label is therefore `not_distinguished_from_other_viruses_at_current_precision`; significant same-direction and opposite-direction contrasts retain separate relative-RSV and context-dependent labels.

## 2026-08-06: fixed-module implementation detail

- The fixed BloodGen3 repertoire was frozen from the official Bioconductor 3.22 `BloodGen3Module` 1.18.0 source archive (183,912 bytes; SHA-256 recorded in the module config). No additional GEO raw data were downloaded.
- Module effect is the median member-gene log2 fold change. Competitive P values are obtained with `limma::cameraPR` using the documented default inter-gene correlation of 0.01, followed by BH correction within each cohort across evaluable modules.
- Modules with less than 50% observed membership remain in the output with `evaluable=FALSE`; they are not silently removed or assigned a P value.
- This is a post-outcome exploratory analysis. Directional module replication provides transcriptomic context and does not establish BFFF exposure or causal pathway modulation.

## 2026-08-06: cell-composition sensitivity stopped by the frozen QC gate

- The non-negative reference estimator passed the locked artificial-mixture check (Spearman rho approximately 0.999).
- No fixed, publicly auditable reference signature validated for pediatric acute-infection whole-blood microarrays was locally available. Adult LM22 and tumor RNA-seq signatures were not substituted because the population and platform reference mismatch cannot be quantified here.
- Per the pre-specified stop rule, no real-data cell fractions, adjusted group effects, or 86-gene attenuation estimates were generated. The BloodGen3 cell-associated modules remain descriptive transcriptomic context and are not called cell fractions.
