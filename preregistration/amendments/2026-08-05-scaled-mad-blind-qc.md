# Pre-registration amendment: standard scaled-MAD blind QC rule

**Approved and recorded:** 2026-08-05, before any case label, case-control effect, model P value, target-set enrichment result, or confirmatory decision was accessed

## Verified implementation issue

The first Task 5 implementation called `stats::mad(..., constant = 1)`. That is an unscaled raw median absolute deviation, whereas the standard R definition uses the normal-consistency scaling constant 1.4826. Consequently, the earlier blind-QC exclusions, retained-sample expression objects, coverage tables, figures, and checksums from commit `5948947` are superseded and cannot be used downstream.

## Frozen replacement rule

For every pre-specified blind QC metric:

- Center is the sample median.
- Scale is `stats::mad(x, center = median(x), constant = 1.4826)`, matching the standard default scaled MAD definition.
- A sample is flagged when `abs(x - median(x)) > 5 * MAD`; the rule is two-sided.
- If MAD is zero and all finite values are equal, no sample is flagged for that metric.
- If MAD is zero but finite values differ from the median, every value with absolute deviation greater than zero is flagged.
- If any base expression, Detection-P, derived metric, or gate input is non-finite, processing stops for the entire cohort. Non-finite values are never silently retained, excluded, or ignored.
- A sample is excluded only when at least two blind metrics are flagged under these rules.

## Outcome and selection boundary

The scaled constant is selected because it is the standard definition used by `stats::mad`, not because of the number or identity of samples it retains. No alternative MAD constant or QC threshold will be compared against case labels, effects, P values, enrichment statistics, or confirmatory decisions. The correction remains entirely blinded.

All prior Task 5 artifacts based on raw MAD (`constant = 1`) are explicitly superseded. They must be regenerated, and only the regenerated scaled-MAD artifacts may enter downstream analysis.
