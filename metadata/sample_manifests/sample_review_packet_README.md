# Task 4 independent sample review

`sample_review_packet.tsv` is an unsigned review packet, not a frozen manifest.

The independent reviewer must inspect every row against the raw evidence columns and edit only the seven review fields below. Do not edit machine evidence or `packet_schema_version`. For all 1,031 rows, fill:

- `review_status`: `completed_independent_review`
- `reviewer`: a stable reviewer identifier other than `sample_audit_analyst`
- `reviewed_at`: ISO-8601 timestamp such as `2026-08-05T12:00:00+0800`
- `review_decision`: `agree`, `disagree`, or `unresolved`
- `final_include`: `TRUE` or `FALSE`
- `final_exclusion_reason`: required whenever `final_include=FALSE`
- `review_notes`: concise rationale, especially for discrepancies

For `agree`, `final_include` and any exclusion reason must reproduce the machine decision. A `disagree` or `unresolved` row is resolved conservatively as excluded; the validator will reject conversion of a machine exclusion into an inclusion.

The seven GSE188427 arm conflicts retain RSV identity and Day 1 eligibility where applicable, but their inpatient/outpatient severity is unresolved and they are never eligible for severity-stratified analysis.

After the completed packet is saved at the same path, rerun `analysis/R/03_audit_samples.R`. The script preserves the review fields, verifies that machine evidence is byte-stable, checks all 1,031 signatures, and only then writes `sample_manifest_frozen.tsv`.
