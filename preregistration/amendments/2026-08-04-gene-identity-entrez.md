# Pre-registration amendment: Entrez-centered gene identity

**Approved and recorded:** 2026-08-04, before any confirmatory cohort outcome was read

## Background

The original Task 2 plan used stable Ensembl gene ID as the biological primary key. Offline `org.Hs.eg.db` inspection showed that a single human Entrez gene may have multiple Ensembl aliases. Expanding those aliases as independent biological genes inflates target-set size and can also obscure identities such as CKMT1A and CKMT1B when incomplete or first-row deduplication is used.

## User-approved replacement

For Task 2, this amendment replaces only the original Ensembl-primary-key rule:

- The biological primary key is a unique human Entrez Gene ID.
- The display symbol is the official `org.Hs.eg.db` SYMBOL associated with that Entrez ID.
- All Ensembl mappings are retained as aliases and never counted as separate genes.
- Later platform probe annotation must resolve to Entrez Gene ID before target-set membership is evaluated.
- TCMSP human UniProt identifiers are mapped offline to Entrez IDs; BATMAN `Gene_ID` is interpreted and audited as Entrez ID.
- Source symbols and official symbols are both retained. Symbol conflicts are audited, not silently overwritten.
- Every allowed upstream row retains `source`, `source_row_id`, source gene identifier, raw and official symbol, herb, compound, and evidence provenance. No first-row or incomplete deduplication may discard a mapped gene identity.

The previously approved source-specific compound-key amendment remains in force. BATMAN primary evidence remains the exact, case-sensitive Score value `known target in HERB`; composite HERB labels remain full-evidence/audit only. Numeric BATMAN scores at least 0.90 enter sensitivity only, while 0.84 through 0.89 remain full-evidence only.

## Outcome boundary and scope

This amendment was approved from upstream identifier-structure inspection only. No confirmatory cohort outcome, disease-target intersection, PPI result, or expression effect informed it. It changes gene identity and counting, not the pharmacological evidence hierarchy.
