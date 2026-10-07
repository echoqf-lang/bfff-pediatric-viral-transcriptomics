# Pre-registration amendment: source-specific compound identifiers

**Recorded:** 2026-08-04, before any confirmatory cohort outcome was read

## Original rule

The frozen plan required compounds to be identified by CID plus a standardized name.

## Discovery and reason

Read-only inspection of the permitted TCMSP target export found no PubChem CID field. It contains TCMSP `MOL_ID` and `molecule_ID`; treating either as a PubChem CID would create an unverifiable identifier. No online CID inference is permitted.

## User-approved amendment

Before any GSE105450 or GSE103842 outcome was read, the user approved a source-specific compound key:

- BATMAN: `PubChem_CID` plus standardized compound name.
- TCMSP: `TCMSP_MOL_ID` plus standardized compound name.

The affected output fields are `compound_id_type`, `compound_id`, and the compound component of the evidence uniqueness key. This amendment does not change the primary or sensitivity gene-evidence rules.

BATMAN primary evidence remains an exact, case-sensitive match of the Score field to `known target in HERB`. Composite labels such as `known target in HIT,HERB` and `known target in DrugBank,HIT,HERB` remain in the complete evidence and audit outputs but are not promoted to primary evidence.

## Outcome boundary

This amendment was made solely from the allowed upstream target exports and `herbs.txt`. No disease intersection, PPI, GEO effect, or confirmatory cohort outcome informed it.
