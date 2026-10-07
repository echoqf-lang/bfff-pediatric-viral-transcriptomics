#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 20 * 1024 * 1024
EXPECTED_GEO = {"GSE38900", "GSE77087", "GSE103842", "GSE97742", "GSE41374", "GSE155925"}
PROHIBITED_SUFFIXES = {".docx", ".doc", ".rds", ".rdata", ".cel", ".idat", ".fastq", ".fq", ".sra", ".tar", ".gz"}
PROHIBITED_NAMES = {".git", ".DS_Store", ".Rhistory", ".RData", "auth.json"}
ABSOLUTE_USER_PATH = re.compile(r"/(?:Users|home)/[^/\s]+/")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    errors = []
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.relative_to(ROOT).parts
    )
    for path in files:
        relative = path.relative_to(ROOT)
        if any(part in PROHIBITED_NAMES for part in relative.parts):
            errors.append(f"prohibited name: {relative}")
        if path.suffix.lower() in PROHIBITED_SUFFIXES:
            errors.append(f"prohibited suffix: {relative}")
        if path.stat().st_size > MAX_FILE_BYTES:
            errors.append(f"file exceeds 20 MiB: {relative}")
        if "GSE105450" in path.name and "analysis" not in relative.parts:
            errors.append(f"deprecated result path: {relative}")
        if path.suffix.lower() in {".md", ".txt", ".tsv", ".csv", ".yml", ".yaml", ".r", ".py", ".cff", ".json"}:
            text = path.read_text(encoding="utf-8", errors="replace")
            if ABSOLUTE_USER_PATH.search(text):
                errors.append(f"absolute local user path: {relative}")

    accession_path = ROOT / "metadata" / "geo_accessions.tsv"
    with accession_path.open(encoding="utf-8", newline="") as handle:
        accessions = {row["accession"] for row in csv.DictReader(handle, delimiter="\t")}
    if accessions != EXPECTED_GEO:
        errors.append(f"GEO accession mismatch: {sorted(accessions)}")

    main_tables = list((ROOT / "results" / "main_tables").glob("*.tsv"))
    supplementary = list((ROOT / "results" / "supplementary_tables").glob("*.tsv"))
    figures = list((ROOT / "results" / "figures").glob("*.pdf"))
    if len(main_tables) != 3:
        errors.append(f"expected 3 main tables, found {len(main_tables)}")
    if len(supplementary) != 8:
        errors.append(f"expected 8 supplementary tables, found {len(supplementary)}")
    if len(figures) != 5:
        errors.append(f"expected 5 figure PDFs, found {len(figures)}")

    frozen = ROOT / "derived_data" / "frozen_candidate_sets" / "bfff_86_gene_set_frozen.tsv"
    bridge = ROOT / "results" / "main_tables" / "Table_3_16_bridge_candidates.tsv"
    stouffer = ROOT / "results" / "supporting_results" / "revised_main_gse77087" / "target_set" / "weighted_stouffer.tsv"
    with frozen.open(encoding="utf-8") as handle:
        frozen_rows = sum(1 for _ in handle) - 1
    with bridge.open(encoding="utf-8") as handle:
        bridge_rows = sum(1 for _ in handle) - 1
    if frozen_rows != 86:
        errors.append("frozen candidate set does not contain 86 rows")
    if bridge_rows != 16:
        errors.append("bridge-candidate table does not contain 16 rows")
    with stouffer.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    if round(float(row["weighted_stouffer_two_sided_p"]), 4) != 0.0675:
        errors.append("weighted-Stouffer value does not round to 0.0675")

    required_boundary = (
        "No analyzed cohort included exposure to Bufei Fanggan Formula",
        "do not establish BFFF efficacy",
    )
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for phrase in required_boundary:
        if phrase not in readme:
            errors.append(f"missing evidence-boundary phrase: {phrase}")

    manifest_path = ROOT / "MANIFEST.sha256"
    manifest_rows = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split(maxsplit=1)
        manifest_rows[relative.lstrip("*")] = digest
    expected_manifest_files = {
        "./" + path.relative_to(ROOT).as_posix()
        for path in files
        if path.name not in {"MANIFEST.sha256", "VALIDATION_REPORT.txt"}
    }
    if set(manifest_rows) != expected_manifest_files:
        errors.append("manifest file inventory mismatch")
    for relative, expected in manifest_rows.items():
        if sha256_file(ROOT / relative) != expected:
            errors.append(f"manifest checksum mismatch: {relative}")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    for phrase in ("given-names: Fan", "given-names: Yuyue", "given-names: Xuehui", "contact: true"):
        if phrase not in citation:
            errors.append(f"missing confirmed citation metadata: {phrase}")
    repository_url = "https://github.com/echoqf-lang/bfff-pediatric-viral-transcriptomics"
    if repository_url not in citation or repository_url not in readme:
        errors.append("confirmed repository URL is missing from citation or README metadata")

    required_license_files = (ROOT / "LICENSE", ROOT / "LICENSES" / "MIT.txt", ROOT / "LICENSES" / "CC-BY-4.0.txt")
    if not all(path.is_file() for path in required_license_files):
        errors.append("dual-license files are incomplete")
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    if zenodo.get("license") != "mit" or [x.get("name") for x in zenodo.get("creators", [])] != ["Qiu, Fan", "Yin, Yuyue", "Yang, Xuehui"]:
        errors.append("Zenodo metadata license or creator order mismatch")

    status = "PASS" if not errors else "FAIL"
    print(f"BFFF_DOI_REPOSITORY_{status}")
    print(f"files={len(files)} main_tables={len(main_tables)} supplementary_tables={len(supplementary)} figures={len(figures)}")
    print("release_ready=true author_metadata_required=false deposit_contact_confirmed=true license_files_present=true license_selection_confirmed=true institutional_release_confirmation_required=false repository_url_required=false external_upload_performed=true doi_reserved_or_published=false")
    for error in errors:
        print(f"ERROR: {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
