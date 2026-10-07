#!/bin/sh
set -eu
repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
manifest="$repo_root/MANIFEST.sha256"
temp_manifest="$repo_root/.MANIFEST.sha256.tmp"
cd "$repo_root"
find . -type f ! -path './.git/*' ! -path './MANIFEST.sha256' ! -path './VALIDATION_REPORT.txt' ! -path './.MANIFEST.sha256.tmp' -print \
  | LC_ALL=C sort \
  | while IFS= read -r file; do shasum -a 256 "$file"; done > "$temp_manifest"
mv "$temp_manifest" "$manifest"
