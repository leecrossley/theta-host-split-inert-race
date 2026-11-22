#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -f "$root/shared-files.json" ]; then
    python3 "$root/tools/shared_files.py"
fi
python3 - "$root" <<'PYTHON'
from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
index = (root / "RESULTS.md").read_text()
for name in sorted(set(re.findall(r"\]\((results/[^/#?)\s]+)", index))):
    if not (root / name / "SHA256SUMS").is_file():
        sys.exit(f"Missing result manifest: {name}/SHA256SUMS")
PYTHON

if command -v sha256sum >/dev/null 2>&1; then
    checksum=sha256sum
else
    checksum="shasum -a 256"
fi
(cd "$root" && $checksum -c COLLECTION.sha256)
if [ -f "$root/shared-sources.sha256" ]; then
    (cd "$root" && $checksum -c shared-sources.sha256)
fi
if [ -e "$root/.git" ]; then
    python3 "$root/tools/collection_manifest.py"
fi

for result in "$root" "$root"/results/*; do
    [ -d "$result" ] || continue
    if [ ! -f "$result/SHA256SUMS" ]; then
        echo "Missing result manifest: $result/SHA256SUMS" >&2
        exit 1
    fi
    printf '\n%s\n' "$result"
    (cd "$result" && $checksum -c SHA256SUMS)
done
