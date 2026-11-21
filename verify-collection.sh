#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -f "$root/shared-files.json" ]; then
    python3 "$root/tools/shared_files.py"
fi
for result in "$root" "$root"/results/*; do
    [ -f "$result/SHA256SUMS" ] || continue
    printf '\n%s\n' "$result"
    if command -v sha256sum >/dev/null 2>&1; then
        (cd "$result" && sha256sum -c SHA256SUMS)
    else
        (cd "$result" && shasum -a 256 -c SHA256SUMS)
    fi
done
