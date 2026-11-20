#!/usr/bin/env bash
set -euo pipefail

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'Use the pinned Mathlib instructions in README.md on this platform.' >&2
  exit 2
fi

repo=$(cd "$(dirname "$0")" && pwd)
(cd "$repo" && sha256sum --check SHA256SUMS)
scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
cd "$scratch"

download() {
  curl --fail --location --silent --show-error --retry 4 --output "$2" "$1"
  printf '%s  %s\n' "$3" "$2" | sha256sum --check -
}
download https://github.com/leanprover/lean4/releases/download/v4.33.1/lean-4.33.1-linux.tar.zst \
  lean.tar.zst 890afd185370f85666025b883914ab4f4b339136f8c96167b69cfb62aecaf235
tar --zstd -xf lean.tar.zst
rm -- lean.tar.zst
export PATH="$scratch/lean-4.33.1-linux/bin:$PATH"

mathlib=0df444a360eaa60ab8c11dca51a86af692955474
git init -q mathlib4
git -C mathlib4 remote add origin https://github.com/leanprover-community/mathlib4.git
git -C mathlib4 fetch --quiet --depth 1 origin "$mathlib"
git -C mathlib4 checkout --quiet --detach FETCH_HEAD
test "$(git -C mathlib4 rev-parse HEAD)" = "$mathlib"
download https://github.com/CocoSgt/riemann-hypothesis-lean/releases/download/v2026-11-20/lean.tar.gz \
  source.tar.gz 9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc
tar -xzf source.tar.gz
rm -- source.tar.gz
source_dir="$scratch/lean/src"
test -d "$source_dir"


cd mathlib4
lake exe cache get
lake env python3 "$repo/verify.py" "$source_dir" "$scratch/build" "$@" | tee "$scratch/check.log"
grep -E '^(Checked [0-9]+ upstream|Compiler commit:|Kernel replay passed:|Axioms for )' \
  "$scratch/check.log" > "$scratch/check.txt"
diff -u "$repo/evidence/check.txt" "$scratch/check.txt"
