#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")" && pwd)
if [[ $# -gt 2 ]]; then
  echo "Usage: bash verify-release.sh [pinned-mathlib-directory [upstream-records-directory]]" >&2
  exit 2
fi
if [[ $# -ge 1 ]]; then
  math=$(cd "$1" && pwd)
fi
if [[ $# -eq 2 ]]; then
  refs=$(cd "$2" && pwd)
fi
cd "$repo"
if command -v sha256sum >/dev/null 2>&1; then
  sum=(sha256sum)
else
  sum=(shasum -a 256)
fi
"${sum[@]}" -c SHA256SUMS
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
export PYTHONDONTWRITEBYTECODE=1
python3 -c 'import flint; assert flint.__version__ == "0.9.0"'
if [[ $# -eq 2 ]]; then
  python3 check_source.py "$refs" >"$tmp/receipt"
else
  python3 check_source.py >"$tmp/receipt"
fi
python3 check_predecessor.py >>"$tmp/receipt"
python3 phase_budget.py >>"$tmp/receipt"
python3 check_budget.py >>"$tmp/receipt"

mathrev=0df444a360eaa60ab8c11dca51a86af692955474
if [[ $# -eq 0 ]]; then
  if [[ "$(uname -s)/$(uname -m)" != Linux/x86_64 ]]; then
    echo "Supply the pinned Mathlib directory and Lean toolchain on this platform." >&2
    exit 2
  fi
  curl -fsSL --retry 4 https://github.com/leanprover/lean4/releases/download/v4.33.1/lean-4.33.1-linux.tar.zst -o "$tmp/lean.tar.zst"
  echo "890afd185370f85666025b883914ab4f4b339136f8c96167b69cfb62aecaf235  $tmp/lean.tar.zst" | "${sum[@]}" -c -
  tar --zstd -xf "$tmp/lean.tar.zst" -C "$tmp"
  export PATH="$tmp/lean-4.33.1-linux/bin:$PATH"
  math="$tmp/mathlib4"
  git init -q "$math"
  git -C "$math" fetch -q --depth 1 https://github.com/leanprover-community/mathlib4.git "$mathrev"
  git -C "$math" checkout -q --detach FETCH_HEAD
  test "$(git -C "$math" rev-parse HEAD)" = "$mathrev"
  (cd "$math" && lake exe cache get)
fi
if [[ "$(git -C "$math" rev-parse HEAD 2>/dev/null)" != "$mathrev" ]]; then
  echo "The Mathlib directory must be a Git checkout at $mathrev." >&2
  exit 2
fi
if ! git -C "$math" diff --quiet HEAD --; then
  echo "The Mathlib checkout has local changes to tracked files." >&2
  exit 2
fi
(cd "$math" && "${sum[@]}" -c - <<'MATHLIB_HASHES'
3bba06a6dcfac62c11b0b06ffe6bbc953a093a92df1a58473b0daa84bd7e07f8  Mathlib.lean
86d3aaad83e64cb64a24167fd4ddb2eeec5ca4a19a9298459e87a57fec151040  Mathlib/NumberTheory/LSeries/DirichletContinuation.lean
1009a101f52a9d5d8ca57666f0256bae91abdb693e0f27b5a4622cffcc38cf03  lake-manifest.json
MATHLIB_HASHES
)
(cd "$math" && lake env python3 "$repo/check_lean.py") >>"$tmp/receipt"
diff -u evidence/check.txt "$tmp/receipt"
cat "$tmp/receipt"
