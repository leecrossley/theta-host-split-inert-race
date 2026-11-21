#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")" && pwd)
if [[ $# -gt 1 ]]; then
  echo "Usage: bash verify-release.sh [pinned-mathlib-directory]" >&2
  exit 2
fi
if [[ $# -ge 1 ]]; then
  math=$(cd "$1" && pwd)
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
python3 real_axis.py >"$tmp/receipt"
python3 zeros.py >>"$tmp/receipt"
python3 independent.py >>"$tmp/receipt"
python3 amplitudes.py >>"$tmp/receipt"
python3 finite_examples.py >>"$tmp/receipt"
python3 deadline_spectrum.py >>"$tmp/receipt"
python3 early_common.py >>"$tmp/receipt"
python3 all_height_common.py >>"$tmp/receipt"
python3 fixed_windows.py >>"$tmp/receipt"
for method in direct split; do
  python3 finite_common.py "$method" >>"$tmp/receipt"
done

python3 three_starts.py > "$tmp/three"
python3 three_recheck.py >> "$tmp/three"
diff -u evidence/three-starts.txt "$tmp/three"
cat "$tmp/three" >> "$tmp/receipt"
python3 three_inverse.py > "$tmp/inverse"
python3 three_inverse_recheck.py >> "$tmp/inverse"
diff -u evidence/three-inverse.txt "$tmp/inverse"
cat "$tmp/inverse" >> "$tmp/receipt"

python3 inverse_monodromy.py > "$tmp/branch"
python3 inverse_branch_recheck.py >> "$tmp/branch"
python3 inverse_resolvent.py >> "$tmp/branch"
python3 inverse_resolvent_recheck.py >> "$tmp/branch"
diff -u evidence/inverse-monodromy.txt "$tmp/branch"
cat "$tmp/branch" >> "$tmp/receipt"

python3 target_box_corners.py > "$tmp/corners.tsv"
diff -u evidence/corner-enclosures.tsv "$tmp/corners.tsv"
python3 target_box_corners_recheck.py > "$tmp/corner-check"
diff -u evidence/corner-check.txt "$tmp/corner-check"
cat "$tmp/corner-check" >> "$tmp/receipt"

for trace in physical limiting; do
  python3 prescribed_anchor.py "$trace" >> "$tmp/prescribed"
  python3 prescribed_anchor_recheck.py "$trace" >> "$tmp/prescribed"
done
diff -u evidence/prescribed-anchor.txt "$tmp/prescribed"
cat "$tmp/prescribed" >> "$tmp/receipt"

python3 fixed_pair_balance.py > "$tmp/fixed-pair"
python3 fixed_pair_recheck.py >> "$tmp/fixed-pair"
diff -u evidence/fixed-pair-balance.txt "$tmp/fixed-pair"
cat "$tmp/fixed-pair" >> "$tmp/receipt"

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
test "$(git -C "$math" rev-parse HEAD)" = "$mathrev"
git -C "$math" diff --quiet HEAD -- Mathlib Mathlib.lean lakefile.lean lake-manifest.json lean-toolchain
(cd "$math" && "${sum[@]}" -c - <<'MATHLIB_HASHES'
3bba06a6dcfac62c11b0b06ffe6bbc953a093a92df1a58473b0daa84bd7e07f8  Mathlib.lean
1009a101f52a9d5d8ca57666f0256bae91abdb693e0f27b5a4622cffcc38cf03  lake-manifest.json
MATHLIB_HASHES
)
(cd "$math" && lake env python3 "$repo/check_lean.py") >>"$tmp/receipt"
diff -u evidence/check.txt "$tmp/receipt"
cat "$tmp/receipt"
