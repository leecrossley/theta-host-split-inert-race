#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")" && pwd)
cd "$repo"
if [[ $# -gt 1 ]]; then
  echo "Usage: bash verify-release.sh [pinned-mathlib-directory]" >&2
  exit 2
fi
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
for file in verify.py square_factors.py independent.py finite_windows.py finite_spectrum.py power_window.py carrier_window.py zero_filter.py separated_zeros.py tilted_moments.py; do
  python3 "$file" >>"$tmp/receipt"
done
for method in direct split; do
  python3 finite_deadlines.py "$method" >>"$tmp/receipt"
done

mathrev=0df444a360eaa60ab8c11dca51a86af692955474
if [[ $# -eq 1 ]]; then
  math=$(cd "$1" && pwd)
else
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
git -C "$math" diff --quiet HEAD --
(cd "$math" && lake env python3 "$repo/check_lean.py") >>"$tmp/receipt"
diff -u evidence/check.txt "$tmp/receipt"
cat "$tmp/receipt"
