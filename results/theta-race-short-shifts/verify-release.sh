#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")" && pwd)
cd "$repo"
if [[ $# -ne 0 && $# -ne 4 ]]; then
  echo "Usage: bash verify-release.sh [mathlib4 upstream-src build-cache theta-prime-trace]" >&2
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
python3 verify.py >"$tmp/receipt"
python3 independent.py >>"$tmp/receipt"
python3 real_axis.py >>"$tmp/receipt"
mathrev=0df444a360eaa60ab8c11dca51a86af692955474
if [[ $# -eq 4 ]]; then
  math=$(cd "$1" && pwd)
  src=$(cd "$2" && pwd)
  mkdir -p "$3"
  build=$(cd "$3" && pwd)
  trace=$(cd "$4" && pwd)
else
  if [[ "$(uname -s)/$(uname -m)" != Linux/x86_64 ]]; then
    echo "Use the four-path replay with the pinned Lean and Mathlib on this platform." >&2
    exit 2
  fi
  fetch() { curl -fsSL --retry 4 "$1" -o "$2"; }
  lean_url=https://github.com/leanprover/lean4/releases/download/v4.33.1/lean-4.33.1-linux.tar.zst
  fetch "$lean_url" "$tmp/lean.tar.zst"
  echo "890afd185370f85666025b883914ab4f4b339136f8c96167b69cfb62aecaf235  $tmp/lean.tar.zst" | "${sum[@]}" -c -
  tar --zstd -xf "$tmp/lean.tar.zst" -C "$tmp"
  export PATH="$tmp/lean-4.33.1-linux/bin:$PATH"
  math="$tmp/mathlib4"
  git init -q "$math"
  git -C "$math" fetch -q --depth 1 https://github.com/leanprover-community/mathlib4.git "$mathrev"
  git -C "$math" checkout -q --detach FETCH_HEAD
  test "$(git -C "$math" rev-parse HEAD)" = "$mathrev"
  (cd "$math" && lake exe cache get)
  fetch https://github.com/CocoSgt/riemann-hypothesis-lean/releases/download/v2026-11-20/lean.tar.gz "$tmp/source.tar.gz"
  echo "9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc  $tmp/source.tar.gz" | "${sum[@]}" -c -
  tar -xzf "$tmp/source.tar.gz" -C "$tmp"
  src="$tmp/lean/src"
  build="$tmp/build"
  trace="$tmp/trace"
  mkdir -p "$trace"
  tracebase=https://raw.githubusercontent.com/leecrossley/theta-prime-trace/7a7909d2dc7bb9d8d0ec7b2dc8c78e70f93782d7
  for file in ThetaPrimeTrace.lean DyadicArithmetic.lean; do
    fetch "$tracebase/$file" "$trace/$file"
  done
fi
test "$(git -C "$math" rev-parse HEAD)" = "$mathrev"

(cd "$math" && lake env python3 "$repo/check_lean.py" "$src" "$build" "$trace" --jobs "${LEAN_JOBS:-4}") >>"$tmp/receipt"
diff -u evidence/check.txt "$tmp/receipt"
cat "$tmp/receipt"
