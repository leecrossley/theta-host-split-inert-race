#!/usr/bin/env bash
set -euo pipefail

root=${RUNNER_TEMP:?}
tools="$root/lean-4.33.1-linux"
math="$root/mathlib4"
export PATH="$tools/bin:$root/short-venv/bin:$PATH"

revision() {
  case "$1" in
    16) echo 6444e730934116ee32554532b71f45e1273bc1b4 ;;
    43) echo 5046457a319c5ebce0db2c2894a64533b75953fb ;;
    45) echo 4a80fe0cc7b5360e713b54941f99fdeb4ac35e11 ;;
    *) return 2 ;;
  esac
}

if [[ "$1" == prepare ]]; then
  for number in 16 43 45; do
    sha=$(revision "$number")
    git fetch --quiet --no-tags --depth=1 origin "$sha"
    git worktree add --quiet --detach "$root/pr-$number" "$sha"
    test "$(git -C "$root/pr-$number" rev-parse HEAD)" = "$sha"
  done
  python3 -m venv "$root/short-venv"
  "$root/short-venv/bin/python" -m pip install -r "$root/pr-16/results/theta-race-short-shifts/requirements.txt"
  cd "$root"
  curl -fsSL --retry 4 -o lean.tar.zst https://github.com/leanprover/lean4/releases/download/v4.33.1/lean-4.33.1-linux.tar.zst
  echo "890afd185370f85666025b883914ab4f4b339136f8c96167b69cfb62aecaf235  lean.tar.zst" | sha256sum -c -
  tar --zstd -xf lean.tar.zst
  rm -- lean.tar.zst
  mathrev=0df444a360eaa60ab8c11dca51a86af692955474
  git init -q "$math"
  git -C "$math" fetch -q --depth 1 https://github.com/leanprover-community/mathlib4.git "$mathrev"
  git -C "$math" checkout -q --detach FETCH_HEAD
  test "$(git -C "$math" rev-parse HEAD)" = "$mathrev"
  (cd "$math" && lake exe cache get)
  curl -fsSL --retry 4 -o short-source.tar.gz https://github.com/CocoSgt/riemann-hypothesis-lean/releases/download/v2026-11-20/lean.tar.gz
  echo "9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc  short-source.tar.gz" | sha256sum -c -
  tar -xzf short-source.tar.gz
  rm -- short-source.tar.gz
  mkdir "$root/short-build"
  exit 0
fi

number=$1
sha=$(revision "$number")
checkout="$root/pr-$number"
test "$(git -C "$checkout" rev-parse HEAD)" = "$sha"
git -C "$checkout" diff --exit-code HEAD -- .
cd "$checkout/results/theta-race-short-shifts"
cmp requirements.txt "$root/pr-16/results/theta-race-short-shifts/requirements.txt"
export LEAN_NUM_THREADS=0 LEAN_JOBS=2
(
  while sleep 180 >/dev/null 2>&1; do
    count=$(find "$root/short-build" -name '*.olean' | wc -l)
    echo "PR $number replay has $count compiled modules in its build directory"
  done
) &
progress=$!
trap 'kill "$progress" 2>/dev/null || true' EXIT
bash verify-release.sh "$math" "$root/lean/src" "$root/short-build" "$checkout/results/theta-prime-trace"
printf 'PR #%s at `%s` passed its original release verifier.\n\n' "$number" "$sha" >> "$GITHUB_STEP_SUMMARY"
