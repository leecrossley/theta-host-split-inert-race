#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
if [[ $# -gt 1 || ( $# -eq 1 && $1 != --rebuild ) ]]; then
  echo "Usage: bash verify-dyadic-starts.sh [--rebuild]" >&2
  exit 2
fi
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
if [[ ${1:-} == --rebuild ]]; then
  c++ -std=c++17 -O3 word_sieve.cpp -o "$work/primary"
  c++ -std=c++17 -O3 word_sieve_recheck.cpp -o "$work/independent"
  python3 dyadic_records.py --rebuild primary --sieve "$work/primary"
  python3 dyadic_records.py --rebuild independent --sieve "$work/independent"
fi
python3 dyadic_records.py > "$work/records"
python3 dyadic_records_recheck.py >> "$work/records"
diff -u evidence/dyadic-starts.txt "$work/records"
cat "$work/records"
