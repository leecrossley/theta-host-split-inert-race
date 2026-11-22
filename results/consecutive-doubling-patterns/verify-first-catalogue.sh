#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
export PYTHONDONTWRITEBYTECODE=1
python3 -c 'import flint; assert flint.__version__ == "0.9.0"'
c++ -std=c++17 -O3 word_sieve.cpp -o "$work/first"
c++ -std=c++17 -O3 word_sieve_recheck.cpp -o "$work/second"
python3 word_cover.py primary "$work/first" > "$work/transcript"
python3 word_cover.py independent "$work/second" >> "$work/transcript"
python3 word_peak.py >> "$work/transcript"
diff -u evidence/first-catalogue.txt "$work/transcript"
cat "$work/transcript"
python3 clock_phases.py > "$work/phases"
python3 clock_phases_recheck.py >> "$work/phases"
diff -u evidence/clock-phases.txt "$work/phases"
cat "$work/phases"
