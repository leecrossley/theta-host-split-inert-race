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
python3 arithmetic_start_grid.py primary "$work/first" > "$work/arithmetic"
python3 arithmetic_start_grid.py independent "$work/second" >> "$work/arithmetic"
python3 recover_thresholds.py >> "$work/arithmetic"
python3 single_field_recovery.py >> "$work/arithmetic"
python3 recover_thresholds_recheck.py >> "$work/arithmetic"
diff -u evidence/arithmetic-starts.txt "$work/arithmetic"
cat "$work/arithmetic"
