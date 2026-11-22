"""Certify first powers of two and all strict-margin records through exponent 34."""

from fractions import Fraction
from pathlib import Path
import argparse
import subprocess

from flint import arb, ctx

from word_cover import Race, LIMIT, real

ROOT = Path(__file__).resolve().parent
EVIDENCE = ROOT / "evidence"
FIELDS = ("547", "571")
SCALE = 10**7
LAST = 34
EDGES = tuple(1 << j for j in range(LAST + 5))


def read_bins(method, label):
    source = EVIDENCE / f"dyadic-{method}-{label}.txt"
    rows = [tuple(map(int, line.split())) for line in source.read_text().splitlines()]
    assert tuple(row[0] for row in rows) == EDGES
    assert all(len(row) == 5 and row[1] <= row[2] and row[3] <= row[4] for row in rows)
    return rows


def rebuild(binary, method):
    divisions, bits = (2**15, 58) if method == "primary" else (2**12, 62)
    nodes = []
    for j in range(divisions + 1):
        lower = ((1 + arb(j) / divisions).log() * 2**bits).floor().unique_fmpz()
        assert lower is not None
        nodes.append(int(lower))
    data = "".join(f"{value}\n" for value in (*nodes, *EDGES))
    for label in FIELDS:
        run = subprocess.run([str(binary), label], text=True, input=data, capture_output=True, check=True)
        expected = (EVIDENCE / f"dyadic-{method}-{label}.txt").read_text()
        if run.stdout != expected:
            raise ArithmeticError(f"{method} sieve differs from the archived bins for {label}")
        print(f"{method} sieve reproduced all {len(EDGES)} integer bins for {label}.")


def primary_moves():
    result = {}
    for field in FIELDS:
        bins = {}
        totals = [0] * 4
        for edge, *values in read_bins("primary", field):
            totals = [a + b for a, b in zip(totals, values, strict=True)]
            bins[edge] = tuple(totals)
        race = Race(bins, 58, 96)
        for n in range(LAST + 2):
            result[field, n] = race.move(Fraction(1 << n), 1)
    return result


def record_table(moves):
    records = []
    last = {}
    first = {}
    for n in range(LAST + 1):
        keys = [(field, n + j) for j in (0, 1) for field in FIELDS]
        if any(moves[key].is_zero() for key in keys):
            assert n in (0, 1, 2)
            continue
        assert all(moves[key] > 0 or moves[key] < 0 for key in keys)
        word = "".join("+" if moves[key] > 0 else "-" for key in keys)
        weakest = [key for key in keys if all(key == other or abs(moves[key]) < abs(moves[other]) for other in keys)]
        assert len(weakest) == 1
        key = weakest[0]
        value = abs(moves[key])
        if value > real(LIMIT):
            first.setdefault(word, n)
        else:
            assert value < real(LIMIT)
        if word in last:
            previous = last[word]
            if key == previous or value < abs(moves[previous]):
                continue
            assert value > abs(moves[previous])
        last[word] = key
        lower = (SCALE * value).floor().unique_fmpz()
        assert lower is not None
        records.append((word, n, *key, int(lower)))
    return sorted(records), first


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", choices=("primary", "independent"))
    parser.add_argument("--sieve", type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 200
    if args.rebuild:
        if args.sieve is None:
            parser.error("--rebuild requires --sieve")
        rebuild(args.sieve.resolve(), args.rebuild)
        return
    moves = primary_moves()
    rows, first = record_table(moves)
    alpha, beta = moves["571", 12], -moves["571", 35]
    assert arb("3.41437510") < alpha < arb("3.41437512")
    assert arb("6.75643706") < beta < arb("6.75643708")
    assert first == {"+--+": 3, "++--": 6, "--++": 7, "++-+": 8, "---+": 11,
                     "-+-+": 12, "+-++": 15, "+++-": 16, "-+--": 18, "-++-": 21,
                     "+---": 22, "-+++": 24, "----": 26, "+-+-": 34}
    final = {row[0]: (row[2], row[3]) for row in rows}
    assert final["-+-+"] == ("571", 12)
    assert all(abs(moves[key]) > alpha for word, key in final.items() if word != "-+-+")
    lines = ["word\tstart_exponent\tfield\tlimiting_exponent\tfloor_margin_times_10000000"]
    lines.extend("\t".join(map(str, row)) for row in rows)
    assert "\n".join(lines) + "\n" == (EVIDENCE / "dyadic-records.tsv").read_text()
    assert len(rows) == 23 and len({(row[2], row[3]) for row in rows}) == 16
    assert len(first) == 14 and first["+-+-"] == 34 and first["-+-+"] == 12
    assert "++++" not in first and "--+-" not in first
    print("35 exact powers of two, 14 strict words, 23 margin records and 16 distinct limiting readings.")
    print("First opposed runs at margin 1.9735: 547 gains at exponent 34; 571 gains at exponent 12.")
    print("The four-gain and --+- patterns do not occur through exponent 34, even at strict margin zero.")


if __name__ == "__main__":
    main()
