"""Reconstruct three octave weights with rational intervals and verify every first-passage cell."""

from csv import DictReader
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path

EVIDENCE = Path(__file__).resolve().parent / "evidence"
FIELDS = {"547": (547, 727), "571": (571, 599)}
SCALE = 10**7


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __add__(self, other):
        other = as_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __mul__(self, other):
        other = as_interval(other)
        products = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        assert self.lo > 0
        return Interval(1 / self.hi, 1 / self.lo)


def as_interval(x):
    if isinstance(x, Interval):
        return x
    return Interval(F(x), F(x))


def constants():
    scale = 10**85
    root = isqrt(2 * scale**2)
    r = Interval(F(root, scale), F(root + 1, scale))
    terms = 90
    lower = sum(F(2, (2*j + 1) * 3**(2*j + 1)) for j in range(terms))
    rest = F(9, 4 * (2*terms + 1) * 3**(2*terms + 1))
    log2 = Interval(lower, lower + rest)
    B = 2 * r * (r + (-1)) * (r + (-1))
    return r, 8 * (B * log2).reciprocal()


def load_octaves(label):
    def load(method):
        path = EVIDENCE / f"dyadic-{method}-{label}.txt"
        return [tuple(map(int, line.split())) for line in path.read_text().splitlines()]
    first, second = load("primary"), load("independent")
    assert len(first) == len(second) == 39
    bins = []
    for j, (u, v) in enumerate(zip(first, second, strict=True)):
        assert len(u) == len(v) == 5 and u[0] == v[0] == 2**j
        for pos, shift in ((1, 4), (3, 6)):
            assert max(u[pos] << shift, v[pos]) <= min(u[pos + 1] << shift, v[pos + 1])
        bins.append((Interval(F(v[1], 2**62), F(v[2], 2**62)),
                     Interval(F(v[3], 2**102), F(v[4], 2**102))))
    return bins


def moves_by_signed_weight():
    root, common = constants()
    result = {}
    for label in FIELDS:
        bins = load_octaves(label)
        for k in range(36):
            y = 2**k
            (a, b), (c, d), (e, f) = bins[k + 1:k + 4]
            first = -2 * a + (2*y) * b
            middle = (1 + root) * c + (-4 + (-2)*root) * y * d
            last = (-root.reciprocal()) * e + 4 * root * y * f
            root_y = as_interval(2**(k//2)) * (root if k%2 else 1)
            result[label, k] = common * root_y.reciprocal() * (first + middle + last)
    return result


def read_tsv(name):
    with (EVIDENCE / name).open() as f:
        return list(DictReader(f, delimiter="\t"))


def check_pivots(keys):
    pivots = read_tsv("dyadic-pivots.tsv")
    assert [(int(row["exponent"]), row["field"]) for row in pivots] == sorted(keys)
    last = -1
    for row in pivots:
        j, field = int(row["exponent"]), row["field"]
        p, f = int(row["prime"]), int(row["residue_degree"])
        u, v = FIELDS[field]
        assert p >= 3 and p % 2 and all(p % t for t in range(3, isqrt(p) + 1, 2))
        omitted = (2, u, v, *((343057911569, 343070636977) if field == "547"
                                else (413718278401, 603339156001)))
        assert p not in omitted and 1 << (j + 2) < p**f < 1 << (j + 3)
        su = {t*t % u for t in range(1, u)}
        sv = {t*t % v for t in range(1, v)}
        degree = next(d for d in (1, 2, 4) if pow(p, d, 16) in (1, 15)
                      and (d % 2 == 0 or (p % u in su) == (p % v in sv)))
        assert degree == f
        assert last == -1 or (j > last and p > 1 << (last + 3))
        last = j
    print("16 prime pivots verify the triangular coefficient minor, from 7 to 137438953567.")


def check_words(moves):
    table = read_tsv("dyadic-records.tsv")
    keys = {(int(row["limiting_exponent"]), row["field"]) for row in table}
    assert len(keys) == 16 and len(table) == 23
    zero = as_interval(0)
    for letters in product("-+", repeat=4):
        word = "".join(letters)
        active = None
        records = {int(row["start_exponent"]): row for row in table if row["word"] == word}
        for n in range(35):
            readings = [(field, n + offset) for offset in (0, 1) for field in FIELDS]
            signed = {key: moves[key] if sign == "+" else -moves[key]
                      for sign, key in zip(letters, readings, strict=True)}
            if n in records:
                row = records[n]
                key = row["field"], int(row["limiting_exponent"])
                limiting = signed[key]
                assert limiting.lo > (active[2].hi if active else 0)
                assert all(limiting.hi < signed[c].lo for c in readings if c != key)
                floor = int(row["floor_margin_times_10000000"])
                assert floor < limiting.lo * SCALE <= limiting.hi * SCALE < floor + 1
                active = key, letters[readings.index(key)], limiting
            else:
                old = active[2] if active else zero
                same_reading = (active is not None and active[0] in readings
                                and letters[readings.index(active[0])] == active[1])
                assert same_reading or any(value.hi <= old.lo for value in signed.values())
    for reading, lower, upper in ((moves["571", 12], "3.41437510", "3.41437512"),
                                  (-moves["571", 35], "6.75643706", "6.75643708")):
        assert F(lower) < reading.lo <= reading.hi < F(upper)
    check_pivots(keys)
    print("All 560 signed start/word tests passed with rational endpoints and independently weighted octaves.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_words(moves_by_signed_weight())
