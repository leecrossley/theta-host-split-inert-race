"""Determine whether each single prime coordinate belongs to the pooled early row spaces."""
from fractions import Fraction
from support_recheck import BAD, characters, hat, prime_data, rank, small_primes


def modular_rank(rows, modulus=65537):
    """Rank by column-wise Gauss-Jordan elimination in a finite field."""
    if not rows:
        return 0
    a = [[x.numerator * pow(x.denominator, -1, modulus) % modulus for x in row]
         for row in rows]
    pivot = 0
    for col in range(len(a[0])):
        row = next((j for j in range(pivot, len(a)) if a[j][col]), None)
        if row is None:
            continue
        a[pivot], a[row] = a[row], a[pivot]
        inv = pow(a[pivot][col], -1, modulus)
        leading = [t*inv % modulus for t in a[pivot]]
        for j in range(pivot+1, len(a)):
            factor = a[j][col]
            a[j] = [(t-factor*u) % modulus for t,u in zip(a[j], leading)]
        pivot += 1
    return pivot


def check():
    residues = [characters(h) for h in (0, 1)]
    fields = [{p: prime_data(p, h, residues[h]) for p in small_primes(1024) if p not in BAD[h]}
              for h in (0, 1)]

    def row(h, z):
        return {p: sign*hat(m,z) for p,(m,sign) in fields[h].items() if z < m < 4*z}

    assert row(0, 8) == {31: Fraction(1,31)}
    assert row(1, 16) == {7: -Fraction(15,49)}
    assert row(0, 16) == {7: -Fraction(15,49), 31: Fraction(30,31), 47: -Fraction(17,47)}
    print("Three exact rows recover log_2(31), log_2(7) and log_2(47).")
    for n in range(8):
        entries = [row(h,1 << j) for j in range(1,n+2) for h in (0,1)]
        columns = sorted(set().union(*entries))
        matrix = [[entry.get(p,0) for p in columns] for entry in entries]
        dimension = rank(matrix)
        assert modular_rank(matrix) == dimension
        recovered = []
        for p in columns:
            extended = matrix + [[int(q == p) for q in columns]]
            extended_rank = rank(extended)
            assert modular_rank(extended) == extended_rank
            if extended_rank == dimension:
                recovered.append(p)
        expected = [] if n < 2 else [31] if n == 2 else [7,31,47]
        assert recovered == expected
        assert len(columns) == (0,0,1,3,8,11,20,36)[n]
        assert dimension == (0,0,1,3,5,7,9,11)[n]
        print(f"Through exponent {n}: {len(columns)} prime columns, rank {dimension}, isolated {recovered}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check()
