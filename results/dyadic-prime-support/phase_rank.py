"""Certify the rank on every phase cell by rational minors and their Bernstein coefficients."""

from fractions import Fraction as F
from itertools import pairwise
from math import comb
from support_recheck import BAD, characters, prime_data, small_primes

LIMIT = F(271,8)
GROUPS = ((0,), (1,), (0,1))


def prime_tables():
    residues = [characters(h) for h in range(2)]
    return [{p: prime_data(p,h,residues[h]) for p in small_primes(271) if p not in BAD[h]}
            for h in range(2)]


def boundaries(fields):
    marks = {F(1), F(2), LIMIT/32}
    for field in fields:
        for norm,_ in field.values():
            if norm <= 8*LIMIT:
                marks.update(F(norm,2**j) for j in range(9) if 1 < F(norm,2**j) < 2)
    return sorted(marks)


def coefficient_rows(fields, hosts, phase):
    last = 5 if 32*phase <= LIMIT else 4
    values = []
    for h in hosts:
        for j in range(1,last+2):
            z = 1 << j
            row = {}
            for p,(norm,sign) in fields[h].items():
                if norm > 8*LIMIT or not z*phase < norm < 4*z*phase:
                    continue
                row[p] = ((2*sign,-F(2*z*sign,norm)) if norm < 2*z*phase
                          else (-sign,F(4*z*sign,norm)))
            if row:
                values.append(row)
    primes = sorted(set().union(*values))
    return primes,[[row.get(p,(0,0)) for p in primes] for row in values]


def at(rows, phase):
    return [[a+b*phase for a,b in row] for row in rows]


def independent_minor(rows):
    """Return original row/column indices of a nonsingular minor and its determinant."""
    m = [list(map(F,row)) for row in rows]
    remaining = list(range(len(m)))
    chosen, columns = [], []
    det = F(1)
    for j in range(len(m[0]) if m else 0):
        k = len(chosen)
        i = next((i for i in range(k,len(m)) if m[i][j]),None)
        if i is None:
            continue
        m[i],m[k] = m[k],m[i]
        remaining[i],remaining[k] = remaining[k],remaining[i]
        chosen.append(remaining[k])
        columns.append(j)
        lead = m[k][j]
        det *= lead
        for i in range(k+1,len(m)):
            if m[i][j]:
                c = m[i][j]/lead
                m[i] = [x-c*y for x,y in zip(m[i],m[k])]
    return chosen,columns,det


def determinant(m):
    chosen, cols, product = independent_minor(m)
    if len(cols) != len(m):
        return F(0)
    inversions = sum(i>j for k,i in enumerate(chosen) for j in chosen[k+1:])
    return (-1)**inversions * product


def bernstein(values):
    """Convert values at i/d to Bernstein coefficients on the unit interval."""
    d = len(values)-1
    if not d:
        return list(values)
    power = [F(0)]*(d+1)
    basis = [F(1)]
    diffs = list(values)
    for j in range(d+1):
        for k,b in enumerate(basis):
            power[k] += b*diffs[0]
        if j<d:
            diffs = [(y-x)*d/(j+1) for x,y in pairwise(diffs)]
            root = F(j,d)
            basis = [x-y*root for x,y in zip([F(0)]+basis,basis+[F(0)])]
    return [sum(power[j]*F(comb(k,j),comb(d,j)) for j in range(k+1)) for k in range(d+1)]


def cell_rank(rows,lo,hi):
    # Rank above r is detected by a minor of degree <= r+1.
    chosen,cols,_ = independent_minor(at(rows,(lo+hi)/2))
    r = len(cols)
    for i in range(r+2):
        sample = lo+(hi-lo)*F(i+1,r+3)
        assert len(independent_minor(at(rows,sample))[1]) == r
    if not r:
        return r
    for index in (F(1,2),F(1,3),F(2,3),F(1,5),F(4,5)):
        chosen,cols,_ = independent_minor(at(rows,lo+(hi-lo)*index))
        assert len(cols) == r
        minor = [[rows[i][j] for j in cols] for i in chosen]
        values = [determinant(at(minor,lo+(hi-lo)*F(j,r))) for j in range(r+1)]
        coeff = bernstein(values)
        if any(coeff) and (min(coeff)>=0 or max(coeff)<=0):
            return r
    raise ArithmeticError(f"No sign-definite pivot minor on ({lo}, {hi})")


def spectrum(primes, rows, get_rank):
    target = get_rank(rows)
    return tuple(p for j,p in enumerate(primes) if get_rank([v[:j]+v[j+1:] for v in rows]) < target)


def expected_rank(hosts,t):
    n=5 if 32*t<=LIMIT else 4
    zeros=(2-int(t>F(31,16)),3-int(t>F(49,32)))
    return sum(n+1-zeros[h] for h in hosts)


def predicted(hosts,t):
    if hosts==(0,):
        return (31,) if t<=F(47,32) or t>F(31,16) else ()
    if hosts==(1,):
        return (7,) if t<=F(79,64) or t>F(49,32) else ()
    if t<=F(79,64) or t>F(31,16):
        return (7,31,47)
    if t<=F(97,64):
        return (31,47)
    if t<=F(49,32):
        return ()
    return (7,)


def shifted_example(fields):
    p,rows=coefficient_rows(fields,(0,1),F(3,2))
    all_rows=at(rows,F(3,2))
    u=all_rows[0]
    v=[a-b for a,b in zip(all_rows[1],all_rows[3])]
    assert [F(31,768)*(46*x-y) for x,y in zip(u,v)] == [int(q==31) for q in p]
    assert [F(47,768)*(14*x-17*y) for x,y in zip(u,v)] == [int(q==47) for q in p]
    print("At phase 3/2 the two pooled level combinations isolate 31 and 47; neither field alone isolates a prime.")


def main():
    fields=prime_tables()
    marks=boundaries(fields)
    assert len(marks)==14
    tests=0
    for hosts in GROUPS:
        for lo,hi in pairwise(marks):
            p,rows=coefficient_rows(fields,hosts,lo)
            assert len(independent_minor(at(rows,lo))[1]) == expected_rank(hosts,lo)
            assert spectrum(p,rows,lambda m:len(independent_minor(at(m,lo))[1])) == predicted(hosts,lo)
            mid=(lo+hi)/2
            p,rows=coefficient_rows(fields,hosts,mid)
            assert cell_rank(rows,lo,hi) == expected_rank(hosts,mid)
            assert spectrum(p,rows,lambda m:cell_rank(m,lo,hi)) == predicted(hosts,mid)
            tests+=len(p)+1
    print(f"Certified {tests} polynomial matrix ranks throughout 13 open phase cells, plus all 13 included endpoints on each of three field choices.")
    print("The pooled lists in phase order are (7,31,47), (31,47), (), (7), (7,31,47).")
    print("All available nonzero readings are independent; the only changes in the initial zero count occur strictly above phases 49/32 and 31/16.")
    shifted_example(fields)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
