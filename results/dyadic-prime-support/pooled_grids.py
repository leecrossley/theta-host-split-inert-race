"""Check every pair of phase cells for the joint two-grid prime classification."""
from csv import DictReader
from fractions import Fraction as F
from itertools import combinations_with_replacement, pairwise, product
from pathlib import Path

from minimum_grids import field_primes
from phase_rank import determinant, independent_minor
from support_recheck import hat

HERE = Path(__file__).resolve().parent
LIMIT = 431
OMITTED_TARGETS = (17, 337, 353, 367)


def reduced(rows):
    """A rational row basis in reduced form, keyed by its leading column."""
    basis = {}
    for source in rows:
        row = list(map(F,source))
        for j in sorted(basis):
            if row[j]:
                c = row[j]
                row = [x-c*y for x,y in zip(row,basis[j])]
        pivot = next((j for j,x in enumerate(row) if x),None)
        if pivot is None:
            continue
        row = [x/row[pivot] for x in row]
        for j,earlier in basis.items():
            if earlier[pivot]:
                basis[j] = [x-earlier[pivot]*y for x,y in zip(earlier,row)]
        basis[pivot] = row
    return basis


def grid_cells(fields):
    masses = {m for field in fields for m,p,s in field}
    marks = {F(1),F(2)}
    for m in masses | {LIMIT}:
        marks.add(F(m,2**(m.bit_length()-1)))
    bounds = sorted(marks)
    return [(v,v) for v in bounds[:-1]] + list(pairwise(bounds))


def entries(field, primes, phase, j):
    z = phase*2**j
    return [256*m*s*hat(m,z) if p in field else 0
            for p in primes for m,s in [field.get(p,(0,0))]]


def levels(cell, fields, primes):
    left,right = cell
    result = []
    for j in range(9):
        if 4*right*2**j > LIMIT:
            break
        for field in fields:
            a,b = entries(field,primes,left,j),entries(field,primes,right,j)
            if any(a) or any(b):
                result.append((a,b))
    return result


def polynomial_minor(blocks, selected, columns, cells):
    # Homogeneous variable t>0 turns each open-cell row into L+tR.
    # Expand by multilinearity and collect equal monomials before testing signs.
    variable = [i for i,row in enumerate(selected) if cells[blocks[row][0]][0]!=cells[blocks[row][0]][1]]
    coefficients = {}
    for bits in product((0,1),repeat=len(variable)):
        choose = dict(zip(variable,bits))
        powers = [0,0]
        matrix = []
        for i,row in enumerate(selected):
            side,ends = blocks[row]
            bit = choose.get(i,0)
            powers[side] += bit
            matrix.append([ends[bit][j] for j in columns])
        terms = tuple(powers)
        coefficients[terms] = coefficients.get(terms,0)+determinant(matrix)
    values = tuple(coefficients.values())
    assert any(values) and (min(values)>=0 or max(values)<=0)


def dependence(blocks, selected, cols, target, cells):
    others = [j for j in cols if j!=target]
    r = len(others)
    assert target in cols and len(set(cols))==len(cols) and len(set(selected))==len(selected)==r
    polynomial_minor(blocks,selected,others,cells)
    relevant = [(side,ends) for side,ends in blocks if any(any(end[c] for c in cols) for end in ends)]
    if len(relevant)<=r:
        return
    degree = [0 if a==b else min(r+1,sum(side==s for side,_ in relevant)) for s,(a,b) in enumerate(cells)]
    for pair in product(*(range(1,n+2) for n in degree)):
        rows = [[l[j]+pair[side]*v[j] for j in cols] for side,(l,v) in relevant]
        assert len(independent_minor(rows)[1])<=r


def verify():
    fields = [field_primes(h,LIMIT-1) for h in range(2)]
    primes = sorted({p for field in fields for m,p,s in field})
    assert len(primes)==16
    # There is a single common norm for each prime wherever it occurs in this prefix.
    assert all(len({m for field in fields for m,q,s in field if q==p})==1 for p in primes)
    cells = grid_cells(fields)
    assert len(cells)==36
    maps = [{p:(m,s) for m,p,s in field} for field in fields]
    rows = [levels(c,maps,primes) for c in cells]
    endpoint = [list(reduced([v for ends in group for v in ends]).values()) for group in rows]
    certificates = {}
    with (HERE/"evidence/pooled-grid-minors.tsv").open() as stream:
        for item in DictReader(stream,delimiter="\t"):
            key = tuple(int(item[k]) for k in ("cell_a","cell_b","target"))
            assert key not in certificates
            certificates[key] = ([primes.index(int(x)) for x in item["dependent_primes"].split(",")],
                                 [int(x) for x in item["minor_rows"].split(",")])
    direct = checked = 0
    for i,k in combinations_with_replacement(range(len(cells)),2):
        enlarged = reduced(endpoint[i]+endpoint[k])
        active_rows = [(s,ends) for s,ix in enumerate((i,k)) for ends in rows[ix]]
        for p in OMITTED_TARGETS:
            j = primes.index(p)
            if j not in enlarged or sum(bool(x) for x in enlarged[j])!=1:
                direct += 1
                continue
            cols,selected = certificates.pop((i,k,p))
            dependence(active_rows,selected,cols,j,(cells[i],cells[k]))
            checked += 1
    assert not certificates and (direct,checked)==(2407,257)
    print("All 2,664 joint exclusions checked across 666 unordered phase-cell pairs: 2,407 endpoint-span tests and 257 polynomial minor certificates.")


def check_recovery():
    fields = [field_primes(h,430) for h in range(2)]
    terms = []
    with (HERE/"evidence/pooled-grid-formula.tsv").open() as stream:
        for row in DictReader(stream,delimiter="\t"):
            assert int(row["prime"])==97
            terms.append((int(row["host"]),F(row["phase"]),int(row["level_exponent"]),F(row["coefficient"])))
    assert {a for h,a,j,c in terms}=={F(3,2),F(7,4)}
    assert max(4*a*2**j for h,a,j,c in terms)<=430
    coeff = {}
    for h,a,j,c in terms:
        for m,p,s in fields[h]:
            coeff[p] = coeff.get(p,0)+c*s*hat(m,a*2**j)
    assert {p:v for p,v in coeff.items() if v}=={97:1}
    print("Two grids recover log_2(97) only when both fields are available. Their five rational terms were checked at every contributing prime.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    verify()
    check_recovery()
