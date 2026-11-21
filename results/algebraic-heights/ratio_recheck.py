"""Independent determinants from exact tent endpoints and GMP prime-certificate replay."""
import csv
from fractions import Fraction
from pathlib import Path

from flint import fmpz, nmod_mat
from height_ratios import SLOPE_COLUMNS, INTERCEPT_COLUMNS

HERE = Path(__file__).resolve().parent
HOSTS = {"A": (547,727), "B": (571,599)}
OMISSIONS = {2,547,727,571,599,343057911569,343070636977,413718278401,603339156001}


def duration_data(p, host):
    u,v = HOSTS[host]
    characters = [1 if pow(p,(q-1)//2,q)==1 else -1 for q in (u,v)]
    for f in (1,2,4):
        if pow(p,f,16) in (1,15) and (characters[0]*characters[1])**f == 1:
            return p**f, characters[0]**f
    raise ArithmeticError((host,p))


def reduced(value,modulus):
    return value.numerator*pow(value.denominator,-1,modulus) % modulus


def single_tent(y,n):
    return max(Fraction(),min(2*(1-y/n),4*y/n-1))


def matrix_entry(cell,p,modulus,h,intercept):
    field,left,right = cell
    n,sign = duration_data(p,field)
    assert all(duration_data(p,i)[0] <= 1_000_000 for i in HOSTS)
    l,r = Fraction(left,8),Fraction(right,8)
    # The values at the knots are continuous. Recover the unique affine
    # coefficients from these two endpoint observations.
    def observation(y):
        return sign*(h*reduced(single_tent(2*y,n),modulus) -
                     reduced(single_tent(y,n),modulus))
    value = observation(l)
    slope = (observation(r)-value) * pow(reduced(r-l,modulus),-1,modulus) % modulus
    return (value-slope*reduced(l,modulus)) % modulus if intercept else slope


def minors():
    with (HERE/"evidence/exact-crossings.tsv").open() as f:
        cells = list(dict.fromkeys((r["host"],int(r["left_eighths"]),int(r["right_eighths"]))
                                   for r in csv.DictReader(f,delimiter="\t")))
    residues = []
    for intercept,columns in ((False,SLOPE_COLUMNS),(True,INTERCEPT_COLUMNS)):
        for modulus in (1_000_000_007,2_147_483_647):
            h = pow(pow(2,(modulus+1)//4,modulus),-1,modulus)
            matrix = [[matrix_entry(c,p,modulus,h,intercept) for p in columns] for c in cells]
            residues.append(int(nmod_mat(matrix,modulus).det()))
    assert residues == [303857281,674797901,290245403,1234467739]
    print("Exact differences of the two original tents at every bounding knot and FLINT determinants "
          "give the same four protected minor residues.")


def check_proth(p,base):
    assert p not in OMISSIONS and p % 16 == 1
    if base == 0:
        assert fmpz(p).is_prime()
    else:
        two_part = (p-1) & (1-p)
        odd_part = (p-1)//two_part
        assert odd_part % 2 and 0 < odd_part < two_part
        assert 1 < base < p and pow(fmpz(base),(p-1)//2,fmpz(p)) == p-1
    traces = [duration_data(p,host) for host in HOSTS]
    assert all(n==p for n,sign in traces)
    return int(traces[0][1] == traces[1][1])


def chains():
    with (HERE/"evidence/ratio-primes.tsv").open() as f:
        small = list(csv.DictReader(f,delimiter="\t"))
    for row in small:
        e = int(row["octave"])
        step = 16 if e<34 else 1 << (e//2+1)
        k = (11*2**e//10)//step + int(row["offset"])
        check_proth(1+k*step,int(row["base"]))
    sequences = ([],[])
    threshold = 2**24
    with (HERE/"evidence/ratio-corner-primes.tsv").open() as f:
        for row in csv.DictReader(f,delimiter="\t"):
            e = threshold.bit_length()-1
            for j,place in enumerate(("lower","middle")):
                anchor = (threshold,57*threshold//50)[j]
                step = 16 if e<40 else 1 << (e//2+1+j)
                for name in ("same","opposite"):
                    key = place+"_"+name
                    p = 1+((anchor-1)//step+int(row[key+"_offset"]))*step
                    direction = check_proth(p,int(row[key+"_base"]))
                    assert direction == (name=="same")
                    sequences[direction].append(p)
            threshold = 13*threshold//10
    assert len(small)==514
    for sequence in sequences:
        ordered = sorted(sequence)
        assert len(set(ordered))==8252 and ordered[-1]>2**1585
        assert all(16*q < 19*p for p,q in zip(ordered,ordered[1:]))
    print("Independent GMP replay verified 514 midpoint primes and 16504 corner primes, "
          "their exact signs, and the narrow chain ratios in both common directions.")


if __name__=="__main__":
    minors()
    chains()
