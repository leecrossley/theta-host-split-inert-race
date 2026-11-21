"""Bridge the two degree-one prime sets to the upper progression estimate."""
import csv
from fractions import Fraction
from pathlib import Path

from height_ratios import HOSTS, is_prime, prime_data
from verify_rank import primes_through

HERE = Path(__file__).resolve().parent
RATE = Fraction(19, 16)
EXCLUDED = {2, 547, 727, 571, 599, 343057911569, 343070636977,
            413718278401, 603339156001}


def bridges():
    by_host = {host: [] for host in HOSTS}
    with (HERE/"evidence/prime-bridges.tsv").open() as source:
        for row in csv.DictReader(source,delimiter="\t"):
            e,j = int(row["octave"]),int(row["bin"])
            step = 16 if e < 34 else 1 << (e//2+1)
            for host in HOSTS:
                k = (j << e)//(5*step) + int(row[host+"_offset"])
                by_host[host].append(1+k*step)
    # The main release replay certifies these 5,140 preceding witnesses.
    with (HERE/"evidence/ratio-primes.tsv").open() as source:
        added = list(csv.DictReader(source,delimiter="\t"))
    assert [int(row["octave"]) for row in added] == list(range(18,532))
    small = 0
    for row in added:
        e = int(row["octave"])
        step = 16 if e < 34 else 1 << (e//2+1)
        k = ((11 << e)//10)//step + int(row["offset"])
        p = 1+k*step
        witness = int(row["base"])
        if e < 34:
            assert witness == 0 and is_prime(p)
            small += 1
        else:
            assert 0 < k < step and k % 2 and 1 < witness < p
            assert pow(witness,(p-1)//2,p) == p-1
        assert p not in EXCLUDED and p % 16 == 1
        for host in HOSTS:
            assert prime_data(host,p,False) == prime_data(host,p,True)
            assert prime_data(host,p,True)[0] == p
            by_host[host].append(p)
    assert small == 16
    primes = [p for p in primes_through(320000) if p not in EXCLUDED]
    final = []
    for host,large in by_host.items():
        initial = [p for p in primes if prime_data(host,p,True)[0] == p]
        chain = sorted(set(initial+large))
        final.append([(p,q) for p,q in zip([1]+chain,chain) if 16*q >= 19*p][-1])
        assert max(large) > 2**531
        assert min(large) < max(initial)
    assert final == [(479,641),(479,577)]
    # For t >= 366 the diagonal log mass is >= Delta/8 - 1.90 (U+L)/t.
    assert Fraction(1,8) > Fraction(190,100) * Fraction(35,3) / 366
    assert RATE**4 < 2 and RATE < 2 and 2*RATE < 4
    print("The 514 additional midpoint primes pass 16 trial divisions and 498 Proth tests. "
          "Both degree-one chains have every ratio below 19/16 after their last gaps "
          "479 to 641 and 479 to 577. The analytic mass is positive from log L = 366.")


def joint_corners():
    # A nonzero diagonal trace pair is on one of two lines. Two such primes
    # have independent trace vectors precisely when they lie on opposite lines.
    previous = [1, 1]
    final_bad = [None, None]
    for p in primes_through(1 << 24):
        if p % 16 not in (1, 15) or p in EXCLUDED:
            continue
        trace = [prime_data(host,p,True) for host in HOSTS]
        if any(norm != p for norm,sign in trace):
            continue
        direction = int(trace[0][1] == trace[1][1])
        if 16*p >= 19*previous[direction]:
            final_bad[direction] = (previous[direction],p)
        previous[direction] = p
    assert final_bad == [(2609,3313),(2351,2879)]
    counts = [0,0]
    lower = 1 << 24
    with (HERE/"evidence/ratio-corner-primes.tsv").open() as source:
        rows = list(csv.DictReader(source,delimiter="\t"))
    assert len(rows) == 4126
    for i,row in enumerate(rows):
        assert int(row["bin"]) == i
        e = lower.bit_length()-1
        for which,anchor in (("lower",lower),("middle",57*lower//50)):
            step = 16 if e < 40 else 1 << (e//2 + (2 if which == "middle" else 1))
            for direction,name in enumerate(("opposite","same")):
                key = which+"_"+name
                k = (anchor-1)//step + int(row[key+"_offset"])
                p = 1+k*step
                assert p not in EXCLUDED and p % 16 == 1
                data = [prime_data(host,p,True) for host in HOSTS]
                assert all(norm == p for norm,sign in data)
                assert int(data[0][1] == data[1][1]) == direction
                assert previous[direction] < p and 16*p < 19*previous[direction]
                witness = int(row[key+"_base"])
                if e < 40:
                    assert witness == 0 and is_prime(p)
                    counts[0] += 1
                else:
                    assert 0 < k < step and k % 2 and 1 < witness < p
                    assert pow(witness,(p-1)//2,p) == p-1
                    counts[1] += 1
                previous[direction] = p
        lower = 13*lower//10
    assert min(previous) > 2**1585 and sum(counts) == 16504
    assert Fraction(6,1097) < Fraction(1,6) < RATE-1
    print(f"The two common degree-one directions have last deficient gaps "
          f"2609 to 3313 and 2351 to 2879 at factor 19/16. "
          f"The joint chains certify {sum(counts)} primes "
          f"({counts[0]} trial divisions and {counts[1]} Proth tests).")


if __name__ == "__main__":
    bridges()
    joint_corners()
