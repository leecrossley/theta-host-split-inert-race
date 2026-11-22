"""Verify the ordinary-prime projection bound and the sharp finite counting cutoff."""
from collections import Counter, defaultdict
from fractions import Fraction as F

from dyadic_witnesses import OMIT, label, prime_divisors, trace
from phase_bridge import direction
from zero_strip import exp_lower


def pencil_bound():
    # Integrate the upper zero count against the decreasing common envelope.
    Q, b, nu = 7_000_000, F(29, 5), F(8, 5)
    R, lower_height = F(193, 20), 1095
    assert F(4, 1) / F(69, 100) < b
    assert F(68, 25)**17 < Q*b < exp_lower(F(88, 5))
    assert exp_lower(F(lower_height) / (R*F(88, 5))) > 625
    assert F(lower_height) / (R*17**2) < 2-nu
    # log(2*pi*e)>2.8, pi>3.14 and e lies strictly between 2.7 and 2.72.
    assert (2*F(314, 100)*F(27, 10))**5 > F(68, 25)**14
    log_qb = F(88, 5)
    integrated = (nu*b/((nu-1)*F(314, 100))*(log_qb-F(14, 5)+1/(nu-1))
                  + F(247, 1000)*(log_qb+1/nu) + F(6894, 1000))
    assert integrated < 93
    assert F(93, 1250) + F(3, 10**6) < F(3, 40)
    assert F(1, 4)-F(14, 160*lower_height) > F(249, 1000)
    assert (F(249, 2000)-F(3, 2)*F(3, 40))/6 == F(1, 500)
    # Strict overlap with the finite pattern chain, log(2) > 0.693.
    assert 1585*F(693, 1000) > lower_height
    print("From log X >= 1095 the paired-zero estimate gives absolute prime-weight error < 3/40 for all 24 twists.")
    print("Every real unit direction sees more than X/(500 log(2X)) ordinary primes in (X,2X).")


def losses(counts):
    survivors = [max(n-1, 0) for n in counts.values()]
    return sum(survivors) - max(survivors, default=0)


def finite_cutoff():
    end = 26000
    events = defaultdict(list)
    for p in prime_divisors(end):
        if p in OMIT or p % 16 not in (1, 15):
            continue
        a,b = trace(label(p))
        if a or b:
            d = direction(a,b)
            events[p].append((1,d))
            events[2*p].append((-1,d))
    counts = Counter()
    last = None
    for x in sorted(set(events) | {end}):
        if x > end:
            break
        for sign,d in events[x]:
            if sign == -1:
                counts[d] -= 1
        if losses(counts) < 2:
            last = (x, "point", sorted(counts.values()))
        for sign,d in events[x]:
            if sign == 1:
                counts[d] += 1
        if x < end and losses(counts) < 2:
            last = (x, "right", sorted(counts.values()))
    assert last == (431, "point", [0,1,2,4])
    # All sixteen sign classes include four, four, two and two ordinary directions.
    assert losses(Counter(direction(*trace(c)) for c in range(16) if any(trace(c)))) == 5
    assert 26000//2 > F(26211, 4)
    print("The final octave with fewer than two forced survivors after direction losses has left endpoint 431/2.")
    print("At that point the four direction counts are 0, 1, 2 and 4. All sixteen sign patterns guarantee five survivors.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    pencil_bound()
    finite_cutoff()
