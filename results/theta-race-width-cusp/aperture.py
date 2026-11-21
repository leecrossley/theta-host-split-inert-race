"""Read one extremum at each isolated prime power, including zero traces."""
from fractions import Fraction
from flint import arb, ctx
from width_kernel import EVIDENCE, PrimeTable, Window, ball, check_inputs, trace
import sys

CASES = (2,3,5,7,11,17,31,47,547,727)
DEN = 10**32


def main():
    check_inputs()
    ctx.prec = 328
    table = PrimeTable(727**2+2)
    lines = []
    for p in CASES:
        for j in (1,2):
            n = p**j
            b = Fraction(n)
            z = Fraction(2*n+1,2*n)
            window = Window(b,table)
            _, result, background = window.measurements(z)
            unit = 4*arb(p).log()*(ball(z).log()/2).sinh()/window.c
            floor = (result*DEN).floor().unique_fmpz()
            assert floor is not None
            enc = arb(2*floor+1,1)/(2*DEN)
            recovered = (enc-background)/unit
            u = trace(p,j)
            assert recovered.unique_fmpz() == u and abs(recovered-u) < arb(10)**-17
            lines.append(f"{p}\t{j}\t{floor}\t{u}")
    data = "\n".join(lines)+"\n"
    if sys.argv[1:] == ["--print-table"]:
        print(data,end="")
        return
    assert not sys.argv[1:] and (EVIDENCE/"weighted-extrema.tsv").read_text() == data
    print("20 weighted extrema at ten primes independently isolate the indicated first and second traces.")
    print("Each archived limiting value has an outward enclosure of width 10^-32. Both ramified odd primes and two are included.")


if __name__ == "__main__":
    main()
