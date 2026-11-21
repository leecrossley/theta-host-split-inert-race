"""Independent original-tent/Arb replay of dilated contacts and deep height brackets."""
import csv
from pathlib import Path

from flint import arb, ctx, fmpz
from crossing_recheck import TOP, field_data, sign

HERE = Path(__file__).resolve().parent
DECIMALS = 72
UNIT = 10**DECIMALS


def rows(name):
    with (HERE/"evidence"/name).open() as file:
        return list(csv.DictReader(file,delimiter="\t"))


def refine(data,catalog):
    refined = []
    target = arb(3947)/2000
    for row in catalog:
        _,move = data["AB".index(row["host"])]
        left,right = (int(row[k]) for k in ("left_eighths","right_eighths"))
        lo,hi = arb(left)/8,arb(right)/8
        gl,gr = lo.sqrt()*move(left,8),hi.sqrt()*move(right,8)
        slope = (gr-gl)/(hi-lo)
        intercept = gl-slope*lo
        orient = int(row["orientation"])
        threshold = int(row["sign"])*target
        start = int(row["lower_billionths"])
        a,b = start*(UNIT//10**9),(start+1)*(UNIT//10**9)
        while b-a>1:
            k = (a+b)//2
            y = arb(k)/UNIT
            if sign(slope*y+intercept-threshold*y.sqrt()) == orient:
                b = k
            else:
                a = k
        assert sign(move(a,UNIT)-threshold)==-orient
        assert sign(move(b,UNIT)-threshold)==orient
        refined.append((row["host"],row["sign"],row["left_eighths"],str(a)))
    return refined


def verify_factor(data):
    records = rows("factor-200-contacts.tsv")
    expected = {(r["first"],r["second"],int(r["relative_sign"]),int(r["left_1600ths"]),
                 int(r["right_1600ths"])):r for r in records}
    found = set()
    for i,(early_knots,first) in enumerate(data):
        for j,(late_knots,second) in enumerate(data):
            knots = sorted({1600,TOP,*(200*k for k in early_knots if 1600<=200*k<=TOP),
                            *(k for k in late_knots if k>=1600)})
            l = knots[0]
            first_l,second_l = first(l,1600),second(l,8)
            for r in knots[1:]:
                first_r,second_r = first(r,1600),second(r,8)
                for eps in (1,-1):
                    left_value = first_l-eps*second_l
                    right_value = first_r-eps*second_r
                    if sign(left_value)==sign(right_value):
                        continue
                    yl,yr = arb(l)/1600,arb(r)/1600
                    gl,gr = yl.sqrt()*left_value,yr.sqrt()*right_value
                    y = yl-gl*(yr-yl)/(gr-gl)
                    assert yl<y<yr
                    key = ("AB"[i],"AB"[j],eps,l,r)
                    record = expected[key]
                    n,m = int(record["height_floor_1e12"]),int(record["amplitude_floor_1e12"])
                    assert n<y*10**12<n+1
                    sg = int(record["first_sign"])
                    if sg:
                        v,w = yl.sqrt()*first_l,yr.sqrt()*first_r
                        value = (v+(w-v)*(y-yl)/(yr-yl))/y.sqrt()
                        assert sign(value)==sg and m<abs(value)*10**12<m+1
                    else:
                        # The first true knot begins this field's support.
                        assert r<=200*early_knots[1] and first_l.contains(0) and first_r.contains(0)
                    found.add(key)
                l,first_l,second_l = r,first_r,second_r
    assert found==set(expected) and len(found)==70
    print("Arb at 416 bits verified all 70 entries of the four ordered factor-200 equations using the original tents at every joint endpoint.")


def main():
    ctx.prec = 416
    ps = [n for n in range(3,TOP,2) if fmpz(n).is_prime()]
    data = [field_data(ps,h) for h in (0,1)]
    refined = refine(data,rows("exact-crossings.tsv"))
    table = rows("height-brackets-72.tsv")
    assert refined==[(v["host"],v["sign"],v["left_eighths"],v["height_floor_1e72"]) for v in table]
    print("Independent prefix sums and signed bisection certified 72 decimal places for all 62 crossing heights.")
    verify_factor(data)


if __name__=="__main__":
    main()
