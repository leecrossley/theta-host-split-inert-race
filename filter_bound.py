"""Certify a positive compact averaging measure and the finite-height bridge."""
import json
from pathlib import Path
from flint import arb, ctx
from check_spectrum import read_brackets

ROOT = Path(__file__).resolve().parent


def prime_error(x):
    return (2*10**6*(-x/4).exp() + 128*(-(x-2)/2).exp()/(1-(-(x-2)).exp())
            + 22000*(-3*x/20).exp())


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    bounds = json.loads((ROOT/"evidence"/"filter.json").read_text())
    chosen = bounds["carrier_row"]-1
    assert bounds["complete_through"] == 48 and 0 <= chosen < 6000
    box = arb(bounds["box_32"])/32
    support = [(arb(k)/16, weight) for k,weight in bounds["atoms_16"]]
    assert all(weight > 0 and u >= 0 for u,weight in support) and box > 0
    total = sum(weight for _,weight in support)
    radius = max(u for u,_ in support) + 2*box
    a, root = arb(2).log(), arb(2).sqrt()
    c = 2*(root-1)**2
    rows = read_brackets()
    gamma = rows[chosen][2]
    period = 2*arb.pi()/gamma
    divisor = arb(10)**bounds["deadline_power"]
    assert 2*radius + period + a < divisor.log()
    heights = [arb(x)/100 for x in (6050,5840,8366,8800,10910,11144)]
    bridge = [arb(0) for _ in heights]
    leakage, carrier = arb(0), arb(0)
    for j,(_,_,t) in enumerate(rows):
        g = (3-2*root*(a*t).cos())/(c*(t*t+arb(1)/4))
        pair = 4*g*(a*t/2).sin()/a
        transfer = (box*t).sinc()**2 * sum(weight*(u*t).cos() for u,weight in support)/total
        if j == chosen:
            carrier = abs(pair*transfer)
        else:
            leakage += abs(pair*transfer)
        for k,x in enumerate(heights):
            bridge[k] -= pair*(t*x).sin()
    tail = (arb(144)/100)/(box*bounds["complete_through"])**2
    assert carrier > arb(4917)/1000 and leakage < arb(285)/100
    assert prime_error(arb(100)) < arb(7)/1000
    assert 2*radius == arb(155)/4
    assert carrier-leakage-tail-prime_error(arb(100)) > arb(19735)/10000
    for i,(value,x) in enumerate(zip(bridge,heights)):
        direction = 1 if i%2 == 0 else -1
        assert x > 55 and direction*value-arb(144)/100-prime_error(x) > 20
    for i,(bottom,top) in enumerate(((9*arb(10).log(),arb(55)),(arb(55),arb(80)),(arb(80),arb(100)))):
        for x in heights[2*i:2*i+2]:
            assert x-3*a/2 > top and x-a/2 < bottom+divisor.log()
    print(f"Positive measure with {len(support)} atom pairs, total integer weight {total}.")
    print("Carrier amplitude > 4.917; all other retained amplitudes total < 2.85.")
    print("Carrier minus all other listed pairs, unlisted modes and prime error exceeds 1.9735.")
    print(f"Two support radii, one carrier period and log 2 fit strictly within log(10^{bounds['deadline_power']}).")
    print("Six fixed prime-score moves clear normalized magnitude 20 and bridge 10^9 through exp(100).")


if __name__ == "__main__":
    main()
