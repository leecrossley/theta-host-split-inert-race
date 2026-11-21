"""Replay the filter as a Laurent polynomial and cover the missing zeros at another displacement."""
from pathlib import Path
import json
from flint import acb, acb_series, arb, ctx, dirichlet_char
from check_spectrum import FAMILY

ROOT = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 192
    data = json.loads((ROOT/"evidence"/"filter.json").read_text())
    assert data["complete_through"] == 48
    b = arb(data["box_32"])/32
    freq = data["carrier_row"]-1
    weights = dict(data["atoms_16"])
    assert 0 <= freq < 6000 and len(weights) == len(data["atoms_16"])
    assert all(isinstance(k,int) and k >= 0 and isinstance(w,int) and w > 0 for k,w in weights.items())
    denominator = sum(weights.values())
    block_size = 12
    blocks = [(start, [weights.get(k,0) for k in range(start+block_size-1,start-1,-1)])
              for start in range(0,max(weights)+1,block_size)]
    table = []
    for line in (ROOT/"evidence"/"zero-brackets.tsv").read_text().splitlines():
        _,_,left,right = map(int,line.split())
        table.append((arb(left)/10**14).union(arb(right)/10**14))
    assert len(table) == 6000
    for x in range(0,49,3):
        s = acb(arb(3)/2,x)
        positive_mass = arb(0)
        for q,label in FAMILY:
            local = acb_series([s,1],2).dirichlet_l(dirichlet_char(q,label))
            positive_mass += ((q/arb.pi()).log()/2 + ((s+1)/2).digamma()/2 + local[1]/local[0]).real
        used = sum(1/(1+(t-x)**2)+1/(1+(t+x)**2) for t in table)
        assert 0 < positive_mass-used < arb(1)/5 < 1/(1+(arb(3)/2)**2)
    print("Seventeen unit-displacement Poisson checks independently cover all ordinates through 48.")
    a = arb(2).log()
    root = arb(2).sqrt()
    wanted, other = arb(0), arb(0)
    for j,t in enumerate(table):
        z = acb(0,t/16).exp()
        value = acb(0)
        for start, coeff in blocks:
            block = acb(0)
            for k in coeff:
                block = block*z+k
            value += block * acb(0,start*t/16).exp()
        triangle = (1-(2*b*t).cos())/(2*(b*t)**2)
        half_sine = (a*t/2).sin()
        kernel = ((root-1)**2 + 4*root*half_sine**2)/(2*(root-1)**2*(t*t+arb(1)/4))
        effect = abs(4*kernel*half_sine*triangle*value.real/(a*denominator))
        if j == freq:
            wanted = effect
        else:
            other += effect
    omitted = (arb(152)/100)/(48*b)**2
    # At x >= 100 the three explicit prime-only remainders total less than 0.007.
    endpoint = arb(100)
    error = (2000000*(-endpoint/4).exp() + 128*(-(endpoint-2)/2).exp()/(1-(-(endpoint-2)).exp())
             + 22000*(-3*endpoint/20).exp())
    assert error < arb(7)/1000
    assert wanted - other - omitted - error > arb(19735)/10000
    extent = arb(max(weights))/16+2*b
    assert (2*extent+2*arb.pi()/table[freq]).exp()*2 < arb(10)**data["deadline_power"]
    print("Laurent-polynomial evaluation and the independent 1.52 tail bound certify both signs and the deadline.")


if __name__ == "__main__":
    main()
