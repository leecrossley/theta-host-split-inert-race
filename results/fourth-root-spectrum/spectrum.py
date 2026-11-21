"""Outward bounds for the fourth-root reconstruction and the smaller zero tests."""
from flint import arb, ctx


def zero_tail(t, q):
    return 80*((2*q*t).log()+1)/(t-2)


def arithmetic_cost():
    # Both moduli are bounded above by the first one. The zero-free theorem
    # is applied to each field's own product of characters.
    x = arb(12493)
    q = max(16*547*727, 16*571*599)
    b, c = arb(10).exp(), arb(25).exp()
    low = 2*10**6*(-x/(20*(q*b).log())).exp()
    middle = 2000*zero_tail(b,q)*(-x/(20*(q*c).log())).exp()
    high = 2000*zero_tail(c,q)
    possible_real = 20/(arb(1).exp()*x)
    trivial = 128*(-(x-2)/2).exp()/(1-(-(x-2)).exp())
    finite_and_higher = 10**6*((-x/10).exp()+(-x/3).exp())
    return low + middle + high + possible_real + trivial + finite_and_higher


def smoothing_costs():
    x = arb(12494)
    # The x-dependent terms decrease beyond this endpoint.
    physical = 10**5*(arb(1)/10**12 + (x+3)*(-x/4).exp()/10**6
                       + (x+3)*(-x/6).exp())
    tail = 10**6*(x/2-x*arb(4).log()).exp()
    radius = (3*(x+2)).sqrt()*(-x/4).exp()/10**6
    cutoff = arb(100).exp()
    projected = radius*(cutoff+1)*10**6 + 100*zero_tail(cutoff,11632)
    assert radius < (-arb(3000)).exp()
    assert physical < arb(11)/10**8
    assert tail < arb(1)/10**10 and projected < arb(1)/10**10
    assert arithmetic_cost() < arb(79)/10**5
    assert physical+tail+arithmetic_cost() < arb(8)/10**4
    return physical, tail, projected


def finite_alternatives(physical, tail, projected):
    cut = arb(10)**12
    q = 11632
    a = (2*q*cut).log()+1
    band = zero_tail(arb(8590),q)
    cases = (
        (1620548, 1620550, 4, 22, 405161),
        (1620548, 1177000, 4, 12, 405161),
        (3895692, 3895692, 2, 30, 973947),
    )
    for reference, ceiling, units, budget, top_log in cases:
        x = arb(ceiling+1)  # all points of the averaging support
        eps = arb(units)/10**8
        p = arb(49)/50 * x/reference
        assert 0 < p < 1
        log_t = arb(ceiling)/4 + arb(8*10**6).log() + arb(ceiling+1).log()/2
        assert log_t < top_log
        low = (2000*eps*(eps*x).exp()+50*((eps*x).cosh()-1))*band
        high = (1000*(eps*x).exp()+50)*160/cut*(a/(1-p)+p/(1-p)**2)
        # The projected unaveraged tail past T is harmless after the cutoff.
        leftover = 50*zero_tail(arb(top_log).exp(),q)
        assert low+high+leftover+tail+projected+physical < arb(budget)/1000
    assert 405161 < 810330/2 and 973947 < 1947900/2


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 160
    physical, tail, projected = smoothing_costs()
    finite_alternatives(physical, tail, projected)
    print("The prime-side smoothing error is below 0.00000011. The unrestricted high odd spectrum costs less than 10^-10.")
    print("The base, trivial-zero, excluded-prime and higher-power cost is below 0.00079 for both fields throughout the required support.")
    print("The finite weighted odd-zero sum approximates every unsmoothed normalized move to absolute error below 0.0008 from log centre 12,494.")
    print("The first-400H region ends at exp(405161). The second window's region ends at exp(973947). The three projected error budgets survive the smoothing.")


if __name__ == "__main__":
    main()
