"""Field data, shared zero certificate, and the explicit prime-only remainder."""
from hashlib import sha256
from pathlib import Path
from flint import arb,ctx,fmpz

ROOT=Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    z=(ROOT/'evidence/dyadic-zeros.tsv').read_bytes()
    assert sha256(z).hexdigest()=='8714dca2f130c861932c69fba5b4962626bd5ebdd77993ae1e8fd7e9330cc7da'
    assert len(z.splitlines())==14000
    for p,q,s,t in ((547,727,343057911569,343070636977),(571,599,413718278401,603339156001)):
        assert all(fmpz(l).is_prime() for l in (p,q,s,t))
        assert s%(16*p*q)==t%(16*p*q)==1
        assert max(s,t)<arb(28).exp()
        assert max(p,q)**8<arb(53).exp() and 16*p*q<=6362704
    x=arb(100)
    err=2*10**6*(-x/4).exp()+128*(-(x-2)/2).exp()/(1-(-(x-2)).exp())+22000*(-3*x/20).exp()
    assert err<arb(7)/1000
    assert 145*arb(2).log()>100 and 144*arb(2).log()<100
    assert 2*arb(100).exp()*arb(4)**4095<arb(10)**2510
    print("The two host frames, four split type primes, the published zero table and the prime-only error below 0.007 from log centre 100 are verified.")


if __name__=="__main__":
    main()
