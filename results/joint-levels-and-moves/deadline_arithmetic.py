"""Certified physical-space comparison bounds for both prime races."""
from functools import lru_cache

from flint import arb

BAD=((2,547,727,343057911569,343070636977),(2,571,599,413718278401,603339156001))
TIERS=(0,24,32,40,48,55,70,100)


def logfloor(tier):
    assert tier in TIERS
    return 9*arb(10).log() if tier==0 else arb(tier)


@lru_cache(maxsize=None)
def costs(field,tier):
    a,r=arb(2).log(),arb(2).sqrt()
    B=2*r*(r-1)**2
    ell=logfloor(tier)
    Z=(ell/2).exp()
    degree=arb(816)/100/B*(2/Z).sqrt()
    # Combine the finite square-root estimate (Theorem 2) with Table 1 for large u.
    # Both give |theta(u;16,a)-u/8| <= eps*u/8 for every u >= Z.
    eps=(8*arb(2072)/1000/Z.sqrt()).max(arb(9)/1000) if Z < 10**10 else arb(9)/1000
    area=6*r-2*(1+r)*(2*r).sqrt()
    square=2+8*eps/B
    move_square=(2*area+2*eps*(2+r))/B
    l4=ell+arb(4).log()
    higher=arb(816)/100/B*(sum((l4/k-ell/2).exp() for k in (3,4))+
                            l4/a*(l4/5-ell/2).exp())
    bad=8*(3*a+sum(arb(p).log() for p in BAD[field][1:]))/(B*Z)
    level=(square+degree-1).max(arb(1))+higher+bad
    move=(move_square+degree+2*(higher+bad))/a
    return level,move
