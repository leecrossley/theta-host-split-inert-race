"""Certify the threshold component by summing the archimedean term in closed form."""
from fractions import Fraction

from flint import arb, ctx, fmpq

SCALE = 10**12
LEFT = 1008796429610
RIGHT = 3221800332553
POWERS = ((4, 2, -2), (16, 2, -2), (31, 31, 8),
          (47, 47, -8), (49, 7, -8), (64, 2, -2))


def ball(q):
    return arb(fmpq(q.numerator, q.denominator))


def knot_term(distance):
    # The defining series is positive and decreases with distance.
    # Its continuous value at distance zero is log 2.
    def endpoint(d):
        if d.is_zero():
            return arb(2).log()
        q = (-d).exp()
        one_minus_q = -(-d).expm1()
        return ((1+q)*q.log1p() + one_minus_q*one_minus_q.log()) / (2*q.sqrt())

    d = abs(distance)
    return endpoint(d.lower().max(0)).union(endpoint(d.upper()))


class Profile:
    def __init__(self):
        self.a = arb(2).log()
        self.r = arb(2).sqrt()
        self.c = 2*(self.r-1)**2
        self.B = self.r*self.c
        self.A = sum((p*m/arb.pi()).log() for p in (547,727) for m in (1,8,16,16))
        self.A -= 8*(arb.const_euler()+self.a)
        self.terms = [(arb(n).log(), m*arb(p).log()/arb(n).sqrt()) for n,p,m in POWERS]
        assert (arb(RIGHT+1)/SCALE+3*self.a/2).exp() < 71
        assert (-arb(LEFT-1)/SCALE+3*self.a/2).exp() < 2

    def kernel(self, x):
        def g(v):
            return 4/self.B*((self.a-abs(v)).max(0)/2).sinh()
        return (g(x+self.a/2)-g(x-self.a/2))/self.a

    def value(self, x):
        prime = sum(weight*self.kernel(x-logn) for logn,weight in self.terms)
        far_right, near_right, near_left, far_left = (
            knot_term(x + k*self.a/2) for k in (3,1,-1,-3))
        jumps = self.r*(far_right-far_left) - (3+self.r)*(near_right-near_left)
        return prime+self.A*self.kernel(x)-8*jumps/(self.a*self.c)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 192
    profile = Profile()
    threshold = arb(19735)/10000
    margin = arb(1)/10**11
    left, right = Fraction(LEFT,SCALE), Fraction(RIGHT,SCALE)
    todo = [(left,right,0)]
    while todo:
        lo,hi,depth = todo.pop()
        box = ball(lo).union(ball(hi))
        if profile.value(box) > -threshold+margin:
            continue
        if depth >= 65:
            raise ArithmeticError(f"Undecided cell [{lo}, {hi}]")
        mid = (lo+hi)/2
        todo.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
    print("Four-knot formula: U + 1.9735 > 10^-11 throughout [1.008796429610, 3.221800332553].")
    endpoint_bounds = ((LEFT-1,-313,-312),(LEFT,81,82),
                       (RIGHT,40,41),(RIGHT+1,-65,-64))
    for index,lo,hi in endpoint_bounds:
        scaled = SCALE*(profile.value(arb(index)/SCALE)+threshold)
        assert lo < scaled < hi
        print(f"x = {index}/10^12: {lo} < 10^12 (U + 1.9735) < {hi}.")
    cin = 2*ball(right-left).exp()
    cout = 2*ball(right-left+Fraction(2,SCALE)).exp()
    # Outward decimal bounds for the component's factor and a strict lower bound for the certified factor.
    print("C = 2 exp(2.213003902943)")
    assert arb("18.286280578160") < cin < arb("18.286280578161") < cout < arb("18.286280578198")
    print("18.286280578160 < certified factor < 18.286280578161; component upper bound < 18.286280578198.")


if __name__ == "__main__":
    main()
