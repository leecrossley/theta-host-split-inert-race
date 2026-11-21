"""Independently certify the extended profile from the prime sum and the real gamma integral."""
from flint import arb, ctx

DEN = 10**12
LEFT = 1008796429610
RIGHT = 3221800332553


class IntegralProfile:
    def __init__(self):
        self.a = arb(2).log()
        self.r = arb(2).sqrt()
        self.B = 2*self.r*(self.r-1)**2
        self.j0 = 3*self.a/2
        self.threshold = arb(19735)/10000
        self.margin = arb(1)/10**11
        self.constant = sum((q/arb.pi()).log() for q in (547,4376,8752,8752,727,5816,11632,11632))
        self.constant += 8*(arb(3)/4).digamma()
        self.powers = []
        for p in range(2,71):
            if any(p%d == 0 for d in range(2,p) if d*d <= p):
                continue
            uv = [1 if pow(p,(q-1)//2,q)==1 else -1 for q in (547,727)]
            n,k = p,1
            while n <= 70:
                multiple = 1 if p == 2 else (4 if pow(p,k,16) in (1,15) else 0)
                value = -multiple*sum(u**k for u in uv)
                if value:
                    self.powers.append((n, value*arb(p).log()))
                n,k = n*p,k+1
        assert sorted(n for n,_ in self.powers)==[4,16,31,47,49,64]
        assert (arb(RIGHT+1)/DEN+self.j0).exp()<71
        assert (-arb(LEFT-1)/DEN+self.j0).exp()<2
        self.branches = [self.coefficients(k) for k in (-2,-1,0,1,2)]
        self.trivial_terms = []
        for k in range(16):
            u = arb(3)/2+2*k
            G = (2*self.r*(self.a*u).cosh()-3)/(2*(self.r-1)**2*(u*u-arb(1)/4))
            weight = 16*G*(u*self.a/2).sinh()/self.a
            assert weight > 0
            self.trivial_terms.append((u,weight))

    def coefficients(self, region):
        # On each smooth branch j(s)=c exp(s/2)+d exp(-s/2).
        c,d = arb(0),arb(0)
        for shift,sign in ((1,1),(-1,-1)):
            side = 2*region+shift
            h = shift*self.a/2
            if -2 < side < 0:
                c += sign*((self.a+h)/2).exp()
                d -= sign*(-(self.a+h)/2).exp()
            elif 0 < side < 2:
                c -= sign*((h-self.a)/2).exp()
                d += sign*((self.a-h)/2).exp()
        return 2*c/(self.B*self.a),2*d/(self.B*self.a)

    def kernel(self, t):
        def g(v):
            return 4/self.B*((self.a-abs(v)).max(0)/2).sinh()
        return (g(t+self.a/2)-g(t-self.a/2))/self.a

    def prime(self, x):
        def score(y):
            value = arb(0)
            for n,weight in self.powers:
                z=y/n
                value += weight*(2*(1-z)).min(4*z-1).max(0)
            return value/(self.B*y.sqrt())
        y=x.exp()/(2*self.r)
        return (score(2*y)-score(y))/self.a

    def left_value(self, x):
        # Exactly integrate the four noncentral pieces after w=exp(-t/2).
        half=(x/2).exp()
        edges=[(-v/2).exp() for v in (self.j0-x,x-self.a/2,x+self.a/2,x+self.j0)]
        edges.append(arb(0))
        assert edges[0]<1
        jx=self.kernel(x)
        def first(w):
            return 2*jx*(w.atan()-w.log1p())
        integral=first(arb(1))-first(edges[0])
        for region,hi,lo in zip((3,2,1,0),edges,edges[1:]):
            c,d=self.branches[region]
            A,B=d/half,c*half
            def primitive(w):
                return ((A-B)/2*(w*w).log1p() + (jx-(A+B)/2)*(-w).log1p()
                        + (-jx-(A+B)/2)*w.log1p() + 2*jx*w.atan())
            integral+=primitive(hi)-primitive(lo)
        return self.prime(x)+self.constant*jx-8*integral

    def right_lower(self, x):
        return self.prime(x)+sum(weight*(-u*x).exp() for u,weight in self.trivial_terms)

    def right_tail(self,x):
        u=arb(3)/2+2*len(self.trivial_terms)
        d=x-self.j0
        return 8*self.r*(-u*d).exp()/(2*self.a*(self.r-1)**2*(u*u-arb(1)/4)*(1-(-2*d).exp()))

    def check_kink(self):
        cells=[]
        for start,stop,step in ((1,128,1),(128,1000,2),(1000,3000,5),(3000,8000,20)):
            for low in range(start,stop,step):
                high=min(low+step,stop)
                t=arb(low+high,high-low)/2000
                cells.append((t,(high-low)*(-t/2).exp()/(2000*t.sinh())))
        assert 4*(self.a/2).cosh()/(self.B*self.a)<14
        assert 8*(self.a/2).sinh()/(self.B*self.a)<9
        error=8*(arb(14)/1000+36*arb(-12).exp()/((arb(3)/2)*(1-arb(-16).exp())))
        for low in range(2070,2080):
            x=arb(2*low+1,1)/4000
            jx=self.kernel(x)
            integral=sum(w*(self.kernel(x-t)+self.kernel(x+t)-2*jx) for t,w in cells)
            assert self.prime(x)+self.constant*jx-8*integral-error > -self.threshold+self.margin


def certify(function,lo,hi,den,threshold,depth=0):
    value=function(arb(lo+hi,hi-lo)/(2*den))
    if value > threshold:
        return
    if depth>=65:
        raise ArithmeticError(f"Inconclusive cell with numerator endpoints {lo}, {hi}")
    certify(function,2*lo,lo+hi,2*den,threshold,depth+1)
    certify(function,lo+hi,2*hi,2*den,threshold,depth+1)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=192
    p=IntegralProfile()
    target=-p.threshold+p.margin
    assert arb(1035)/1000 < p.j0 < arb(1040)/1000
    certify(p.left_value,LEFT,1035*10**9,DEN,target)
    p.check_kink()
    certify(p.right_lower,1040*10**9,RIGHT,DEN,target)
    print("Real integral, ten breakpoint cells and positive trivial terms independently certify the full interval.")
    assert p.left_value(arb(LEFT-1)/DEN) < -p.threshold-p.margin
    x=arb(RIGHT+1)/DEN
    assert p.right_tail(x)<arb(1)/10**25
    assert p.right_lower(x)+p.right_tail(x) < -p.threshold-p.margin
    print("Independent exterior checks put both neighbouring points strictly outside the superlevel set.")


if __name__ == "__main__":
    main()
