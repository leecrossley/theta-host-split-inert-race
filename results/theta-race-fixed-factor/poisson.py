"""A first-derivative check of both inverse-square bounds from the positive Poisson sum."""
from flint import acb, acb_series, arb, ctx, dirichlet_char

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def main():
    ctx.prec=160
    sigma=arb(1)/1000
    g=arb(6574356596173,1)/(2*10**12)
    for t,m,cap in ((arb(0),0,937),(g,1,217)):
        s=acb(arb(1)/2+sigma,t)
        total=arb(0)
        for p in (547,727):
            for k,b in ((1,1),(8,5),(16,5),(16,13)):
                q=p*k
                label=next(a for a in range(p-1,q,p) if a%k==b%k)
                chi=dirichlet_char(q,label)
                v=acb_series([s,1],2).dirichlet_l(chi)
                total+=((q/arb.pi()).log()/2+((s+1)/2).digamma()/2+v[1]/v[0]).real
        weight=total/sigma-m/sigma**2
        assert 0<weight and sigma**2*weight<1
        assert weight/(1-sigma**2*weight)<cap, (m,weight)
        print(f"regularized Poisson sum at real displacement 1/1000 independently gives inverse-square sum < {cap}")


if __name__=="__main__":
    main()
