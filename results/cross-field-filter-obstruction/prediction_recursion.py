"""A second proof by rational Toeplitz prediction recursion and a row-sum perturbation."""
from flint import acb, arb, ctx, fmpq
from radii import BOUNDS, brackets


def certify(moments, field, kind):
    scale = 2**55
    integer = [int(float(val.mid())*scale) for val in moments]
    # Each entry in the true Toeplitz matrix is within err[j] of integer[j]/scale.
    err = [abs(scale*val-i) for val, i in zip(moments, integer)]
    margin = int((err[0]+2*sum(err[1:])).upper().ceil().unique_fmpz())
    integer[0] -= margin
    row = [fmpq(i,scale) for i in integer]
    variance = row[0]
    predictors = []
    for m in range(1, len(row)+1):
        assert variance > 0
        if m in BOUNDS:
            threshold = fmpq(BOUNDS[m][2*field+kind],1000)
            assert variance > threshold**2, (field,kind,m,float(variance))
        if m == len(row):
            break
        reflection = (row[m]-sum(b*row[m-1-j] for j,b in enumerate(predictors)))/variance
        predictors = [b-reflection*c for b,c in zip(predictors,reversed(predictors))]+[reflection]
        variance *= 1-reflection**2


def main():
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    ctx.prec = 224
    n = max(BOUNDS)
    a = arb(2).log()
    s = arb(2).sqrt()
    for field in (0, 1):
        moments = [[arb(0) for _ in range(n)] for _ in range(2)]
        for lo, hi in brackets(field):
            middle = arb(lo+hi)/(2*10**14)
            radius = arb(hi-lo)/(2*10**14)
            t = arb(lo+hi,hi-lo)/(2*10**14)
            half_sine = (t*a/2).sin()
            level = (1+4*s*half_sine**2/(s-1)**2)/(2*(arb(1)/4+t*t))
            weights = [2*level**2, 8*(level*half_sine/a)**2]
            unit = acb(0,a*middle).exp()
            phase = acb(1)
            for j in range(n):
                cosine = phase.real + arb(0,j*a*radius)
                for kind in (0,1):
                    moments[kind][j] += weights[kind]*cosine
                phase *= unit
        for kind, row in enumerate(moments):
            certify(row, field, kind)
    print('Exact rational prediction recursion, with midpoint phases and an outward row-sum correction, verifies all 32 lower bounds.')


if __name__ == '__main__':
    main()
