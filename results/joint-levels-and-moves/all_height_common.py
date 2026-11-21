"""Enclose all mixed harmonics in two common positive probability tests."""
from fractions import Fraction
from flint import arb,ctx

from deadline_arithmetic import costs
from deadline_spectrum import ordinate,tables

CARRIERS=((0,11632,1453,236949851779030,236949851779034),
          (1,9136,5709,231411042788551,231411042788555))
M=8
BOX=3500


def verify(method):
    ctx.prec=144 if method=='sinc' else 160
    data=tables()
    assert all(data.count(r)==1 for r in CARRIERS)
    g,h=map(ordinate,CARRIERS)
    a,r=arb(2).log(),arb(2).sqrt()
    def response(t):
        if method=='sinc':
            return (BOX*t).sinc()**2
        return abs((1-(2*BOX*t).cos())/(2*(BOX*t)**2))
    shifts=[]
    mass_error=arb(0)
    for k in range(-M,M+1):
        for l in range(-M,M+1):
            nu=k*g+l*h
            weight=arb(M+1-abs(k))*(M+1-abs(l))/(M+1)**2
            shifts.append((nu,weight,(k,l)))
            if k or l:
                assert abs(nu)>0
                mass_error+=weight*response(nu)
    assert M*(g+h)<38 and mass_error<arb(1)/100
    missing=arb(3)/10*(M+1)**2/(BOX*(48-M*(g+h)))**2
    for i in (0,1):
        main_level=main_move=lost_level=lost_move=arb(0)
        for row in data:
            if row[0]!=i:
                continue
            t=ordinate(row)
            u=a*t/2
            num=(3-2*r*(2*u).cos() if method=='sinc' else (r-1)**2+4*r*u.sin()**2)
            level=num/(2*(r-1)**2*(arb(1)/4+t*t))
            move=abs(2*level*u.sin()/a)
            leaked=arb(0)
            for nu,weight,harmonic in shifts:
                if row==CARRIERS[i] and harmonic==((1,0) if i==0 else (0,1)):
                    main_level,main_move=2*weight*level,2*weight*move
                else:
                    delta=t-nu
                    assert abs(delta)>0
                    leaked+=weight*response(delta)
            lost_level+=2*level*leaked
            lost_move+=2*move*leaked
        al,ad=costs(i,100)
        level_amp=main_level/(1+mass_error)
        move_amp=main_move/(1+mass_error)
        el=(lost_level+18*missing+arb(1)/1000)/(1-mass_error)+al
        ed=(lost_move+50*missing+arb(1)/1000)/(1-mass_error)+ad
        fraction=1-((el/level_amp)**2+(ed/move_amp)**2).sqrt()
        assert fraction*level_amp>arb(117)/100 and fraction*move_amp>arb(245)/100, (i,fraction*level_amp,fraction*move_amp)
        assert move_amp-ed>arb(7)/2,(i,move_amp-ed)
        oscillation=move_amp-(lost_move+50*missing+arb(1)/1000)/(1-mass_error)
        assert oscillation>arb(585)/100,(i,oscillation)
    print(f'{method}: 289 mixed harmonics certify both midpoint-level/move ellipses with axes 1.17 and 2.45, the two-move square of radius 3.5 and an odd-spectrum margin above 5.85.')


def high():
    ctx.prec=144
    c,v=arb(8510),arb(8552)
    denom=(v*v-c*c-arb(1)/4).sqrt()
    assert 11632*8590<10**8 and 8590-38==8552
    assert denom>700 and v/(2*denom)<7 and 3*c/700<37 and c-7>8500
    assert 81*37*10**6*(-arb(41)).exp()<arb(1)/10000
    assert 81*48*(-arb(150)).exp()/(1-(-arb(200)).exp())<arb(1)/10000
    assert 14002+Fraction(2*(14000+84),16999)<14004
    assert 100+80*arb(2).log()<14004
    print('All unverified high zeros are suppressed throughout the strip. Every completion is before exp(14004) times H^(17001/16999).')


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    high()
    for method in ('sinc','cosine'):
        verify(method)
