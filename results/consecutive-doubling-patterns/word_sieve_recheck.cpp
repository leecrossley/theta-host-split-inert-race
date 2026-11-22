// Independent odd-integer sieve and positive atanh series for both prime-race hosts.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using U = std::uint64_t;
using Wide = __uint128_t;
using Signed = __int128_t;

constexpr U unit = U(1) << 62;
constexpr U subdivisions = U(1) << 12;

struct Bin {
    Signed low = 0, high = 0, reciprocalLow = 0, reciprocalHigh = 0;
};

std::string number(Signed value) {
    if (value == 0) return "0";
    Wide t = value < 0 ? Wide(-value) : Wide(value);
    std::string out;
    while (t) {
        out.push_back(char('0' + t % 10));
        t /= 10;
    }
    if (value < 0) out.push_back('-');
    std::reverse(out.begin(),out.end());
    return out;
}

U lowerProduct(U a, U b) { return U(Wide(a)*b/unit); }
U upperProduct(U a, U b) { return U((Wide(a)*b+unit-1)/unit); }

std::array<Wide,2> logBounds(U p, const std::vector<U>& nodes) {
    unsigned power = 63u - __builtin_clzll(p);
    U scale = U(1) << power;
    U j = (p-scale)*subdivisions/scale;
    U actual = p*subdivisions, base = (subdivisions+j)*scale;
    U h = U(Wide(actual-base)*unit/(actual+base));
    U hUp = h+1;
    U h2 = lowerProduct(h,h), h2Up = upperProduct(hUp,hUp);
    U h3 = lowerProduct(h,h2), h3Up = upperProduct(hUp,h2Up);
    U h5 = lowerProduct(h3,h2), h5Up = upperProduct(h3Up,h2Up);
    // t < 1/8192. Terms beyond the fifth power add less than one unit.
    Wide low = Wide(power)*nodes.back()+nodes[j]+2*(Wide(h)+h3/3+h5/5);
    Wide high = Wide(power)*(nodes.back()+1)+nodes[j]+1+2*(Wide(hUp)+(h3Up+2)/3+(h5Up+4)/5)+1;
    return {low,high};
}

std::vector<int> character(unsigned q) {
    std::vector<int> values(q);
    for (unsigned n=1;n<q;++n) {
        unsigned exponent=(q-1)/2, base=n, result=1;
        while (exponent) {
            if (exponent%2) result=result*base%q;
            base=base*base%q;
            exponent/=2;
        }
        values[n]=result==1 ? 1 : -1;
    }
    return values;
}

#ifndef THETA_REUSE_LOG_BOUNDS
int main(int argc,char** argv) {
    if (argc!=2) return 2;
    std::string host(argv[1]);
    if (host!="547" && host!="571") return 2;
    unsigned u=host=="547" ? 547 : 571, v=host=="547" ? 727 : 599;
    auto chi=character(u), psi=character(v);
    std::vector<U> nodes(subdivisions+1),edges;
    for (U& n:nodes) if (!(std::cin>>n)) return 2;
    U t;
    while (std::cin>>t) edges.push_back(t);
    if (edges.empty() || !std::is_sorted(edges.begin(),edges.end())) return 2;
    U cap=edges.back();
    if (cap<16 || cap>(U(1)<<38)) return 2;
    std::vector<Bin> sums(edges.size());
    std::vector<U> primes;
    for (U n=3;n*n<=cap;n+=2) {
        bool ok=true;
        for (U p:primes) {
            if (p*p>n) break;
            if (n%p==0) { ok=false; break; }
        }
        if (ok) primes.push_back(n);
    }
    constexpr U count=U(1)<<19;
    for (U lo=3;lo<=cap;lo+=2*count) {
        U size=std::min(count,(cap-lo)/2+1), last=lo+2*(size-1);
        std::vector<char> flags(size,1);
        for (U p:primes) {
            if (p*p>last) break;
            U first=std::max(p*p,((lo+p-1)/p)*p);
            if (first%2==0) first+=p;
            for (U j=(first-lo)/2;j<size;j+=p) flags[j]=0;
        }
        for (U k=0;k<size;++k) if (flags[k]) {
            U p=lo+2*k;
            if (p==u || p==v) continue;
            U rem=p%16, square=rem*rem%16;
            int fieldDegree=(rem==1 || rem==15) ? 1 : (square==1 || square==15) ? 2 : 4;
            if (chi[p%u]*psi[p%v]<0 && fieldDegree==1) fieldDegree=2;
            Wide norm=p;
            for (int d=1; d<fieldDegree && norm<=cap; ++d) norm*=p;
            if (norm>cap) continue;
            U n=U(norm);
            auto bounds=logBounds(p,nodes);
            Signed l=Signed(bounds[0]), h=Signed(bounds[1]);
            Signed ll=Signed((bounds[0]<<40)/n), hh=Signed(((bounds[1]<<40)+n-1)/n);
            Bin& s=sums[std::lower_bound(edges.begin(),edges.end(),n)-edges.begin()];
            if (fieldDegree==1 && chi[p%u]<0) {
                s.low+=l; s.high+=h; s.reciprocalLow+=ll; s.reciprocalHigh+=hh;
            } else {
                s.low-=h; s.high-=l; s.reciprocalLow-=hh; s.reciprocalHigh-=ll;
            }
        }
    }
    for (std::size_t j=0;j<edges.size();++j) {
        const auto& s=sums[j];
        std::cout<<edges[j]<<' '<<number(s.low)<<' '<<number(s.high)<<' '
                 <<number(s.reciprocalLow)<<' '<<number(s.reciprocalHigh)<<'\n';
    }
}
#endif
