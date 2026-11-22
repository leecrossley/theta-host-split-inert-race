// Segmented sieve on the two degree-one progressions. Logarithm enclosures use fixed-point integers.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr unsigned bits = 15;
constexpr unsigned bins = 1u << bits;
constexpr unsigned block = 1u << 18;

struct Bounds {
    i128 a = 0, b = 0, c = 0, d = 0;
};

std::string decimal(i128 x) {
    bool negative = x < 0;
    u128 v = negative ? u128(-x) : u128(x);
    std::string s;
    do { s.push_back(char('0' + v % 10)); v /= 10; } while (v);
    if (negative) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}


// Fixed-point Taylor bounds for log(1+z), 0 <= z <= 2^-15, in units of 2^-58.
std::pair<u64, u64> fineCorrection(u64 remainder, u64 base) {
    constexpr u64 unit = u64(1) << 58;
    u64 z = (u128(remainder) * unit) / base;
    auto down = [&](u64 a, u64 b) { return u64(u128(a)*b / unit); };
    auto up = [&](u64 a, u64 b) { return u64((u128(a)*b + unit-1) / unit); };
    u64 low2 = down(z,z), high2 = up(z+1,z+1);
    u64 low3 = down(low2,z), high3 = up(high2,z+1);
    u64 high4 = up(high3,z+1);
    i128 lower = i128(z) - (high2+1)/2 + low3/3 - (high4+3)/4;
    u64 upper = (z+1) - low2/2 + (high3+2)/3;
    return {u64(std::max<i128>(0,lower)), upper};
}

struct Sums {
    std::vector<u64> logTable, edges;
    std::vector<Bounds> sums;

    void add(std::size_t bucket, u64 p, u64 norm, int sign) {
        unsigned k = 63u - unsigned(__builtin_clzll(p));
        u64 j = ((p - (u64(1) << k)) << bits) >> k;
        u64 base = (bins + j) << k;
        u64 remainder = (p << bits) - base;
        auto range = fineCorrection(remainder, base);
        u64 lo = k*logTable[bins] + logTable[j] + range.first;
        u64 hi = k*(logTable[bins]+1) + logTable[j]+1 + range.second;
        constexpr unsigned shift = 38;
        i128 lowerRec = (u128(lo) << shift) / norm;
        i128 upperRec = ((u128(hi) << shift) + norm-1) / norm;
        auto &s = sums[bucket];
        if (sign > 0) {
            s.a += lo; s.b += hi; s.c += lowerRec; s.d += upperRec;
        } else {
            s.a -= hi; s.b -= lo; s.c -= upperRec; s.d -= lowerRec;
        }
    }
};

std::vector<int> legendre(unsigned p) {
    std::vector<int> a(p, -1);
    for (unsigned x=1; x<p; ++x) a[(x*x)%p] = 1;
    a[0] = 0;
    return a;
}

struct Strike {
    u64 p, first, second;
};

#ifndef THETA_REUSE_LOG_BOUNDS
int main(int argc, char** argv) {
    Sums acc;
    if (argc != 2) return 2;
    std::string host(argv[1]);
    if (host != "547" && host != "571") return 2;
    const unsigned u0 = host == "547" ? 547 : 571;
    const unsigned v0 = host == "547" ? 727 : 599;
    acc.logTable.resize(bins+1);
    for (auto &entry : acc.logTable) if (!(std::cin >> entry)) return 2;
    u64 x;
    while (std::cin >> x) acc.edges.push_back(x);
    if (acc.edges.empty() || !std::is_sorted(acc.edges.begin(), acc.edges.end())) return 2;
    const u64 cap = acc.edges.back();
    if (cap > (u64(1)<<38) || cap < 16) return 2;
    acc.sums.resize(acc.edges.size());
    unsigned root = 1;
    while (u64(root+1)*(root+1) <= cap) ++root;
    std::vector<bool> prime(root+1, true);
    prime[0] = prime[1] = false;
    for (unsigned p=2; u64(p)*p<=root; ++p) if (prime[p])
        for (unsigned j=p*p; j<=root; j+=p) prime[j] = false;
    auto t = legendre(u0), u = legendre(v0);
    std::vector<Strike> strikes;
    for (u64 p=3; p<=root; p+=2) if (prime[p]) {
        // First multiple p*m >= p^2 in each progression 16k+1, 16k+15.
        u64 inv = ((p%16)*(p%16)*(p%16))%16;
        auto start = [&](u64 r) {
            u64 m = p + ((r*inv + 16 - p%16)%16);
            return (p*m-r)/16;
        };
        strikes.push_back({p,start(1),start(15)});
        if (p==u0 || p==v0) continue;
        bool f1 = (p%16==1 || p%16==15) && t[p%u0]==u[p%v0];
        if (f1) continue;
        u64 sq = p*p;
        bool f2 = sq%16==1 || sq%16==15;
        if (!f2 && sq>cap/sq) continue;
        u64 norm = f2 ? sq : sq*sq;
        auto pos = std::lower_bound(acc.edges.begin(),acc.edges.end(),norm);
        if (pos != acc.edges.end()) acc.add(pos-acc.edges.begin(),p,norm,-1);
    }
    std::vector<unsigned char> a(block), b(block);
    std::size_t bucket = 0;
    for (u64 first=0; first<=cap/16; first+=block) {
        auto count = unsigned(std::min<u64>(block,cap/16-first+1));
        std::fill_n(a.begin(),count,1);
        std::fill_n(b.begin(),count,1);
        for (auto &s : strikes) {
            if (s.p*s.p>16*(first+count)+15) break;
            auto mark = [&](u64 &k, std::vector<unsigned char> &mask) {
                for (; k<first+count; k+=s.p) mask[k-first] = 0;
            };
            mark(s.first,a); mark(s.second,b);
        }
        for (unsigned j=0; j<count; ++j) {
            if (!(a[j] || b[j])) continue;
            auto enter = [&](bool possible, u64 p) {
                if (!possible || p<2 || p>cap || p==u0 || p==v0) return;
                int sign = t[p%u0];
                if (sign != u[p%v0]) return;
                while (acc.edges[bucket] < p) ++bucket;
                acc.add(bucket,p,p,-sign);
            };
            enter(a[j],16*(first+j)+1);
            enter(b[j],16*(first+j)+15);
        }
    }
    for (std::size_t j=0; j<acc.edges.size(); ++j) {
        const auto &s=acc.sums[j];
        std::cout << acc.edges[j] << ' ' << decimal(s.a) << ' ' << decimal(s.b)
                  << ' ' << decimal(s.c) << ' ' << decimal(s.d) << '\n';
    }
}
#endif
