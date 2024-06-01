#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>

using namespace std;
using namespace boost::multiprecision;

typedef cpp_int ll;

class Functions {
    private:
        vector<ll> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};
        ll p, a, g;

    public:
        Functions(ll _n, ll _a);
        ll fexp(ll a, ll b, ll m);
        bool composite(ll n, ll a, ll d, ll s);
        bool millerRabin(ll pp);
        void NextPrime();
        tuple<vector<ll>, vector<ll>, bool> primeFact(ll n, long long lim);
        void Generator();
        ll discLogBrute(ll g, ll a, ll p);
        ll discLogBabyGiantStep(ll g, ll a, ll p);
        ll mod_inv(ll a, ll m);
        pair<ll, ll> congPair(ll g, ll a, ll p, ll q, ll e, ll e1, ll e2);
        ll chinese_remainder(vector<pair<ll, ll>> congruences);
        ll discLogPohligHellman(ll g, ll a, ll p);
        void DiscreteLogarithm();
        ll PollardRho(ll n);
        tuple<vector<ll>, vector<ll>, bool> factorization(ll n);
};


#endif // FUNCTIONS_H
