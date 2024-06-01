#include "../include/functions.h"

// Return a factor of n
ll Functions::PollardRho(ll n) {
    if(!(n&1)) return 2;
    auto start = chrono::high_resolution_clock::now();

    ll x = 2, y = 2, c = rand()%1000 - 1;
    ll d = 1;
    while(d==1) {
        x = (x*x+c)%n;
        y = (y*y+c)%n; y = (y*y+c)%n;
        d = x>=y? x-y : y-x;
        d = gcd(n,d);
        auto now = chrono::high_resolution_clock::now();
        
        if(chrono::duration_cast<chrono::seconds>(now-start).count()>20) return -1;
    }
    return d;
}

// Returns the factorization of n
tuple<vector<ll>, vector<ll>, bool> Functions::factorization(ll n) {
    vector<ll> fact, exp;
    bool p = false;
    while(n>1) {
        if(millerRabin(n)) {
            fact.emplace_back(n);
            exp.emplace_back(1);
            break;
        }
        ll factor = PollardRho(n);
        while(factor==n){
            factor = PollardRho(n);
        }
        if(factor==-1) {
            fact.emplace_back(n);
            exp.emplace_back(1);
            p = true;
            break;
        }
        while(n%factor==0) {
            n/=factor;
            if(fact.empty() || fact.back()!=factor) {
                fact.emplace_back(factor);
                exp.emplace_back(1);
            } else exp.back()++;
        }
    }

    return {fact,exp,p};
}