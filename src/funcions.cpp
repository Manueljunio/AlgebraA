#include "../include/functions.h"

template <
    class result_t   = chrono::milliseconds,
    class clock_t    = chrono::steady_clock,
    class duration_t = chrono::milliseconds
>
auto since(chrono::time_point<clock_t, duration_t> const& start){
    return chrono::duration_cast<result_t>(clock_t::now() - start);
}

// Constructor
Functions::Functions(ll _n, ll _a):
    p(_n), a(_a), g(1) {}

// Fast exponentiation modulo m
// a = primo pequeno, b = d abaixo, m = possivel primo
ll Functions::fexp(ll a, ll b, ll m){
    ll ans = 1; a %= m; // a = a%m
    while(b > 0){
        if(b&1) ans = ans*a%m;
        a = a*a%m, b >>= 1;
    }
    return ans;
}

// Test if a number is composite
// n = possivel primo, a = primo pequeno
// d = (n-1)/2^s, s = menor expoente tal que d é ímpar
bool Functions::composite(ll n, ll a, ll d, ll s){
    ll x = fexp(a, d, n);
    if(x == 1 or x == n-1) return false;
    for(ll i = 1; i < s; i++){
        x = x*x%n;
        if(x == n-1) return false;
    }
    return true;
}

// Miller Rabin function that returns if a number is prime
bool Functions::millerRabin(ll pp){
    ll s = 0, d = pp-1;
    while(!(d&1)) d >>= 1, s++;

    bool found = true;
    for(auto x: primes){
        if(pp == x){ found = true; break; }
        if(composite(pp, x, d, s)) found = false;
    }
    return found;
}

// Returns the next prime after n
void Functions::NextPrime(){
    if(!(p&1)) p--;

    int cnt = 0;
    while(true){
        p += 2, cnt++;
        if(millerRabin(p)) break;
    }
    
    cout << "p: " << p << ", iterations_cnt: " << cnt << endl;
}

tuple<vector<ll>, vector<ll>, bool> Functions::primeFact(ll n, long long lim){
    vector<ll> fact, exp;
    for (long long i = 2; i*i <= n and (lim == -1 or i < lim); i++) if(n%i == 0) {
        fact.emplace_back(i), exp.emplace_back(0);
        while(n%i == 0) n /= i, exp[fact.size()-1]++;
    }
    bool partial = (n != 1 and !millerRabin(n));
    if(n > 1) fact.emplace_back(n), exp.emplace_back(1);

    return {fact, exp, partial};
}

// Returns a number that is a possible generator for Zp and its order interval
void Functions::Generator(){
    ll phi = p-1, n = phi;
    auto [fact, exp, partial] = primeFact(n, 1e7);

    vector<ll> ord(fact.size()); ll min_order = 1;
    for(size_t i = 0; i < fact.size(); i++){
        ll a = 2+rand()%(p-2);
        while(fexp(a, phi/fact[i], p) == 1) a = 2+rand()%(p-2);
        g *= fexp(a, phi/(fexp(fact[i], exp[i], p)), p), g %= p;
        if(!partial or i < fact.size()-1) min_order *= fexp(fact[i], exp[i], p);
    }
    if(partial) min_order *= 10000019; // first prime greater than 1e7
    
    cout << "g: " << g << endl;
    cout << "ord(g, p) is in: [" << min_order << ", " << phi << "]" << endl;
}

ll Functions::discLogBrute(ll g, ll a, ll p){
    for(int x = 0; x < p; x++) if(fexp(g, x, p) == a) return x;
    return -1;
}

ll Functions::discLogBabyGiantStep(ll g, ll a, ll p){
    a %= p;
    ll n = sqrt(p)+1, ans = -1;
    map<ll, ll> vals;
    for (ll pp = 1; pp <= n; pp++)
        vals[fexp(g, pp * n, p)] = pp;
    for (ll q = 0; q <= n; ++q) {
        ll cur = (fexp(g, q, p) * a) % p;
        if (vals.count(cur)) return vals[cur] * n - q;
    }

    return -1;
}

// Returns modular inverse of a modulo m, i.e., mod_inv * a = 1 mod m
ll Functions::mod_inv(ll a, ll m) {
    return a <= 1 ? a : m - (ll)(m/a) * mod_inv(m%a, m) % m;
}

pair<ll, ll> Functions::congPair(ll g, ll a, ll p, ll q, ll e, ll e1, ll e2){
    ll inv = mod_inv(e1, p), x = 0;
    for(int i = 1; i < e+1; i++){
        ll a = fexp(e1, fexp(q, e-1, p), p);
        ll b = fexp(e2*fexp(inv, x, p), fexp(q, e-i, p), p);
        x += discLogBabyGiantStep(a, b, p)*fexp(q, i-1, p);
    }

    return {x, fexp(q, e, p)};
}

// Returns the solution for a system of congruences x = a_i % m_i
ll Functions::chinese_remainder(vector<pair<ll, ll>> congruences){
    ll M = 1;
    for(auto const& c : congruences) M *= c.second;

    ll solution = 0;
    for(auto [a_i, m_i] : congruences)
        solution += a_i * (M/m_i)%M * mod_inv(M/m_i, m_i), solution %= M;

    return solution;
}

ll Functions::discLogPohligHellman(ll g, ll a, ll p){
    ll phi = p-1;

    auto [fact, exp, paraial] = primeFact(phi, -1);
    vector<pair<ll, ll>> cong;

    for(size_t i = 0; i < fact.size(); i++){
        ll e1 = fexp(g, phi/(fexp(fact[i], exp[i], p)), p);
        ll e2 = fexp(a, phi/(fexp(fact[i], exp[i], p)), p);
        cong.emplace_back(congPair(a, g, p, fact[i], exp[i], e1, e2));
    }

    return (chinese_remainder(cong)+p)%p;
}

// Returns the discrete logarithm of a modulo p in base g
void Functions::DiscreteLogarithm(){
    auto start = chrono::steady_clock::now();

    ll ans = -1;
    ans = discLogBrute(g, a, p);
    cout << "brute: " << ans << endl;
    ans = discLogBabyGiantStep(g, a, p);
    cout << "babygiantstep: " << ans << endl;
    ans = discLogPohligHellman(g, a, p);

    cout << "Elapsed(ms)=" << since(start).count() << std::endl;
    cout << "Discrete Logarithm of 'a' module 'p' in base 'g' = " << ans << endl;
}

// Return a factor of n
ll Functions::PollardRho(ll n) {
    if(!(n&1)) return 2;
    auto start = chrono::high_resolution_clock::now();

    ll x = 2, y = 2, c = rand()%10 - 1;
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
