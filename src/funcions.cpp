#include "../include/functions.h"
#include <chrono>
#include <stdexcept>

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
    prime(_n), a(_a), g(1) {}

// Fast exponentiation modulo m
ll Functions::fexp(ll x, ll b, ll m){
    ll ans = 1; x %= m;
    while(b > 0){
        if(b&1) ans = ans*x%m;
        x = x*x%m, b >>= 1;
    }
    return ans;
}

// Test if a number is composite
bool Functions::composite(ll n, ll r, ll d, ll s){
    ll x = fexp(r, d, n);
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
    if(!(prime&1)) prime--;

    int cnt = 0;
    while(true){
        prime += 2, cnt++;
        if(millerRabin(prime)) break;
    }
    
    cout << "p: " << prime << ", iterations: " << cnt << endl;
}

// Returns the factorization of n using trivial algorithm
void Functions::primeFact(ll n, long long lim){
    fact.clear(), exp.clear(), partial = 0;
    for (long long i = 2; i*i <= n and (lim == -1 or i < lim); i++) if(n%i == 0) {
        fact.emplace_back(i), exp.emplace_back(0);
        while(n%i == 0) n /= i, exp[fact.size()-1]++;
    }
    partial = (n != 1 and !millerRabin(n));
    if(n > 1) fact.emplace_back(n), exp.emplace_back(1);
}

// Return a factor of n
ll Functions::PollardRho(ll n) {
    if(!(n&1)) return 2;

    auto s = chrono::high_resolution_clock::now();

    ll x = (rand()%n)+1, y = x, c = rand()%n+1;
    ll d = 1;
    while(d==1) {
        auto now = chrono::high_resolution_clock::now();
        if(chrono::duration_cast<chrono::seconds>(now-s).count()>10.0) return 0;

        x = (x*x+c+n)%n;
        y = (y*y+c+n)%n; y = (y*y+c+n)%n;
        d = x>=y? x-y : y-x;
        d = gcd(n,d);
    }
    return d;
}

// Factorize n and compute the facts array with all the divisors of n
void Functions::factorize(ll n){
    if(n == 1) return;

    if(millerRabin(n)){
        fact.emplace_back(n);
        return;
    }

    ll d = PollardRho(n);
    if(d==0) {
        partial = true;
        return;
    }    

    factorize(d), factorize(n/d);
}

// Returns the factorization of n using pollard rho algorithm
void Functions::primeFactRho(ll n){
    fact.clear(), exp.clear(), partial = 0;
    factorize(n);
    sort(fact.begin(), fact.end());
    fact.erase(unique(fact.begin(), fact.end()), fact.end());
    for(size_t i = 0; i < fact.size(); i++){
        exp.emplace_back(0);
        while(n%fact[i] == 0) exp[i]++, n /= fact[i];
    }
}

// Returns a number that is a possible generator for Zp and its order interval
void Functions::Generator(){
    ll phi = prime-1, n = phi;
    primeFactRho(n);

    vector<ll> ord(fact.size()); ll min_order = 1;
    for(size_t i = 0; i < fact.size(); i++){
        ll b = 2+rand()%(prime-2);
        while(fexp(b, phi/fact[i], prime) == 1) b = 2+rand()%(prime-2);
        g *= fexp(b, phi/(fexp(fact[i], exp[i], prime)), prime), g %= prime;
        if(!partial or i < fact.size()-1) min_order *= fexp(fact[i], exp[i], prime);
    }
    
    cout << "g: " << g << ", minimum order: " << min_order << endl;
}

// Brute the discrete logarithm testing if g^x % p == a
ll Functions::discLogBrute(ll g, ll h, ll p){
    h %= p; ll f = 1;
    for(int x = 0; x < p; x++){
        if(f == h) return x;
        f = f*g%p;
    }
    return -1; // Discrete log not found
}

ll Functions::discLogBabyGiantStep(ll g, ll h, ll p) {//
    auto s = chrono::high_resolution_clock::now();
    h %= p;
    ll n = sqrt(p) + 1;
    map<ll, ll> vals;
    for (ll i = 1; i <= n; i++) {
        auto now = chrono::high_resolution_clock::now();
        if(chrono::duration_cast<chrono::seconds>(now-s).count()>20) throw length_error("");
        vals[fexp(g, i * n, p)] = i;
    }
    for (ll j = 0; j <= n; j++) {
        ll cur = (fexp(g, j, p) * h) % p;
        if (vals.count(cur))
            return vals[cur] * n - j;
        auto now = chrono::high_resolution_clock::now();
        if(chrono::duration_cast<chrono::seconds>(now-s).count()>20) throw length_error("");
    }
    return -1;
}

// Returns modular inverse of a modulo m, i.e., mod_inv * a = 1 mod m
ll Functions::mod_inv(ll x, ll m) {
    ll m0 = m, t, q;
    ll x0 = 0, x1 = 1;
    if (m == 1) return 0;

    while (x > 1) {
        q = x / m;
        t = m;
        m = x % m;
        x = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

    if (x1 < 0) x1 += m0;

    return x1;
}

// Returns the solution for a system of congruences x = a_i % m_i
ll Functions::chinese_remainder(vector<pair<ll, ll>> congruences){
    ll M = 1;
    for(auto const& c : congruences) M *= c.second;

    ll solution = 0;
    for(auto [a_i, m_i] : congruences)
        solution += a_i * (M/m_i) * mod_inv(M/m_i, m_i), solution %= M;

    return solution;
}

pair<ll, ll> Functions::congPair(ll p, ll q, ll e, ll e1, ll e2) {
    ll inv = mod_inv(e1, p);
    ll x = 0;
    ll q_pow_i = 1;
    for (int i = 0; i < e; i++) {
        ll b = (e2 * fexp(inv, x, p)) % p;
        ll c = fexp(e1, q_pow_i, p);
        ll dlog = discLogBabyGiantStep(c, b, p);
        x = (x + dlog * q_pow_i) % (p - 1);
        q_pow_i = (q_pow_i * q) % (p - 1);
    }
    return {x, fexp(q, e, p)};
}

ll Functions::discLogPohligHellman(ll g, ll h, ll p) {
    g %= p, h %= p;

    ll phi = p - 1;
    if(partial) primeFactRho(phi);
    vector<pair<ll, ll>> cong;

    for (size_t i = 0; i < fact.size(); i++) {
        ll q = fact[i];
        ll e = exp[i];
        ll e1 = fexp(g, phi / (fexp(q, e, p)), p);
        ll e2 = fexp(h, phi / fexp(q, e, p), p);
        cong.emplace_back(congPair(p, q, e, e1, e2));
    }

    return chinese_remainder(cong);
}

// Returns the discrete logarithm of a modulo p in base g
void Functions::DiscreteLogarithm(){
    start = chrono::high_resolution_clock::now();

    ll ans = -1;
    bool b = false;
    try {
        if(prime < (ll)1e6) ans = discLogBrute(g, a, prime);
        else if(prime < (ll)1e12) ans = discLogBabyGiantStep(g, a, prime);
        else ans = discLogPohligHellman(g, a, prime);
        b = true;
    } catch (length_error &e) {cout << "i: time limit\n";}
    if(b) cout << "i: " << ans << endl;
    
    auto now = chrono::high_resolution_clock::now();
    cout << "t: " << chrono::duration_cast<chrono::seconds>(now-start).count() << "s" << std::endl;
}
