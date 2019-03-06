#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <vector>
#include <tuple>
#include <set>

// Outputs a pair of whatever to a stream.
template<typename T, typename U>
std::ostream& operator<<(std::ostream& stream, const std::pair<T, U>& p) {
    std::cout << "(" << p.first << ", " << p.second << ")";
}

// Outputs a vector to a stream.
template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& container) {
    if (container.size() == 0) {
        return stream << "{}";
    } else {
        auto it = container.begin();
        stream << '{' << *it++;
        while (it != container.end()) {
            stream << ", " << *it++;
        }
        return stream << '}';
    }
}

// Outputs a set to a stream.
template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::set<T>& container) {
    if (container.size() == 0) {
        return stream << "{}";
    } else {
        auto it = container.begin();
        stream << '{' << *it++;
        while (it != container.end()) {
            stream << ", " << *it++;
        }
        return stream << '}';
    }
}

// Determines if n is a power of 2
bool is_power_of_two(int64_t n) {
    return (n & (n - 1)) == 0;
}

// Computes a^b in log(b) time.
// Assumes b is non-negative.
int64_t ipow(int64_t a, int64_t b) {
    if (b == 0) return 1;
    else if (b & 1) return a*ipow(a, b - 1);
    else return ipow(a*a, b >> 1);
}

// Computes GCD(a, b)
int64_t gcd(int64_t a, int64_t b) {
    while (b) {
        int64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Computes a*b (mod m) safely.
// This will be needed if you want to do
// computations with moduli greater than 2^30 - 1.
int64_t mod_mul(int64_t a, int64_t b, int64_t m) {
    int64_t result = 0;
    a %= m;
    b %= m;
    while (b) {
        if (b & 1) {
            result = (result + a) % m;
        }
        a = (a + a) % m;
        b >>= 1;
    }
    return result;
}

// Computes a^b (mod m)
int64_t mod_pow(int64_t a, int64_t b, int64_t m) {
    int64_t result = 1;
    a %= m;
    while (b) {
        if (b & 1) {
            result = (result * a) % m;//mod_mul(result, a, m);
        }
        a = (a*a) % m;//mod_mul(a, a, m);
        b >>= 1;
    }
    return result;
}

// Extended GCD Algorithm
std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b) {
    if (a == 0) {
        return std::make_tuple(b, 0ll, 1ll);
    } else {
        auto t = extended_gcd(b % a, a);
        auto g = std::get<0>(t);
        auto y = std::get<1>(t);
        auto x = std::get<2>(t);
        return std::make_tuple(g, x - (b/a)*y, y);
    }
}

// Computes a^(-1) (mod n)
int64_t mod_inv(int64_t a, int64_t n) {
    auto result = extended_gcd(a, n);
    if (std::get<0>(result) == 1) {
        return (n + std::get<1>(result)) % n;
    } else {
        throw std::runtime_error("mod_inv args not coprime.");
    }
}

// Computes the Legendre Symbol
int64_t legendre_symbol(int64_t a, int64_t p) {
    return mod_pow(a, (p - 1)/2, p);
}

// Determines if a is a quadratic residue modulo p
int64_t is_quadratic_residue(int64_t a, int64_t p) {
    return legendre_symbol(a, p) == 1;
}

// Computes the square root of n modulo p
// using the Tonelli-Shanks algorithm.
std::vector<int64_t> tonelli_shanks(int64_t n, int64_t p) {
    if (!is_quadratic_residue(n, p)) {
        return {};
    }

    int64_t q = p - 1;
    int64_t s = 0;
    while (~q & 1) {
    	q >>= 1;
    	s += 1;
    }

    // p = 3 (mod 4)
    // Hence, the solutions are trivial.
    if (s == 1) {
    	auto x = mod_pow(n, (p + 1)/4, p);
    	return {x, p - x};
    }

    // Select a quadratic non-residue (mod p)
    // This runs in expected logarithmic time
    // given Lagrange's theorem on the number of
    // quadratic residues modulo p.
    int64_t z = 0;
    for (int64_t k = 1; k < p; ++k) {
        if (!is_quadratic_residue(k, p)) {
            z = k;
            break;
        }
    }

    int64_t c = mod_pow(z, q, p);
    int64_t r = mod_pow(n, (q + 1)/2, p);
    int64_t t = mod_pow(n, q, p);
    int64_t m = s;

    while (t != 1) {
        int i = 1;
        int64_t x = (t*t) % p;
        while (x != 1) {
            x = (x*x) % p;
            i += 1;
        }
        int64_t b = mod_pow(c, (1ll << (m - i - 1)), p);
        // You could use mod_mul to ensure safety when
        // handling very large numbers.
        r = (r*b) % p;
        c = (b*b) % p;
        t = (t*c) % p;
        m = i;
    }

    return {r, p - r};
}

// Cartesian Product helper function
template<typename T>
void cartesian_product_helper(
    std::vector<std::vector<T>> const& v,
    std::vector<std::vector<T>>& result,
    std::vector<T>& path,
    int i)
{
    if (i == v.size()) {
        result.push_back(path);
    } else {
        for (int j = 0; j < v[i].size(); ++j) {
            path.push_back(v[i][j]);
            cartesian_product_helper(v, result, path, i + 1);
            path.pop_back();
        }
    }
}

// Generate the cartesian product of a bunch of vectors.
template<typename T>
std::vector<std::vector<T>> cartesian_product(std::vector<std::vector<T>> const& v) {
    std::vector<std::vector<T>> result;
    std::vector<T> path;
    result.reserve(50);
    cartesian_product_helper(v, result, path, 0);
    return result;
}

// Chinese Remainder Theorem on a sequence (ai, pi).
// This will compute x such that x = ai (mod pi) for all i
int64_t chinese_remainder_theorem(const std::vector<std::pair<int64_t, int64_t>>& pr) {
    int64_t x = 0;
    int64_t m = 1;

    for (int i = 0; i < pr.size(); ++i) {
        m *= pr[i].second;
    }

    for (int i = 0; i < pr.size(); ++i) {
        int64_t a = pr[i].first;
        int64_t pk = pr[i].second;

        int64_t y0 = m / pk;
        int64_t y1 = mod_inv(y0, pk);

        //x += mod_mul(mod_mul(a, y0, m), y1, m);
        x += ((a * y0 % m) * y1) % m;
        if (x >= m) x -= m;
    }

    return x;
}

// Determines if an integer is prime using
// the Miller-Rabin primality test. This is
// the deterministic variant constructed
// for integers that can fit in 64 bits.
bool is_prime(int64_t n) {
    static std::vector<int64_t> witnesses = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };

    if (n <= 2) return n == 2;
    if (~n & 1) return false;

    int64_t s = 0;
    int64_t d = n - 1;
    while (~d & 1) {
        d >>= 1;
        s += 1;
    }

    for (auto a : witnesses) {
        if (mod_pow(a, d, n) != 1) {
            bool ok = true;
            for (int64_t r = 0; r < s; ++r) {
                if (mod_pow(a, d * (1ll << r), n) == n - 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return false;
        }
    }

    return true;
}

// Finds the smallest prime factor of n
int64_t smallest_prime_factor(int64_t n) {
    static const int64_t sieve_limit = 20000000;
    static std::vector<int64_t> smf;
    static std::vector<int64_t> primes;

    // Perform the sieve on the first
    // call to this function.
    if (smf.size() == 0) {
        smf.resize(sieve_limit + 1);
        std::fill(smf.begin(), smf.end(), sieve_limit);
        std::vector<bool> composite(sieve_limit + 1, false);
        composite[0] = true;
        composite[1] = true;
        smf[0] = 1;
        smf[1] = 1;
        for (int64_t i = 2; i <= sieve_limit; ++i) {
            if (!composite[i]) {
                primes.push_back(i);
                smf[i] = i;
                for (int64_t j = 2*i; j <= sieve_limit; j += i) {
                    composite[j] = true;
                    smf[j] = std::min(smf[j], i);
                }
            }
        }
    }

    if (n <= sieve_limit) {
        return smf[n];
    } else {
        // Try to quickly terminate
        // in case n happens to be prime.
        if (is_prime(n)) {
            return n;
        }

        // Try small primes.
        for (auto p : primes) {
            if (n % p == 0) {
                return p;
            }
        }

        // In this case, n is a composite number divisible
        // by some prime greater than the sieve limit. It
        // will have one prime factor less than or equal
        // to its square root, and we'll find it by slow
        // trial division. This can be replaced by a careful
        // implementation of Pollard's Rho algorithm if needed.
        int p = sieve_limit + (sieve_limit & 1? 0 : 1);
        while (n % p) p += 2;
        return p;
    }
}

// Generates the prime factorization of n
std::vector<std::pair<int64_t, int64_t>> factorize(int64_t n) {
    std::vector<std::pair<int64_t, int64_t>> result;
    result.reserve(9);
    while (n > 1) {
        int64_t p = smallest_prime_factor(n);
        int64_t e = 0;
        while (n % p == 0) {
            n /= p;
            e += 1;
        }
        result.push_back(std::make_pair(p, e));
    }
    return result;
}

namespace std {
    template<>
    class hash<std::pair<int64_t, int64_t>> {
    public:
        std::size_t operator()(const std::pair<int64_t, int64_t>& t) const {
            static auto h = std::hash<int64_t>();
            return h(t.first) ^ h(t.second);
        }
    };
}

// Computes the square roots of a modulo n.
// This will also memoize results.
std::set<int64_t> sqrt_mod_n(int64_t a, int64_t n) {
    static std::unordered_map<std::pair<int64_t, int64_t>, std::set<int64_t>> memo;

    auto pr = std::make_pair(a, n);
    if (memo[pr].size() > 0) {
        return memo[pr];
    }

    if (n == 1) {
        return memo[pr] = {0ll};
    } else if (n > 1) {
        a %= n;
        if (gcd(a, n) == 1) {
            if (is_power_of_two(n)) {
                if (a % std::min(n, 8ll) == 1) {
                    int64_t k = 0;
                    int64_t t = n;
                    while (t > 1) {
                        t >>= 1;
                        k += 1;
                    }
                    if (k == 1) {
                        return memo[pr] = {1ll};
                    } else if (k == 2) {
                        return memo[pr] = {1ll, 3ll};
                    } else if (k == 3) {
                        return memo[pr] = {1ll, 3ll, 5ll, 7ll};
                    } else {
                        // Small optimization for the case of a == 1.
                        if (a == 1) {
                            return memo[pr] = {1ll, (n >> 1) - 1ll, (n >> 1) + 1ll, n - 1ll};
                        } else {
                            std::set<int64_t> roots;
                            for (auto x : sqrt_mod_n(a, n >> 1)) {
                                int64_t i = (((x*x - a) >> (k - 1)) % 2 == 1? 1 : 0);
                                int64_t r = x + i*(1 << (k - 2));
                                roots.insert(r);
                                roots.insert(n - r);
                            }
                            return memo[pr] = roots;
                        }
                    }
                }
            } else if (is_prime(n)) {
                std::set<int64_t> roots;
                for (auto r : tonelli_shanks(a, n)) {
                    roots.insert(r);
                }
                return memo[pr] = roots;
            } else {
                auto pe = factorize(n);

                // In the case of n being just an odd prime power.
                if (pe.size() == 1) {
                    int64_t p = pe[0].first;
                    int64_t k = pe[0].second;

                    // Since n = p^k, we have to solve
                    // the equation x^2 = a (mod p), then
                    // use Hensel's Lifting Lemma.
                    auto roots = tonelli_shanks(a, p);
                    int64_t pk = p;
                    int64_t pi = p*p;
                    for (int i = 2; i <= k; ++i) {
                        int64_t x = roots[0];
                        int64_t y = mod_inv(2, pk) * mod_inv(x, pk) % pk;
                        roots[0] = (pi + x - ((((x * x % pi) - a + pi)*y) % pi)) % pi;
                        roots[1] = pi - roots[0];
                        pk *= p;
                        pi *= p;
                    }
                    return memo[pr] = {roots[0], roots[1]};
                } else {
                    // Construct solutions for prime powers.
                    std::vector<std::vector<std::pair<int64_t, int64_t>>> solutions(pe.size());
                    for (int i = 0; i < pe.size(); ++i) {
                        auto m = ipow(pe[i].first, pe[i].second);
                        auto r = sqrt_mod_n(a, m);
                        solutions[i].reserve(r.size());
                        for (auto&& r0 : r) {
                            solutions[i].push_back(std::make_pair(r0, m));
                        }
                    }

                    // Construct all the possible square roots using
                    // the Chinese Remainder Theorem.
                    auto cp = cartesian_product(solutions);
                    std::set<int64_t> roots;
                    for (auto&& p : cp) {
                        roots.insert(chinese_remainder_theorem(p));
                    }
                    return memo[pr] = roots;
                }
            }
        }
    }
    // No solutions.
    return {};
}

int main() {
    // Example usage.
    std::cout << sqrt_mod_n(1, 5777) << '\n';
    std::cout << sqrt_mod_n(1, 19937) << '\n';
    std::cout << sqrt_mod_n(1, 9001) << '\n';
    return 0;
}
