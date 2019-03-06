# Square Root Modulo N

This is a C++ implementation of a generalized Square Root function. I wrote it while working on some Project Euler problems.

# Algorithm

To compute the square root of a non-zero `a` modulo an odd prime `p`, we use the [Tonelli-Shanks Algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm).
This square root exists if and only if `a` is quadratic residue modulo `p`. To test this, 
we use the [Legendre Symbol](https://en.wikipedia.org/wiki/Legendre_symbol) in conjunction with [Euler's Criterion](https://en.wikipedia.org/wiki/Euler%27s_criterion).
In the case of square roots modulo a composite number `n`, we first factor `n` into a product of prime powers.
The next step is to find the square roots modulo each of those prime powers. To do this, we can use Tonelli-Shanks 
to compute the square roots modulo the primes, and then [Hensel's Lifting Lemma](https://en.wikipedia.org/wiki/Hensel%27s_lemma) 
to raise the solutions to the appropriate prime power moduli. Finally, we have to construct the solution using the 
[Chinese Remainder Theorem](https://en.wikipedia.org/wiki/Chinese_remainder_theorem). Of course, square roots modulo odd primes exist in pairs, 
so when constructing square roots, we have to consider the cartesian product of the solutions for the prime power moduli.

For the case of powers of 2, we handle them specially. If our modulus is `2`, the solution is `1` if and only if `a` is 1.
In the case of `4`, the only quadratic residue is `1`, and the solution is `{1, 3}`. For the higher powers of `2`, there are 
always exactly four unique solutions if and only if `a = 1 (mod 8)`. Otherwise, there are no solutions. These solutions can be constructed 
inductively starting from the base case of `8` corresponding to the solutions `{1, 3, 5, 7}`. Note that this is because `1` is the only 
quadratic residue in `(Z/8Z)*`.
