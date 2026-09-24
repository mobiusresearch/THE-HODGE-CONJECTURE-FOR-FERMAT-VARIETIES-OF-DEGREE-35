"""Exact checks for Section 8 of the paper
"The Hodge conjecture for Fermat varieties of degree 35" (the level m = 420).

Requirements: python3 (>= 3.8) and the file gap.py of this directory.  Everything is exact integer
arithmetic in pure Python; python-flint is optional and used only for an independent computation of
lattice indices (Hermite and Smith normal forms), as in gap.py.

    python3 level420.py          (about 1 minute)

Notation as in gap.py and Section 5 of the paper: R_m, K_m, B_m, D_m, S_m, u(a), len(x), and
rho_q^(m) = (m/q) rho_q.  Here m = 420.  For x = sum_y c_y (y) in R_420 and a divisor N > 1 of 420,
n_N(x) is the sum of the c_y over the y of order N; every y of order N is y = (420/N) u_y with a unique
u_y in (Z/N)^x.  For a Dirichlet character chi of conductor f dividing N put
    s_N(chi; x) = sum_{ord y = N} c_y chi(u_y).
Characters are evaluated exactly in Z[zeta_12] = Z[i, zeta_3] (basis 1, z, z^2, z^3, z = zeta_12,
z^4 = z^2 - 1).  The characters named in Section 8 are
    chi_4 (mod 4), chi_3 (mod 3), eps_5, eps_7 (the quadratic characters mod 5 and 7),
    lambda (mod 5, lambda(2) = i), theta (mod 7, theta(3) = zeta_3 = z^4),
and the order classes R = {2, 3, 5, 6, 7, 10, 14, 21, 30, 42, 70, 210}.
The two characters of dimension 14 used in Section 8 are
    w_1 = (56,69,176,191,193,199,204,209,210,236,243,253,257,267,296,301),
    w_2 = (1,17,28,132,141,168,191,199,236,253,257,267,296,371,387,416).

Part A  (Lemma 8.2, Proposition 8.3) Q(420) = {12,15,20,21,28,35,420}; the generators of S_420 lie in
        B_420 and nu_15 (the parity of the number of entries of order 15) vanishes on all of them and on
        the six rho_q^(420); w_1, w_2 are Hodge characters with nu_15 = 1;
        B_420 = S_420 + sum_q Z u(rho_q^(420)) + Z u(w_1), and 2 u(w_1) lies in S_420.  Hence L''_420 =
        S_420 + sum_q Z u(rho_q^(420)) = {x in B_420 : nu_15(x) = 0}, of index 2 in B_420.
Part B  (Lemma 8.4, (8.2)-(8.12)) for every odd character chi mod 420 and every y != 0 in Z/420, the sum
        sum_t conj(chi)(t) (2<ty> - 420) is given by the formula in the proof of Lemma 8.4; B_{1,chi} != 0;
        the eleven relations (8.2)-(8.12) have the stated coefficients and vanish on the generators of B_420.
Part C  (Lemmas 8.5-8.8, Theorem 8.1) the parity relations, the congruence for n'_420(x), the six-class
        congruence on unit entries, the norm computation for (8.11), and the count pattern of w_1, w_2;
        cross-check by enumeration: every multiset of 2 or 4 units satisfying the seven conditions
        s_420(chi) = 0 (chi odd primitive of conductor 420) has n'_420 even.
Part D  (cross-check) every Hodge character of level 420 with four entries, one of which divides 420,
        has nu_15 = 0.

Exit status 0 if and only if every check passes.
"""
import sys
import time
from itertools import combinations_with_replacement
from math import gcd

import gap as G

M = 420
checks = []


def check(name, cond):
    checks.append((name, bool(cond)))
    print(("PASS  " if cond else "FAIL  ") + name)
    return bool(cond)


# ====================================================================== Z[zeta_12]
def red(bins):
    """sum_k bins[k] z^k (k = 0..11) as a vector in the basis 1, z, z^2, z^3 (z^4 = z^2 - 1)"""
    v = list(bins) + [0] * (12 - len(bins))
    for k in range(11, 3, -1):          # z^k = z^(k-2) - z^(k-4)
        c = v[k]
        if c:
            v[k] = 0
            v[k - 2] += c
            v[k - 4] -= c
    return tuple(v[:4])


def zp(e):
    b = [0] * 12
    b[e % 12] = 1
    return red(b)


def zmul(a, b):
    bins = [0] * 12
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    bins[i + j] += x * y
    return red(bins)


def zadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def zscal(k, a):
    return tuple(k * x for x in a)


ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)
I = zp(3)            # i = z^3
ZETA3 = zp(4)        # zeta_3 = z^4


def zconj(a):
    """complex conjugation z -> z^-1 = z^11"""
    out = ZERO
    for k, x in enumerate(a):
        if x:
            out = zadd(out, zscal(x, zp(-k)))
    return out


# ====================================================================== characters mod 420
UNITS = G.units(M)
IND5 = {1: 0, 2: 1, 4: 2, 3: 3}                 # 2^k mod 5
IND7 = {1: 0, 3: 1, 2: 2, 6: 3, 4: 4, 5: 5}     # 3^k mod 7


def crt(u):
    return (0 if u % 4 == 1 else 1, 0 if u % 3 == 1 else 1, IND5[u % 5], IND7[u % 7])


def cexp(chi, u):
    """chi = (k4, k3, k5, k7): chi(u) = z^(6 k4 e4 + 6 k3 e3 + 3 k5 a5 + 2 k7 b7) for u a unit mod 420 with
    u = (-1)^e4 mod 4, (-1)^e3 mod 3, 2^a5 mod 5, 3^b7 mod 7"""
    e4, e3, a5, b7 = crt(u)
    k4, k3, k5, k7 = chi
    return (6 * k4 * e4 + 6 * k3 * e3 + 3 * k5 * a5 + 2 * k7 * b7) % 12


def conductor(chi):
    k4, k3, k5, k7 = chi
    return (4 if k4 % 2 else 1) * (3 if k3 % 2 else 1) * (5 if k5 % 4 else 1) * (7 if k7 % 6 else 1)


def conj(chi):
    k4, k3, k5, k7 = chi
    return (k4 % 2, k3 % 2, (-k5) % 4, (-k7) % 6)


def lift(r, f):
    """a unit mod 420 congruent to r mod f (r prime to f)"""
    x = r % f
    while gcd(x, M) != 1:
        x += f
    return x


def chi_at(chi, r):
    """chi(r) for r prime to the conductor f of chi, as an element of Z[zeta_12]; 0 if gcd(r, f) > 1"""
    f = conductor(chi)
    if gcd(r, f) != 1:
        return ZERO
    return zp(cexp(chi, lift(r, f)))


ALL = [(k4, k3, k5, k7) for k4 in range(2) for k3 in range(2) for k5 in range(4) for k7 in range(6)]
ODD = [c for c in ALL if cexp(c, M - 1) == 6]
CHI4, CHI3, EPS5, LAM, EPS7, THETA = (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 2, 0), (0, 0, 1, 0), (0, 0, 0, 3), (0, 0, 0, 2)


def mulc(*cs):
    return tuple(sum(c[j] for c in cs) % (2, 2, 4, 6)[j] for j in range(4))


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


DIVS = [d for d in G.divisors(M) if d > 1]


def u_of(y):
    """(N, u_y) for y != 0 in Z/420: y has order N and y = (420/N) u_y"""
    N = G.order(y, M)
    return N, (y // (M // N)) % N


def euler_coeff(chi, N):
    """(phi(420)/phi(N)) prod_{p | N, p not dividing f} (1 - conj(chi)(p)), or None if f does not divide N"""
    f = conductor(chi)
    if N % f:
        return None
    c = (phi(M) // phi(N), 0, 0, 0)
    for p in G.prime_factors(N):
        if f % p:
            c = zmul(c, zadd(ONE, zscal(-1, zconj(chi_at(chi, p)))))
    return c


def fB1(chi):
    """f B_{1,chi} = sum_{a = 1}^{f} chi(a) a, f the conductor (chi nontrivial)"""
    f = conductor(chi)
    out = ZERO
    for a in range(1, f + 1):
        out = zadd(out, zscal(a, chi_at(chi, a)))
    return out


def s_value(chi, x, N):
    """s_N(chi; x)"""
    out = ZERO
    for y, c in x.items():
        NN, uy = u_of(y)
        if NN == N:
            out = zadd(out, zscal(c, chi_at(chi, uy)))
    return out


def relation_value(chi, x, coeffs):
    """sum_N coeffs[N] s_N(chi; x)"""
    out = ZERO
    for N, a in coeffs.items():
        out = zadd(out, zmul(a, s_value(chi, x, N)))
    return out


def n_class(x, orders):
    return sum(c for y, c in x.items() if G.order(y, M) in orders)


def zfmt(a):
    """an element of Z[zeta_12] = Z[i, zeta_3], written in the basis 1, i, zeta_3, i zeta_3
    (z = zeta_12 = -i zeta_3, z^2 = 1 + zeta_3, z^3 = i)"""
    c1, cz, cz2, cz3 = a
    terms = []
    for coef, name in ((c1 + cz2, ""), (cz3, "i"), (cz2, "zeta_3"), (-cz, "i zeta_3")):
        if coef:
            terms.append((str(coef) if (coef not in (1, -1) or not name) else ("-" if coef < 0 else "")) + name)
    return " + ".join(terms).replace("+ -", "- ") if terms else "0"


def fmt(a):
    return "(" + ",".join(map(str, a)) + ")"


W1 = (56, 69, 176, 191, 193, 199, 204, 209, 210, 236, 243, 253, 257, 267, 296, 301)
W2 = (1, 17, 28, 132, 141, 168, 191, 199, 236, 253, 257, 267, 296, 371, 387, 416)
RSET = {2, 3, 5, 6, 7, 10, 14, 21, 30, 42, 70, 210}


# ====================================================================== Part A
def part_A():
    print("Part A. The class of xi_420 (Lemma 8.2, Proposition 8.3)")
    Q = [d for d in G.divisors(M) if d > 1 and G.mu_prime(d) == 1]
    check("Q(420) = %s: 420 has the four divisors %s in P, and |Q(420)| = 2^3 - 1 = 7" % (Q, G.P_divisors(M)),
          Q == [12, 15, 20, 21, 28, 35, 420] and G.P_divisors(M) == [3, 4, 5, 7])
    Sg = G.S_generators(M)
    reps = {q: G.inflate(G.REPS[q], q, M) for q in Q if q != 420}
    check("the %d generators of S_420 lie in B_420, and the six rho_q^(420) are Hodge characters: %s"
          % (len(Sg), "; ".join("rho_%d^(420) = %s" % (q, fmt(a)) for q, a in reps.items())),
          all(G.in_B(g, M) for g in Sg) and all(G.is_hodge(a, M) for a in reps.values()))
    check("15 divides 420, gcd(15, 420/15) = 1, and nu_15 vanishes on all %d generators of S_420 (Lemma 8.2)" % len(Sg),
          gcd(15, 28) == 1 and all(G.nu(g, 15, M) == 0 for g in Sg))
    ent = {q: [y for y in a if G.order(y, M) == 15] for q, a in reps.items()}
    check("nu_15 vanishes on the six rho_q^(420); their entries of order 15 are %s" % ent,
          all(len(e) % 2 == 0 for e in ent.values()) and ent[15] == [28, 364])
    for nm, w in (("w_1", W1), ("w_2", W2)):
        e15 = [y for y in w if G.order(y, M) == 15]
        check("%s = %s is a Hodge character of level 420 with %d entries (dimension %d), with exactly one entry of "
              "order 15, namely %s" % (nm, fmt(w), len(w), len(w) - 2, e15),
              G.is_hodge(w, M) and len(w) == 16 and len(e15) == 1)
    gens = Sg + [G.u(a, M) for a in reps.values()] + [G.u(W1, M)]
    ok, r, E, D, primes = G.prove_equal_to_B(M, gens)
    check("B_420 = S_420 + sum_q Z u(rho_q^(420)) + Z u(w_1): the echelon basis of L+ has %d rows, pivot product %d%s"
          % (r, D, (", full rank modulo %s" % " and ".join(map(str, primes))) if primes else ""),
          ok and r == (M - 1) - len(G.half_units(M)))
    ES = G.Echelon()
    for g in Sg:
        ES.insert(g)
    Lpp = ES.copy()
    for a in reps.values():
        Lpp.insert(G.u(a, M))
    check("2 u(w_1) lies in S_420 and u(w_1) does not lie in L''_420 = S_420 + sum_q Z u(rho_q^(420)); so "
          "[B_420 : L''_420] = 2 and L''_420 = {x in B_420 : nu_15(x) = 0}",
          ES.contains(G.scal(2, G.u(W1, M))) and not Lpp.contains(G.u(W1, M)))
    if G.flint is not None:
        base = Sg + [G.u(a, M) for a in reps.values()]
        res = G.flint_indices(M, [Sg, base, base + [G.u(W1, M)], base + [G.u(W1, M), {M // 2: 1}]])
        check("FLINT: [K_420 : L] = %s for L = S_420, L''_420, L''_420 + Z u(w_1), and the latter + Z (210) "
              "(expected [256, 4, 2, 1]); Hermite pivot products %s agree with the pure-Python echelon bases"
              % ([x[1] for x in res], [x[0] for x in res]),
              [x[1] for x in res] == [256, 4, 2, 1]
              and [x[0] for x in res] == [ES.pivot_product(), Lpp.pivot_product(), E.pivot_product(), D])
    return gens


# ====================================================================== Part B
# The relations (8.2)-(8.12) as printed in Section 8: character, normalising factor kappa, and coefficients.
# Each asserts sum_N coeffs[N] s_N(chi; x) = 0 for x in B_420; the coefficients computed from Lemma 8.4 must
# equal kappa * coeffs.  Coefficients are elements of Z[zeta_12]; integers are written as (k, 0, 0, 0).
def Zi(a, b=0):
    return zadd((a, 0, 0, 0), zscal(b, I))


RELATIONS = [
    ("(8.2)", "chi_3", CHI3, 24, {3: Zi(2), 6: Zi(4), 12: Zi(2), 15: Zi(1), 30: Zi(2), 60: Zi(1)}),
    ("(8.3)", "chi_3 eps_5", mulc(CHI3, EPS5), 4, {15: Zi(3), 105: Zi(1)}),
    ("(8.4)", "eps_7", EPS7, 8, {7: Zi(2), 21: Zi(2), 35: Zi(1), 105: Zi(1)}),
    ("(8.5)", "eps_5 eps_7", mulc(EPS5, EPS7), 4, {35: Zi(1), 70: Zi(2), 140: Zi(1)}),
    ("(8.6)", "chi_4", CHI4, 16, {4: Zi(3), 12: Zi(3), 28: Zi(1), 84: Zi(1)}),
    ("(8.7)", "chi_4 eps_5", mulc(CHI4, EPS5), 12, {20: Zi(1)}),
    ("(8.8)", "chi_4 chi_3 eps_7", mulc(CHI4, CHI3, EPS7), 4, {84: Zi(1)}),
    ("(8.9)", "chi_4 chi_3 eps_5 eps_7", mulc(CHI4, CHI3, EPS5, EPS7), 1, {420: Zi(1)}),
    ("(8.10)", "lambda", LAM, 2,
     {5: Zi(12), 10: Zi(12, 12), 15: Zi(6, -6), 20: Zi(6, 6), 30: Zi(12), 35: Zi(2, 2), 60: Zi(6),
      70: Zi(0, 4), 105: Zi(2), 140: Zi(0, 2), 210: Zi(2, 2), 420: Zi(1, 1)}),
    ("(8.11)", "chi_4 eps_5 theta", mulc(CHI4, EPS5, THETA), 1,
     {140: Zi(2), 420: zadd(ONE, zscal(-1, zconj(ZETA3)))}),
    ("(8.12)", "chi_4 chi_3 lambda theta", mulc(CHI4, CHI3, LAM, THETA), 1, {420: Zi(1)}),
]


def part_B(gens):
    print()
    print("Part B. The character relations (Lemma 8.4 and (8.2)-(8.12))")
    t0 = time.time()
    bad = 0
    nz = True
    for chi in ODD:
        f = conductor(chi)
        cb = conj(chi)
        fb = fB1(cb)                                   # f B_{1, conj chi}
        nz &= fB1(chi) != ZERO
        ce = [(-cexp(chi, t)) % 12 for t in UNITS]     # exponents of conj(chi)(t)
        for y in range(1, M):
            bins = [0] * 12
            for e, t in zip(ce, UNITS):
                bins[e] += 2 * ((t * y) % M) - M
            lhs = zscal(f, red(bins))
            N, uy = u_of(y)
            a = euler_coeff(chi, N)
            rhs = ZERO if a is None else zmul(zmul(zscal(2 * M, fb), a), chi_at(chi, uy))
            bad += lhs != rhs
    check("Lemma 8.4: for all %d odd characters chi mod 420 and all 419 elements y != 0 of Z/420, "
          "sum_t conj(chi)(t) (2<ty> - 420) = 2*420 B_{1,conj chi} (phi(420)/phi(N)) prod_{p | N, p not dividing f} "
          "(1 - conj(chi)(p)) chi(u_y) if f | N, and = 0 otherwise (N the order of y)" % len(ODD), bad == 0)
    check("B_{1,chi} != 0 for all %d odd characters chi mod 420" % len(ODD), nz)
    for tag, name, chi, kappa, coeffs in RELATIONS:
        f = conductor(chi)
        comp = {N: euler_coeff(chi, N) for N in DIVS if euler_coeff(chi, N) not in (None, ZERO)}
        same = set(comp) == set(coeffs) and all(comp[N] == zscal(kappa, coeffs[N]) for N in coeffs)
        odd = cexp(chi, M - 1) == 6
        van = all(relation_value(chi, g, coeffs) == ZERO for g in gens)
        check("%s: chi = %s is odd of conductor %d; Lemma 8.4 gives %d times the coefficients %s; the relation "
              "vanishes on all %d generators of B_420" % (tag, name, f, kappa,
                                                         "{" + ", ".join("%d: %s" % (N, zfmt(coeffs[N])) for N in sorted(coeffs)) + "}", len(gens)),
              odd and same and van)
    print("  (%.1f s)" % (time.time() - t0), file=sys.stderr)


# ====================================================================== Part C
def b_count(x):
    """b(x): the number of unit entries u with u = +-2 mod 5"""
    return sum(c for y, c in x.items() if gcd(y, M) == 1 and y % 5 in (2, 3))


def six_classes(x):
    """(n_{A,C}) for A in {{+-1}, {+-2}} mod 5 and C in {{+-1}, {+-3}, {+-2}} mod 7 (theta = 1, zeta_3, zeta_3^2)"""
    out = {}
    for A in ((1, 4), (2, 3)):
        for B in ((1, 6), (3, 4), (2, 5)):
            out[(A, B)] = sum(c for y, c in x.items() if gcd(y, M) == 1 and y % 5 in A and y % 7 in B)
    return out


def part_C(gens):
    print()
    print("Part C. The lemmas of Section 8.3 and the count pattern")
    par = [("n_15 + n_60", [15, 60]), ("n_15 + n_105", [15, 105]), ("n_35 + n_105", [35, 105]),
           ("n_35 + n_140", [35, 140]), ("n_4 + n_12 + n_28 + n_84", [4, 12, 28, 84]), ("n_20", [20]),
           ("n_84", [84]), ("n_420", [420]), ("n_R + n_15", sorted(RSET | {15})), ("len", DIVS)]
    check("Lemma 8.5: %s are even on all generators of B_420" % ", ".join(p[0] for p in par),
          all(n_class(g, p[1]) % 2 == 0 for g in gens for p in par))
    check("Lemma 8.6: n'_420(x) + n_15(x) is even on all generators of B_420 (n'_420 = number of unit entries = +-2 mod 5)",
          all((b_count(g) + n_class(g, [15])) % 2 == 0 for g in gens))
    ok6 = True
    for g in gens:
        sc = six_classes(g)
        for A in ((1, 4), (2, 3)):
            vals = [sc[(A, B)] % 2 for B in ((1, 6), (3, 4), (2, 5))]
            ok6 &= len(set(vals)) == 1
    check("Lemma 8.7: on all generators of B_420, for each A in {{+-1}, {+-2}} mod 5 the three counts n_{A,C}, "
          "C in {{+-1}, {+-2}, {+-3}} mod 7, have the same parity", ok6)
    # (8.11): the norm argument
    w = zadd(ONE, zscal(-1, zconj(ZETA3)))              # 1 - conj(zeta_3)
    normw = zmul(w, zconj(w))
    sixth = [zp(2 * k) for k in range(6)]
    check("(8.11): conj(chi)(3) = conj(zeta_3) for chi = chi_4 eps_5 theta; |1 - conj(zeta_3)|^2 = %s = 3; "
          "the values of chi are sixth roots of unity" % (normw,),
          chi_at(conj(mulc(CHI4, EPS5, THETA)), 3) == zconj(ZETA3) and normw == (3, 0, 0, 0)
          and all(chi_at(mulc(CHI4, EPS5, THETA), u) in sixth for u in range(1, 140) if gcd(u, 140) == 1))
    # cross-check of Lemma 8.7 + 8.6 at the level of unit parts, by enumeration
    prim = [c for c in ODD if conductor(c) == M]
    tab = {uu: [cexp(c, uu) for c in prim] for uu in UNITS}

    def ok_units(ms):
        for k in range(len(prim)):
            bins = [0] * 12
            for uu in ms:
                bins[tab[uu][k]] += 1
            if red(bins) != ZERO:
                return False
        return True

    found = {}
    for size in (1, 2, 3, 4):
        found[size] = [(1,) + rest for rest in combinations_with_replacement(UNITS, size - 1) if ok_units((1,) + rest)]
    evenb = all(sum(1 for uu in ms if uu % 5 in (2, 3)) % 2 == 0 for s in (2, 4) for ms in found[s])
    check("cross-check: the %d odd primitive characters of conductor 420; the multisets of 1, 2, 3, 4 units containing "
          "1 with s_420(chi) = 0 for all of them number %s; n'_420 is even on all of them (multiplication by a unit changes "
          "n'_420 by an even number on multisets of even size)" % (len(prim), [len(found[s]) for s in (1, 2, 3, 4)]),
          len(prim) == 7 and not found[1] and not found[3] and evenb)
    orbits = [{tuple(sorted((t * y) % M for y in w_)) for t in UNITS} for w_ in (W1, W2)]
    check("the Galois orbits [w_1], [w_2] have %d and %d elements and are disjoint (w_1, w_2 are not Galois conjugate)"
          % (len(orbits[0]), len(orbits[1])), len(orbits[0]) == len(orbits[1]) == 96 and not orbits[0] & orbits[1])
    for nm, w_ in (("w_1", W1), ("w_2", W2)):
        x = G.u(w_, M)
        pat = (n_class(x, [420]), n_class(x, [15]), n_class(x, [35]), n_class(x, [60]), n_class(x, [105]),
               n_class(x, [140]), n_class(x, RSET))
        rest = n_class(x, [4, 12, 20, 28, 84])
        check("%s: (n_420, n_15, n_35, n_60, n_105, n_140, n_R) = %s, n_4 + n_12 + n_20 + n_28 + n_84 = %d, n'_420 = %d; "
              "this is the pattern (6, 1, 1, 1, 3, 3, 1) of Theorem 8.1, with 6 + 1 + 1 + 1 + 3 + 3 + 1 = 16 entries"
              % (nm, pat, rest, b_count(x)), pat == (6, 1, 1, 1, 3, 3, 1) and rest == 0 and b_count(x) % 2 == 1)


# ====================================================================== Part D
def part_D():
    print()
    print("Part D. Cross-check: Hodge characters of level 420 with four entries")
    f4 = G.hodge_with_divisor(M, 4)
    check("nu_15 vanishes on the %d Hodge characters of level 420 with four entries having an entry dividing 420 "
          "(nu_15 is invariant under (Z/420)^x, and every orbit contains such a character)" % len(f4),
          all(G.nu(G.u(a, M), 15, M) == 0 for a in f4) and len(f4) > 0)


def main():
    t0 = time.time()
    gens = part_A()
    part_B(gens)
    part_C(gens)
    part_D()
    bad = [nm for nm, ok in checks if not ok]
    print()
    print("all %d checks passed" % len(checks) if not bad else "%d of %d checks FAILED" % (len(bad), len(checks)))
    print("(%s; %.1f s)" % ("python-flint %s" % G.flint.__version__ if G.flint else "python-flint not available: the "
                            "independent index computation was skipped", time.time() - t0), file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
