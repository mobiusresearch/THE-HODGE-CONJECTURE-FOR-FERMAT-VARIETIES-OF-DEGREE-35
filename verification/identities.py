"""Exact checks of the finite computations used in the paper
"The Hodge conjecture for Fermat varieties of degree 35".

Requirements: python3 (>= 3.8) and sympy (sympy is used only in Part I).

    python3 identities.py

Part I   the explicit polynomial and determinant identities of Section 3
         ((3.1), (3.2), (3.4)-(3.6), Lemmas 3.3, 3.5, 3.9, 3.10);
Part II  the characters alpha and kappa (Lemmas 3.1, 4.1, Theorem 4.3, Table 2, Remark 1.2(e));
Part III the identity for gamma (Section 7.1, Propositions 7.1, 7.2, (B.1));
Part IV  the identity for delta (Section 7.2, Propositions 7.3, 7.4, Section B.3);
Part V   exceptionality of beta, gamma, delta (Remark 5.9(a), Lemma 5.6, Remark 7.5).
The computations with Aoki's gap group (Section 5) are in gap.py.
All checks are exact (integer arithmetic, or polynomial identities in sympy).  One line is printed
per check; the exit status is 0 if and only if every check passes.
"""
import math
import sys
from collections import Counter
from itertools import permutations

import sympy as sp

checks = []


def check(name, cond):
    checks.append((name, bool(cond)))
    print(("PASS  " if cond else "FAIL  ") + name)


def is_zero(expr):
    return sp.expand(sp.together(expr)) == 0 or sp.simplify(expr) == 0


def det(rows):
    """determinant of a square array of sympy expressions (Leibniz formula)"""
    n = len(rows)
    tot = 0
    for perm in permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        term = (-1) ** inv
        for i in range(n):
            term *= rows[i][perm[i]]
        tot += term
    return tot


# ====================================================================== Part I
print("Part I. Polynomial and determinant identities of Section 3")
s, t = sp.symbols("s t")
p = sp.symbols("p0:5")
q0, q1, b = sp.symbols("q0 q1 b")
c1, c2, c3, c4 = sp.symbols("c1:5")
P = s**5 + p[4] * s**4 * t + p[3] * s**3 * t**2 + p[2] * s**2 * t**3 + p[1] * s * t**4 + p[0] * t**5
Q = s**2 + q1 * s * t + q0 * t**2
R = s**2 + b * s * t + t**2
u = [t**3 * P * Q, c1 * t * P * Q**2, c2 * P**2, c3 * R**5, c4 * s**7 * t * Q]
check("(3.1): u_0, ..., u_4 are binary forms of degree 10",
      all(sp.Poly(sp.expand(x), s, t).is_homogeneous and sp.Poly(sp.expand(x), s, t).total_degree() == 10 for x in u))
E = sp.expand(sum(u).subs(t, 1))
check("(3.2): E_y(s) = (u_0 + ... + u_4)(s, 1) has degree 10 in s, with leading coefficient c2 + c3",
      sp.Poly(E, s).degree() == 10 and sp.expand(E.coeff(s, 10) - (c2 + c3)) == 0)
check("Lemma 3.3: R(0,1) = R(1,0) = 1", R.subs({s: 0, t: 1}) == 1 and R.subs({s: 1, t: 0}) == 1)

# (3.4)-(3.5) with s, t, R, P, Q as independent symbols
S_, T_, R_, P_, Q_ = sp.symbols("S T R_ P_ Q_")
U = [T_**3 * P_ * Q_, c1 * T_ * P_ * Q_**2, c2 * P_**2, c3 * R_**5, c4 * S_**7 * T_ * Q_]
A_ = S_**6 * T_ * R_**3 * P_ * Q_
Psi = c1**2 * c2**16 * c3**21 * c4**30
alpha = (1, 2, 16, 21, 30)
check("(3.5): u_0 u_1^2 u_2^16 u_3^21 u_4^30 = Psi(y) A^35, A = s^6 t R^3 P Q, Psi = c1^2 c2^16 c3^21 c4^30",
      sp.expand(sp.Mul(*[x**e for x, e in zip(U, alpha)]) - Psi * A_**35) == 0)
AA = sp.Poly(sp.expand(s**6 * t * R**3 * P * Q), s, t)
check("(3.4): A = s^6 t R^3 P Q is a binary form of degree 20 = deg(u3 u4), so A/(u3 u4) has degree 0",
      AA.is_homogeneous and AA.total_degree() == 20 == sp.Poly(sp.expand(u[3] * u[4]), s, t).total_degree())

# (3.6): Delta_v as a 3 x 3 determinant, and eta_v = A Delta_v / (u3 u4) ds, in the chart t = 1
Pv, Qv, dP, dQ, dPv, dQv, dc1, dc2 = sp.symbols("Pv Qv dP dQ dPv dQv dc1 dc2")  # P, Q, P', Q', dot P, dot Q, dot c1, dot c2
Delta = (dc1 / c1 + dQv / Qv) * (dP / Pv - dQ / Qv) - (dc2 / c2 + dPv / Pv - dQv / Qv) * dQ / Qv
dot_log = [dPv / Pv + dQv / Qv, dc1 / c1 + dPv / Pv + 2 * dQv / Qv, dc2 / c2 + 2 * dPv / Pv]   # dot u_k / u_k, k = 0, 1, 2
der_log = [dP / Pv + dQ / Qv, dP / Pv + 2 * dQ / Qv, 2 * dP / Pv]                                # u_k' / u_k,   k = 0, 1, 2
check("(3.6): Delta_v = det(1 1 1 / dot u_k/u_k / u_k'/u_k), k = 0, 1, 2",
      is_zero(det([[1, 1, 1], dot_log, der_log]) - Delta))
sv, Rv = sp.symbols("sv Rv")
Aval = sv**6 * Rv**3 * Pv * Qv
u3u4 = c3 * Rv**5 * c4 * sv**7 * Qv
check("(3.6): A Delta_v / (u3 u4) = P Delta_v / (c3 c4 s R^2)  (Definition 3.4)",
      is_zero(Aval * Delta / u3u4 - Pv * Delta / (c3 * c4 * sv * Rv**2)))
check("(3.4) and Lemma 3.10: sum_k alpha_k (exponent of c_k in u_k) gives Psi",
      sp.expand(sp.Mul(*[(1 if k == 0 else [c1, c2, c3, c4][k - 1]) ** alpha[k] for k in range(5)]) - Psi) == 0)

# Lemma 3.5: weights and the contracted Euler identity
check("Lemma 3.5(a): sum(alpha_k - 1) = 65 and 65 + 3 - 2*34 = 0", sum(a - 1 for a in alpha) == 65 and 65 + 3 - 68 == 0)
x = sp.symbols("x0:5")
V = sp.symbols("V0:5")
W = sp.symbols("W0:5")
Fsym = sum(xk**35 for xk in x)
Fi = [sp.diff(Fsym, xk) for xk in x]
e = [[1 if i == j else 0 for j in range(5)] for i in range(5)]


def dV(*vecs):
    return det([list(v) for v in vecs])


def Omega2(i, j, v, w):
    """Omega_ij(v, w) = dV(Euler, d_i, d_j, v, w)"""
    return dV(x, e[i], e[j], v, w)


def Omega1(i, j, k, v):
    """Omega_ijk(v) = (iota_{d_k} Omega_ij)(v) = dV(Euler, d_i, d_j, d_k, v)"""
    return dV(x, e[i], e[j], e[k], v)


ok = True
for i in range(5):
    for j in range(i + 1, 5):
        for k in range(j + 1, 5):
            lhs = Fi[i] * Omega2(j, k, V, W) - Fi[j] * Omega2(i, k, V, W) + Fi[k] * Omega2(i, j, V, W)
            dFV = sum(Fi[l] * V[l] for l in range(5))
            dFW = sum(Fi[l] * W[l] for l in range(5))
            rhs = dFV * Omega1(i, j, k, W) - dFW * Omega1(i, j, k, V) + 35 * Fsym * dV(e[i], e[j], e[k], V, W)
            ok &= sp.expand(lhs - rhs) == 0
check("Lemma 3.5(b): F_i Om_jk - F_j Om_ik + F_k Om_ij = dF ^ Om_ijk + 35 F i_k i_j i_i dV (all 10 triples)", ok)
check("Lemma 3.5(a): Omega_ij(Euler, .) = 0", all(sp.expand(Omega2(i, j, x, W)) == 0 for i in range(5) for j in range(5) if i != j))

# Lemma 3.9: the contraction
a_, d_ = sp.symbols("a0:5"), sp.symbols("d0:5")   # dot u_k/u_k and u_k'/u_k
Vv = [x[k] * a_[k] / 35 for k in range(5)]
Ws = [x[k] * d_[k] / 35 for k in range(5)]
check("Lemma 3.9: dV(Euler, d_3, d_4, dx/dv, dx/ds) = x0 x1 x2 / 35^2 * det(1 1 1 / a_0 a_1 a_2 / d_0 d_1 d_2)",
      sp.expand(dV(x, e[3], e[4], Vv, Ws) - x[0] * x[1] * x[2] / 35**2 * det([[1, 1, 1], list(a_[:3]), list(d_[:3])])) == 0)
check("Lemma 3.9: x^(alpha-1) x0 x1 x2 / (35^2 x3^34 x4^34) / 35^2 = 35^-4 x^alpha / (x3^35 x4^35)",
      sp.simplify(sp.Mul(*[x[k] ** (alpha[k] - 1) for k in range(5)]) * x[0] * x[1] * x[2] / (35**2 * x[3]**34 * x[4]**34) / 35**2
                  - sp.Rational(1, 35**4) * sp.Mul(*[x[k] ** alpha[k] for k in range(5)]) / (x[3]**35 * x[4]**35)) == 0)


# ====================================================================== characters
def units(m):
    return [t for t in range(1, m) if math.gcd(t, m) == 1]


def deg(a, m, t=1):
    """|t a| = (sum of the representatives in {1, ..., m-1} of the entries of t a) / m"""
    tot = sum((t * y) % m for y in a)
    assert tot % m == 0
    return tot // m


def in_A(a, m):
    return all(y % m for y in a) and sum(a) % m == 0


def in_B(a, m):
    """Theorem 2.1(iii): n = len(a) - 2 even and |t a| = n/2 + 1 for all units t"""
    n = len(a) - 2
    return in_A(a, m) and n % 2 == 0 and all(deg(a, m, t) == n // 2 + 1 for t in units(m))


def in_B_odd(a, m):
    """Theorem 2.1(iv): n = len(a) - 2 = 2p + 1 and p + 1 <= |t a| <= p + 2 for all units t"""
    n = len(a) - 2
    pp = (n - 1) // 2
    return in_A(a, m) and n % 2 == 1 and all(pp + 1 <= deg(a, m, t) <= pp + 2 for t in units(m))


def in_D(a, m):
    """a = (d0, -d0, ..., dr, -dr) up to order, all entries nonzero"""
    cnt = Counter(y % m for y in a)
    if cnt[0] or len(a) % 2:
        return False
    for y in list(cnt):
        if (2 * y) % m == 0:
            if cnt[y] % 2:
                return False
        elif cnt[y] != cnt[(-y) % m]:
            return False
    return True


def same(a, c):
    return Counter(a) == Counter(c)


def scale(k, a, m):
    return tuple((k * y) % m for y in a)


def lift(k, a, m):
    """Lemma 2.3(c): the character of level k m with entries k <a_i>"""
    return tuple(k * (y % m) for y in a)


def conjugates(a, m):
    return {scale(t, a, m) for t in units(m)}


def std(pr, i, m):
    """Aoki's standard character sigma_{p,i} (Section 2.2)"""
    d = m // pr
    if pr == 2:
        return (i, i + d, (m - 2 * i) % m, d)
    return tuple([(i + k * d) % m for k in range(pr)] + [(m - pr * i) % m])


# ====================================================================== Part II
print()
print("Part II. alpha, kappa, beta (Lemmas 3.1, 4.1, Theorem 4.3, Table 2, Remark 1.2(e))")
m = 35
T = [t for t in units(m) if deg(alpha, m, t) == 2]
check("alpha in A^3_35 and in B^3_35", in_A(alpha, m) and in_B_odd(alpha, m))
check("|t alpha| = 2 for t in T, 3 otherwise; T = {1,2,3,4,6,9,11,12,18,19,22,27}",
      T == [1, 2, 3, 4, 6, 9, 11, 12, 18, 19, 22, 27] and all(deg(alpha, m, t) == 3 for t in units(m) if t not in T))
check("h^{2,1} = h^{1,2} = 12", len(T) == 12 and len(units(m)) - len(T) == 12)
check("[alpha] has 24 distinct elements (dim M_alpha = 24) and -alpha = 34 alpha", len(conjugates(alpha, m)) == 24 and scale(34, alpha, m) in conjugates(alpha, m))
pair_sums = [(alpha[i] + alpha[j]) % m for i in range(5) for j in range(i + 1, 5)]
check("Lemma 3.1(c) and Section B.1: the ten pair sums are 3,17,22,31,18,23,32,2,11,16 (none is 0)",
      pair_sums == [3, 17, 22, 31, 18, 23, 32, 2, 11, 16] and 0 not in pair_sums)
cols = {"s": (0, 0, 0, 0, 7), "t": (3, 1, 0, 0, 1), "R": (0, 0, 0, 5, 0), "P": (1, 1, 2, 0, 0), "Q": (1, 2, 0, 0, 1)}
check("before (3.4): in u_0 u_1^2 u_2^16 u_3^21 u_4^30 the forms s,t,R,P,Q occur with exponents 210,35,105,35,35",
      [sum(y * z for y, z in zip(alpha, cols[g])) for g in "stRPQ"] == [210, 35, 105, 35, 35])
check("the exponents are those of (3.1): u_k = prod g^(e_g[k]) up to c_k",
      all(sp.expand(U[k] / sp.Mul(*[sym ** cols[g][k] for g, sym in zip("stRPQ", (S_, T_, R_, P_, Q_))])
                    - [1, c1, c2, c3, c4][k]) == 0 for k in range(5)))

kappa = (17, 22, 31)
beta = alpha + kappa
check("kappa in A^1_35", in_A(kappa, m))
check("Lemma 4.1: |t kappa| = 2 for t in T and 1 otherwise; |t alpha| + |t kappa| = 4",
      all(deg(kappa, m, t) == (2 if t in T else 1) for t in units(m)) and all(deg(alpha, m, t) + deg(kappa, m, t) == 4 for t in units(m)))
check("Theorem 4.3: beta = alpha*kappa in B^6_35, dim M_beta = 24", in_B(beta, m) and len(conjugates(beta, m)) == 24)
TABLE2 = {1: (1, 2, 16, 21, 30, 17, 22, 31), 2: (2, 4, 32, 7, 25, 34, 9, 27), 3: (3, 6, 13, 28, 20, 16, 31, 23),
          4: (4, 8, 29, 14, 15, 33, 18, 19), 6: (6, 12, 26, 21, 5, 32, 27, 11), 9: (9, 18, 4, 14, 25, 13, 23, 34),
          11: (11, 22, 1, 21, 15, 12, 32, 26), 12: (12, 24, 17, 7, 10, 29, 19, 22), 18: (18, 1, 8, 28, 15, 26, 11, 33),
          19: (19, 3, 24, 14, 10, 8, 33, 29), 22: (22, 9, 2, 7, 30, 24, 29, 17), 27: (27, 19, 12, 7, 5, 4, 34, 32),
          8: (8, 16, 23, 28, 30, 31, 1, 3), 13: (13, 26, 33, 28, 5, 11, 6, 18), 16: (16, 32, 11, 21, 25, 27, 2, 6),
          17: (17, 34, 27, 7, 20, 9, 24, 2), 23: (23, 11, 18, 28, 25, 6, 16, 13), 24: (24, 13, 34, 14, 20, 23, 3, 9),
          26: (26, 17, 31, 21, 10, 22, 12, 1), 29: (29, 23, 9, 14, 30, 3, 8, 24), 31: (31, 27, 6, 21, 20, 2, 17, 16),
          32: (32, 29, 22, 7, 15, 19, 4, 12), 33: (33, 31, 3, 28, 10, 1, 26, 8), 34: (34, 33, 19, 14, 5, 18, 13, 4)}
check("Table 2: the 24 rows t alpha | t kappa are correct",
      sorted(TABLE2) == units(m) and all(scale(t, beta, m) == TABLE2[t] for t in units(m)))
b1, b2 = (1, 1, 1, 2, 31, 34), (16, 17, 21, 22, 30, 34)
check("Remark 1.2(e): b1, b2 in A^4_35 and b1*b2 = beta*(1,34,1,34) up to order",
      in_A(b1, m) and in_A(b2, m) and same(b1 + b2, beta + (1, 34, 1, 34)))

# ====================================================================== Part III
print()
print("Part III. gamma on X^4_70 (Section 7.1, Propositions 7.1, 7.2, (B.1))")
gamma = (1, 20, 24, 42, 61, 62)
beta2 = lift(2, beta, 35)
varsigma1, varsigma2 = (1, 18, 53, 68), (17, 26, 36, 61)
sigma55 = std(5, 5, 35)
sigma3 = lift(2, sigma55, 35)
eps = (2, 68, 4, 66, 10, 60, 17, 53, 18, 52, 26, 44, 32, 38, 34, 36)
check("2beta = (2,4,32,42,60,34,44,62) (Lemma 2.3(c), e = 2)", beta2 == (2, 4, 32, 42, 60, 34, 44, 62))
check("gamma in B^4_70; dim M_gamma = 24", in_B(gamma, 70) and len(conjugates(gamma, 70)) == 24)
check("|t gamma| = 3 for all t in (Z/70)^x", all(deg(gamma, 70, t) == 3 for t in units(70)))
check("2beta in B^6_70", in_B(beta2, 70))
check("varsigma1, varsigma2 in B^2_70 (|t varsigma_i| = 2 for all t)", in_B(varsigma1, 70) and in_B(varsigma2, 70))
check("sigma_{5,5} = (5,12,19,26,33,10) at level 35 (p = 5, d = 7, i = 5), gcd(5,7) = 1",
      sigma55 == (5, 12, 19, 26, 33, 10) and math.gcd(5, 7) == 1 and 0 < 5 < 7)
check("sigma_{5,5} in B^4_35 (|t sigma_{5,5}| = 3 for all t)", in_B(sigma55, 35))
check("sigma3 = 2 sigma_{5,5} = (10,24,38,52,66,20)", sigma3 == (10, 24, 38, 52, 66, 20))
check("epsilon in D^14_70; its pairs are {2,68},{4,66},{10,60},{17,53},{18,52},{26,44},{32,38},{34,36}",
      in_D(eps, 70) and len(eps) == 16 and all((eps[2 * k] + eps[2 * k + 1]) % 70 == 0 for k in range(8)))
check("gamma*epsilon = 2beta*varsigma1*varsigma2*sigma3 up to order", same(gamma + eps, beta2 + varsigma1 + varsigma2 + sigma3))
check("(B.1): the common multiset is {1,2,4,10,17,18,20,24,26,32,34,36,38,42,44,52,53,60,61,62,66,68}",
      sorted(gamma + eps) == [1, 2, 4, 10, 17, 18, 20, 24, 26, 32, 34, 36, 38, 42, 44, 52, 53, 60, 61, 62, 66, 68])
check("Proposition 7.2: 2beta*varsigma1 in B^10, 2beta*varsigma1*varsigma2 in B^14, 2beta*varsigma1*varsigma2*sigma3 in B^20 (level 70)",
      in_B(beta2 + varsigma1, 70) and in_B(beta2 + varsigma1 + varsigma2, 70) and in_B(beta2 + varsigma1 + varsigma2 + sigma3, 70)
      and len(beta2 + varsigma1 + varsigma2 + sigma3) == 22)
check("Proposition 7.2: r = 4, s = 14 are even (len gamma = 6, len epsilon = 16)", len(gamma) == 6 and len(eps) == 16)
check("Section 7.1: varsigma1*(36,34,35,35) = sigma_{2,18}*sigma_{2,1} and varsigma2*(52,18,35,35) = sigma_{2,26}*sigma_{2,17} "
      "up to order (level 70); (36,34,35,35), (52,18,35,35) are decomposable",
      same(varsigma1 + (36, 34, 35, 35), std(2, 18, 70) + std(2, 1, 70)) and same(varsigma2 + (52, 18, 35, 35), std(2, 26, 70) + std(2, 17, 70))
      and in_D((36, 34, 35, 35), 70) and in_D((52, 18, 35, 35), 70))
check("Section 7.1: sigma3 = sigma_{5,10} up to order (level 70), and gcd(10, 14) = 2", same(sigma3, std(5, 10, 70)) and math.gcd(10, 14) == 2)

# ====================================================================== Part IV
print()
print("Part IV. delta on X^4_210 (Section 7.2, Propositions 7.3, 7.4, Section B.3)")
delta = (2, 9, 129, 142, 168, 180)
gamma3 = lift(3, gamma, 70)
tau = (24, 94, 138, 164)
eps2 = (24, 186, 72, 138)
d47 = scale(47, delta, 210)
check("3gamma = (3,60,72,126,183,186) (Lemma 2.3(c), e = 3) and 3gamma in B^4_210",
      gamma3 == (3, 60, 72, 126, 183, 186) and in_B(gamma3, 210))
check("delta in B^4_210 (|t delta| = 3 for all t)", in_B(delta, 210))
check("tau in B^2_210 (|t tau| = 2 for all t); -tau = (46,72,116,186) up to order",
      in_B(tau, 210) and same(scale(-1, tau, 210), (46, 72, 116, 186)))
check("epsilon' in D^2_210, pairs {24,186}, {72,138}", in_D(eps2, 210) and (24 + 186) % 210 == 0 and (72 + 138) % 210 == 0)
check("47 delta = (94,3,183,164,126,60) mod 210 and gcd(47,210) = 1", d47 == (94, 3, 183, 164, 126, 60) and math.gcd(47, 210) == 1)
check("3gamma*tau = 47delta*epsilon' up to order = {3,24,60,72,94,126,138,164,183,186}",
      same(gamma3 + tau, d47 + eps2) and sorted(gamma3 + tau) == [3, 24, 60, 72, 94, 126, 138, 164, 183, 186])
check("3gamma*tau in B^8_210", in_B(gamma3 + tau, 210))
check("Section B.3: tau = sigma_{3,24} = (24, 94, 164, 138) up to order (level 210)", same(tau, std(3, 24, 210)))
check("Proposition 7.4: 2t = 2 forces t = 1 mod 105, 9t = 9 forces t = 1 mod 70; the 48 characters t delta are distinct",
      all(t % 105 == 1 for t in units(210) if (2 * t) % 210 == 2) and all(t % 70 == 1 for t in units(210) if (9 * t) % 210 == 9)
      and len(units(210)) == 48 and len(conjugates(delta, 210)) == 48)


# ====================================================================== Part V
print()
print("Part V. Exceptionality of beta, gamma, delta (Remark 5.9(a), Lemma 5.6, Remark 7.5)")


def standard_characters(level):
    """all sigma_{p,i} at the given level: primes p | level, level > p, 0 < i < level/p"""
    out = {}
    for pr in range(2, level):
        if level % pr or any(pr % r == 0 for r in range(2, pr)):
            continue
        for i in range(1, level // pr):
            out[(pr, i)] = std(pr, i, level)
    return out


def parity(a, level, support):
    lam = {y % level for z in support for y in (z, -z)}
    return sum(1 for y in a if y % level in lam) % 2


def order(y, level):
    return level // math.gcd(y % level, level)


st35 = standard_characters(35)
check("level 35: the standard characters are sigma_{5,i} (1 <= i <= 6) and sigma_{7,i} (1 <= i <= 4)",
      sorted(st35) == [(5, i) for i in range(1, 7)] + [(7, i) for i in range(1, 5)])
hit = {k: sum(1 for y in v if y % 35 in {1, 34, 6, 29}) for k, v in st35.items()}
check("only sigma_{5,1}, sigma_{5,6}, sigma_{7,1}, sigma_{7,4} have entries in {+-1, +-6}, two each",
      {k: n for k, n in hit.items() if n} == {(5, 1): 2, (5, 6): 2, (7, 1): 2, (7, 4): 2}
      and st35[(5, 1)] == (1, 8, 15, 22, 29, 30) and st35[(5, 6)] == (6, 13, 20, 27, 34, 5)
      and st35[(7, 1)] == (1, 6, 11, 16, 21, 26, 31, 28) and st35[(7, 4)] == (4, 9, 14, 19, 24, 29, 34, 7))
check("varpi vanishes on decomposable characters (the support is symmetric) and varpi(beta) = 1",
      parity(beta, 35, [1, 6]) == 1 and all(parity((d, -d % 35), 35, [1, 6]) == 0 for d in range(1, 35)))
for name, lev, target, entry in (("beta", 35, beta, 21), ("gamma", 70, gamma, 42), ("delta", 210, delta, 168)):
    st = standard_characters(lev)
    nu5 = lambda a: sum(1 for y in a if order(y, lev) == 5) % 2
    check("level %d: gcd(5, %d/5) = 1; nu_5 (parity of the number of entries of order 5) vanishes on all %d standard "
          "characters and on all pairs (y, -y); the only entry of order 5 of %s is %d, so nu_5(%s) = 1"
          % (lev, lev, len(st), name, entry, name),
          math.gcd(5, lev // 5) == 1 and all(nu5(v) == 0 for v in st.values())
          and all(nu5((d, (-d) % lev)) == 0 for d in range(1, lev))
          and [y for y in target if order(y, lev) == 5] == [entry] and nu5(target) == 1)

# ====================================================================== summary
bad = [nm for nm, okk in checks if not okk]
print()
print("all %d checks passed" % len(checks) if not bad else "%d of %d checks FAILED" % (len(bad), len(checks)))
sys.exit(1 if bad else 0)
