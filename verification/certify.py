"""Computer-verified certificate for Proposition A.1 of the paper
"The Hodge conjecture for Fermat varieties of degree 35".

Requirements: python3 (>= 3.8), sympy, python-flint (Arb ball arithmetic), mpmath.

    python3 certify.py          the point y* of Appendix A

The family (Section 3.2 of the paper).  Unknowns y = (p0,...,p4, q0, q1, b, c1,...,c4) and
    P = s^5 + p4 s^4 + p3 s^3 + p2 s^2 + p1 s + p0,   Q = s^2 + q1 s + q0,   R = s^2 + b s + 1,
    E_y(s) = P Q + c1 P Q^2 + c2 P^2 + c3 R^5 + c4 s^7 Q = E_0 + E_1 s + ... + E_10 s^10.

What is proved for the point y* (every inequality below is a statement about Arb balls, each of
which contains the exact value of the quantity it encloses):
  (K) Krawczyk test.  With b fixed at the exact dyadic value b* of data/point.json, the eleven
      equations E_0 = ... = E_10 = 0 in the other eleven unknowns have exactly one solution y*
      in an explicit complex box of half-width 2^-100, and the Jacobian J_b(y*) is invertible.
      The Jacobian is the exact symbolic derivative (sympy) of the integer polynomials E_k.
  (N) The seven quantities of Definition 3.2 are nonzero at y*:
        c1, c2, c3, c4, disc R, Res(R,P), Res(R,Q)
      Resultants are Sylvester determinants.
  (T) The tangent vector v* = v(y*) (normalised by v*_b = 1) is enclosed.
  (I) The residue sum I(y*, v*) of Definition 3.4 is enclosed and is nonzero.
The script also compares the enclosures with the numbers printed in the paper (Table 1,
Proposition A.1(b),(c) and Section A.3).  Two further checks:
  (J) the eleven polynomials E_k and their partial derivatives, which the certificate evaluates, agree
      exactly (as polynomials) with an independent forward-mode differentiation of E_y(s);
  a floating-point evaluation of I(y*, v*) by numerical contour integration (mpmath), independent
      of sympy and Arb.  This one is a cross-check only and is NOT part of the certificate.
Exit status 0 if and only if every check passes.
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import mpmath
import sympy as sp
from flint import acb, acb_mat, arb, ctx

ctx.prec = 400
HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ exact polynomials
NAMES = ["p0", "p1", "p2", "p3", "p4", "q0", "q1", "b", "c1", "c2", "c3", "c4"]
SYM = sp.symbols(" ".join(NAMES))
p0, p1, p2, p3, p4, q0, q1, b, c1, c2, c3, c4 = SYM
s = sp.Symbol("s")
P_S = s**5 + p4 * s**4 + p3 * s**3 + p2 * s**2 + p1 * s + p0
Q_S = s**2 + q1 * s + q0
R_S = s**2 + b * s + 1
E_S = sp.expand(P_S * Q_S + c1 * P_S * Q_S**2 + c2 * P_S**2 + c3 * R_S**5 + c4 * s**7 * Q_S)
assert sp.Poly(E_S, s).degree() == 10
COEFFS = [sp.expand(E_S.coeff(s, k)) for k in range(11)]
JAC = [[sp.diff(f, v) for v in SYM] for f in COEFFS]             # exact partial derivatives, 11 x 12
IB = NAMES.index("b")
REST = [k for k in range(12) if k != IB]


def compile_poly(f):
    """list of (integer coefficient, exponent tuple) of a polynomial in the 12 unknowns"""
    if f == 0:
        return []
    out = []
    for e, c in sp.Poly(f, *SYM).terms():
        assert c == int(c)
        out.append((int(c), e))
    return out


C_COEFFS = [compile_poly(f) for f in COEFFS]
C_JAC = [[compile_poly(f) for f in row] for row in JAC]


def ev(terms, y):
    """evaluation of a compiled polynomial at a vector of acb balls (inclusion isotone)"""
    tot = acb(0)
    for c, e in terms:
        mon = acb(c)
        for v, k in zip(y, e):
            if k:
                mon *= v ** k
        tot += mon
    return tot


def F_and_J(y):
    return [ev(t, y) for t in C_COEFFS], [[ev(t, y) for t in row] for row in C_JAC]


# ------------------------------------------------------------------ input
def exact_dyadic(text):
    fr = Fraction(text)
    assert fr.denominator & (fr.denominator - 1) == 0, "b must be a dyadic rational"
    x = arb(fr.numerator) / arb(fr.denominator)
    assert x.is_exact()
    return x


def exact_b(pair):
    z = acb(exact_dyadic(pair[0]), exact_dyadic(pair[1]))
    assert z.real.is_exact() and z.imag.is_exact()
    return z


def as_acb(pair):
    return acb(arb(pair[0]), arb(pair[1]))


def newton(z, bval, iters=40):
    """non-rigorous refinement of the eleven unknowns other than b (midpoints only; not part of the proof)"""
    z = [w.mid() for w in z]
    for _ in range(iters):
        y = z[:IB] + [bval] + z[IB:]
        Fv, J = F_and_J(y)
        Jr = acb_mat([[J[i][k] for k in REST] for i in range(11)])
        dz = Jr.solve(acb_mat([[f] for f in Fv]), algorithm="approx")
        z = [(z[k] - dz[k, 0]).mid() for k in range(11)]
    return z


def box(z, r):
    return acb(arb(z.real.mid(), r), arb(z.imag.mid(), r))


def krawczyk(zmid, bval, r):
    """Lemma A.2 on the box zmid +- r (real and imaginary parts), with b = bval exact."""
    ym = zmid[:IB] + [bval] + zmid[IB:]
    Fm, Jm = F_and_J(ym)
    Kinv = acb_mat([[Jm[i][k] for k in REST] for i in range(11)]).inv().mid()      # approximate inverse
    X = [box(w, r) for w in zmid]
    _, JX = F_and_J(X[:IB] + [bval] + X[IB:])
    JXr = acb_mat([[JX[i][k] for k in REST] for i in range(11)])
    ident = acb_mat([[1 if i == j else 0 for j in range(11)] for i in range(11)])
    A = ident - Kinv * JXr
    dX = acb_mat([[X[k] - zmid[k]] for k in range(11)])
    Kv = acb_mat([[w] for w in zmid]) - Kinv * acb_mat([[f] for f in Fm]) + A * dX
    K = [Kv[k, 0] for k in range(11)]
    inside = all(X[k].contains_interior(K[k]) for k in range(11))
    rowsum = [sum((A[i, j].abs_upper() for j in range(11)), arb(0)) for i in range(11)]
    q = max(float(rs.upper()) for rs in rowsum)
    return inside and all(rs < 1 for rs in rowsum), K, q


# ------------------------------------------------------------------ truncated power series in h = s - rho
def pl(expr, y):
    """coefficients (low -> high) in s of a sympy polynomial in s and the unknowns, evaluated at y"""
    pol = sp.Poly(expr, s)
    return [ev(compile_poly(sp.expand(pol.coeff_monomial(s**k))), y) for k in range(pol.degree() + 1)]


def taylor(coeffs, rho, n):
    out = []
    for j in range(n):
        tot = acb(0)
        for i in range(j, len(coeffs)):
            tot += coeffs[i] * math.comb(i, j) * rho ** (i - j)
        out.append(tot)
    return out


def smul(a, c):
    return [sum((a[i] * c[k - i] for i in range(k + 1)), acb(0)) for k in range(len(a))]


def sinv(a):
    out = [1 / a[0]]
    for k in range(1, len(a)):
        out.append(-sum((a[i] * out[k - i] for i in range(1, k + 1)), acb(0)) / a[0])
    return out


def spow(a, e):
    if e < 0:
        return sinv(spow(a, -e))
    r = [acb(1)] + [acb(0)] * (len(a) - 1)
    for _ in range(e):
        r = smul(r, a)
    return r


def sadd(a, c):
    return [x + z for x, z in zip(a, c)]


def ssub(a, c):
    return [x - z for x, z in zip(a, c)]


def sconst(c, n):
    return [c] + [acb(0)] * (n - 1)


def residue_sum(y, v):
    """I(y, v) of Definition 3.4: the sum over the two roots rho of R of Res_{s=rho} eta_v,
    eta_v = P Delta_v / (c3 c4 s R^2) ds."""
    k = 2                                                  # pole order at the roots of R
    Pc, Qc = pl(P_S, y), pl(Q_S, y)
    dP, dQ = [v[0], v[1], v[2], v[3], v[4]], [v[5], v[6]]
    Pd = [Pc[i] * i for i in range(1, len(Pc))]
    Qd = [Qc[i] * i for i in range(1, len(Qc))]
    bb = y[IB]
    sq = (bb * bb - 4).sqrt()
    roots = [(-bb + sq) / 2, (-bb - sq) / 2]
    cc1, cc2, cc3, cc4 = y[8], y[9], y[10], y[11]
    dc1, dc2 = v[8], v[9]
    total = acb(0)
    for idx in (0, 1):
        rho, rho2 = roots[idx], roots[1 - idx]
        n = k
        Ps, Qs = taylor(Pc, rho, n), taylor(Qc, rho, n)
        iP, iQ = sinv(Ps), sinv(Qs)
        LP, LQ = smul(taylor(Pd, rho, n), iP), smul(taylor(Qd, rho, n), iQ)       # P'/P, Q'/Q
        MP, MQ = smul(taylor(dP, rho, n), iP), smul(taylor(dQ, rho, n), iQ)       # dot P/P, dot Q/Q
        Delta = ssub(smul(sadd(sconst(dc1 / cc1, n), MQ), ssub(LP, LQ)),
                     smul(ssub(sadd(sconst(dc2 / cc2, n), MP), MQ), LQ))
        spart = ([rho, acb(1)] + [acb(0)] * n)[:n]                                 # s = rho + h
        other = ([rho - rho2, acb(1)] + [acb(0)] * n)[:n]                          # s - rho2 = (rho - rho2) + h
        f = smul(spow(spart, -1), spow(other, -k))         # 1 / (s (s - rho2)^2)
        f = smul(f, Ps)
        f = smul(f, Delta)
        total += f[k - 1] / (cc3 * cc4)
    return total


# ------------------------------------------------------------------ the seven nonvanishing quantities
def sylvester_det(f, g):
    """resultant of two polynomials (coefficients low -> high) as the Sylvester determinant"""
    m, n = len(f) - 1, len(g) - 1
    fr, gr = list(reversed(f)), list(reversed(g))
    rows = [[acb(0)] * i + fr + [acb(0)] * (n - 1 - i) for i in range(n)]
    rows += [[acb(0)] * i + gr + [acb(0)] * (m - 1 - i) for i in range(m)]
    return acb_mat(rows).det()


def seven(y):
    Pc, Qc, Rc = pl(P_S, y), pl(Q_S, y), pl(R_S, y)
    return [
        ("c1", y[8]), ("c2", y[9]), ("c3", y[10]), ("c4", y[11]),
        ("disc R", Rc[1] ** 2 - 4 * Rc[0] * Rc[2]),
        ("Res(R,P)", sylvester_det(Rc, Pc)), ("Res(R,Q)", sylvester_det(Rc, Qc)),
    ]


# ------------------------------------------------------------------ floating-point cross-checks (mpmath; not part of the proof)
def mp_c(ball):
    return mpmath.mpc(mpmath.mpf(ball.real.mid().str(60, radius=False)), mpmath.mpf(ball.imag.mid().str(60, radius=False)))


def jacobian_exact_check():
    """(J) exact check of the Jacobian: forward-mode differentiation (dual numbers with sympy
    coefficients) of an independent implementation of E, compared with the sp.diff entries."""
    def mul(f, g):
        out = [(0, 0)] * (len(f) + len(g) - 1)
        for i, (a0, a1) in zip(range(len(f)), f):
            for j, (b0, b1) in zip(range(len(g)), g):
                c0, c1_ = out[i + j]
                out[i + j] = (c0 + a0 * b0, c1_ + a0 * b1 + a1 * b0)
        return out

    def scal(c, f):
        return [(c[0] * a0, c[0] * a1 + c[1] * a0) for a0, a1 in f]

    def add(*fs):
        out = [(0, 0)] * 11
        for f in fs:
            for i, (a0, a1) in zip(range(len(f)), f):
                out[i] = (out[i][0] + a0, out[i][1] + a1)
        return out

    for k in range(12):
        y = [(v, 1 if i == k else 0) for i, v in zip(range(12), SYM)]
        one, zero = (1, 0), (0, 0)
        Pd = [y[0], y[1], y[2], y[3], y[4], one]
        Qd = [y[5], y[6], one]
        Rd = [one, y[7], one]
        R5 = mul(mul(mul(mul(Rd, Rd), Rd), Rd), Rd)
        PQ = mul(Pd, Qd)
        Ed = add(PQ, scal(y[8], mul(PQ, Qd)), scal(y[9], mul(Pd, Pd)), scal(y[10], R5), scal(y[11], [zero] * 7 + Qd))
        for i in range(11):
            if sp.expand(Ed[i][0] - COEFFS[i]) != 0 or sp.expand(Ed[i][1] - JAC[i][k]) != 0:
                return False
    return True


def mp_residue_check(ymid, vmid):
    """I(y, v) by numerical integration of eta_v over small circles centred at the two roots of R"""
    mpmath.mp.dps = 40
    y = [mp_c(w) for w in ymid]
    v = [mp_c(w) for w in vmid]
    pp = y[0:5]; qq = y[5:7]; bb = y[7]; k1, k2, k3, k4 = y[8:12]
    P = lambda z: z**5 + pp[4] * z**4 + pp[3] * z**3 + pp[2] * z**2 + pp[1] * z + pp[0]
    dP = lambda z: 5 * z**4 + 4 * pp[4] * z**3 + 3 * pp[3] * z**2 + 2 * pp[2] * z + pp[1]
    Q = lambda z: z**2 + qq[1] * z + qq[0]
    dQ = lambda z: 2 * z + qq[1]
    R = lambda z: z**2 + bb * z + 1
    Pv = lambda z: v[0] + v[1] * z + v[2] * z**2 + v[3] * z**3 + v[4] * z**4
    Qv = lambda z: v[5] + v[6] * z

    def eta(z):
        D = (v[8] / k1 + Qv(z) / Q(z)) * (dP(z) / P(z) - dQ(z) / Q(z)) - (v[9] / k2 + Pv(z) / P(z) - Qv(z) / Q(z)) * dQ(z) / Q(z)
        return P(z) * D / (k3 * k4 * z * R(z) ** 2)

    sq = mpmath.sqrt(bb * bb - 4)
    roots = [(-bb + sq) / 2, (-bb - sq) / 2]
    others = [mpmath.mpc(0)] + list(mpmath.polyroots([1, qq[1], qq[0]])) + list(mpmath.polyroots([1, pp[4], pp[3], pp[2], pp[1], pp[0]]))
    total = mpmath.mpc(0)
    for i in (0, 1):
        rho = roots[i]
        rad = min(abs(rho - w) for w in others + [roots[1 - i]]) / 3
        f = lambda th: eta(rho + rad * mpmath.expj(th)) * 1j * rad * mpmath.expj(th)
        total += mpmath.quad(f, mpmath.linspace(0, 2 * mpmath.pi, 9)) / (2j * mpmath.pi)
    return total


# ------------------------------------------------------------------ one point
def certify(bval, start):
    z0 = [start[nm] for nm in NAMES if nm != "b"]
    zmid = newton(z0, bval)
    ok_k, K, q = krawczyk(zmid, bval, arb(2) ** (-100))
    ystar = K[:IB] + [bval] + K[IB:]                    # every coordinate ball contains the exact zero y*
    _, JK = F_and_J(ystar)
    Jr = acb_mat([[JK[i][k] for k in REST] for i in range(11)])
    w = Jr.solve(acb_mat([[-JK[i][IB]] for i in range(11)]))      # rigorous enclosure (raises if not provably invertible)
    vstar = [w[k, 0] for k in range(IB)] + [acb(1)] + [w[k, 0] for k in range(IB, 11)]
    nd = seven(ystar)
    ok_n = all(val.abs_lower() > 0 for _, val in nd)
    Ival = residue_sum(ystar, vstar)
    ok_i = bool(Ival.abs_lower() > 0)
    return dict(ok=ok_k and ok_n and ok_i, ok_k=ok_k, ok_n=ok_n, ok_i=ok_i, q=q, zmid=zmid,
                K=K, ystar=ystar, vstar=vstar, nd=nd, I=Ival)


def headline():
    data = json.load(open(os.path.join(HERE, "data", "point.json")))
    bval = exact_b(data["b"])
    start = {k: as_acb(v) for k, v in data["start"].items()}
    table = {k: as_acb(v) for k, v in data["table"].items()}
    r = certify(bval, start)
    ok = True

    def check(label, cond):
        nonlocal ok
        ok &= bool(cond)
        print("  [%s] %s" % ("ok" if cond else "FAIL", label))

    print("Point y* of Appendix A; b* = %s + %s i (exact dyadic)" % tuple(data["b"]))
    print()
    print("(J) the coefficients E_k and the Jacobian (sympy.diff) agree exactly with forward-mode")
    print("    differentiation of an independent implementation of E_y(s)")
    check("exact agreement, all 11 x 12 entries", jacobian_exact_check())
    print()
    print("(K) Krawczyk test on the box of half-width 2^-100 (Lemma A.2)")
    check("K(X) lies in the interior of X and ||I - K J||_inf < 1", r["ok_k"])
    print("      row-sum bound ||I - K J||_inf <= %.2e" % r["q"])
    check("row-sum bound <= 4.2e-27 (Section A.3)", r["q"] <= 4.2e-27)
    radmax = max(float(x.rad()) for x in r["K"])
    print("      largest radius of the enclosure of y*: %.1e" % radmax)
    check("enclosure radii < 1e-29 (Section A.3)", radmax < 1e-29)
    for key, tol in (("start", data["start_tolerance"]), ("table", data["table_tolerance"])):
        ref = start if key == "start" else table
        close = all(abs((r["ystar"][i] - ref[nm]).real).upper() < arb(tol) and
                    abs((r["ystar"][i] - ref[nm]).imag).upper() < arb(tol) for i, nm in zip(range(12), NAMES) if nm != "b")
        check("y* agrees with the %s values of data/point.json within %s%s" % (key, tol, " (Table 1)" if key == "table" else ""), close)
    print()
    print("(N) the seven quantities of Definition 3.2 at y* (balls)")
    for nm, val in r["nd"]:
        print("      %-9s %s" % (nm, val.str(15, radius=True)))
    check("all seven balls exclude 0", r["ok_n"])
    lows = dict((nm, val.abs_lower()) for nm, val in r["nd"])
    stated = [("c1", "1.974"), ("c2", "0.156"), ("c3", "0.156"), ("c4", "2.073"), ("disc R", "3.766"),
              ("Res(R,P)", "2.567"), ("Res(R,Q)", "0.419")]
    check("lower bounds of Proposition A.1(b)", all(lows[nm] >= arb(v) for nm, v in stated))
    diff = r["ystar"][9] + r["ystar"][10]
    check("c2 + c3 = 0 (coefficient of s^10): ball contains 0", diff.contains(0))
    print()
    print("(T) tangent vector v* (v*_b = 1)")
    for nm, val in zip(NAMES, r["vstar"]):
        print("      %-3s %s" % (nm, val.str(20, radius=True)))
    trad = max(float(x.rad()) for x in r["vstar"])
    check("enclosure radii < 1e-48 (Section A.3); largest %.1e" % trad, trad < 1e-48)
    print()
    print("(I) residue sum I(y*, v*) (Definition 3.4)")
    print("      I =", r["I"].str(25, radius=True))
    check("the ball excludes 0", r["ok_i"])
    Iref = acb(arb("29.588332077128010399"), arb("-6.7688877382028338180"))
    err = r["I"] - Iref
    check("agrees with Proposition A.1(c) to within 1e-18", abs(err.real).upper() < arb("1e-18") and abs(err.imag).upper() < arb("1e-18"))
    check("|I| >= 30.35", r["I"].abs_lower() >= arb("30.35"))
    print()
    print("Floating-point cross-check with mpmath (independent of sympy and Arb; NOT part of the proof)")
    Ic = mp_residue_check([x.mid() for x in r["ystar"]], [x.mid() for x in r["vstar"]])
    dev = abs(Ic - mp_c(r["I"].mid()))
    check("I by contour integration = %s (deviation %.1e)" % (mpmath.nstr(Ic, 20), float(dev)), dev < mpmath.mpf(10) ** -25)
    print()
    print("y* (enclosure)")
    for nm, val in zip(NAMES, r["ystar"]):
        print("      %-3s %s" % (nm, val.str(30, radius=True)))
    return ok


def main():
    t0 = time.time()
    ok = headline()
    print()
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    print("(python-flint %s, sympy %s, mpmath %s; %.1f s)" % (__import__("flint").__version__, sp.__version__, mpmath.__version__, time.time() - t0),
          file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
