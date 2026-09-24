# Verification files

These files accompany the paper *The Hodge conjecture for Fermat varieties of degree 35*. They contain the computer-verified part of the proof of Theorem 1.1(A), which is Proposition A.1; exact checks of the identities and memberships used in Sections 3, 4 and 7 and Appendix B; the exact computations with Aoki's gap group in Section 5 (Lemma 5.7, Proposition 5.8, Remarks 5.9, 5.10 and 7.5); and the computations for the level 420 in Section 8 (Theorem 8.1).

## Requirements

- Python 3.8 or later
- `python-flint`, which provides Arb ball arithmetic and FLINT
- `sympy`
- `mpmath`

Install them with `pip install python-flint sympy mpmath`. Tested with Python 3.9.6, python-flint 0.6.0, sympy 1.14.0 and mpmath 1.3.0. Nothing else is needed; the scripts read only the files in `data/`.

`gap.py` and `level420.py` need only the Python standard library. If `python-flint` is installed, it also computes the lattice indices independently with FLINT.

## Files

| File | Contents |
|---|---|
| `certify.py` | The certificate for Proposition A.1 (Appendix A). |
| `identities.py` | Exact checks of the identities and memberships in Sections 3, 4 and 7 and Appendix B. |
| `gap.py` | Exact computations with the gap group `B_m/S_m` (Section 5). |
| `level420.py` | Exact computations for the level 420 (Section 8); it imports `gap.py`. |
| `data/point.json` | The exact value `b*` and the other eleven coordinates of `y*`: 38 decimals as starting values, and the 26 decimals printed in Table 1. |
| `expected_output/*.txt` | The output of the four commands below. |

## `python3 certify.py` (a few seconds)

This checks the family of Section 3.2 at one point:

```
E_y(s) = P Q + c1 P Q^2 + c2 P^2 + c3 R^5 + c4 s^7 Q
P = s^5 + p4 s^4 + ... + p0
Q = s^2 + q1 s + q0
R = s^2 + b s + 1
```

Every inequality is proved with Arb balls at 400 bits. Each ball contains the exact value of the quantity it stands for. The script makes the following checks:

- **(J)** The eleven coefficients `E_0, ..., E_10` have integer coefficients, and their partial derivatives are computed exactly with `sympy`. Both are compared, as exact polynomials, with a separate forward-mode differentiation of `E_y(s)`.
- **(K)** The Krawczyk test (Lemma A.2), with `b = b*` fixed at an exact dyadic value and a box of half-width `2^-100`:
  - `E_0 = ... = E_10 = 0` has exactly one solution `y*` in the box;
  - `J_b(y*)` is invertible, with the row-sum bound `4.18e-27`;
  - `y*` agrees with Table 1 to within `1e-26`.
- **(N)** The seven quantities of Definition 3.2 are nonzero at `y*`, each printed as a ball: `c1, c2, c3, c4, disc R, Res(R,P), Res(R,Q)`.
  - The lower bounds stated in Proposition A.1(b) are also checked.
- **(T)** An enclosure of the tangent vector `v* = v(y*)`, with `v*_b = 1`.
- **(I)** The residue sum `I(y*, v*)` of Definition 3.4:
  - it equals `29.58833207712801039873... - 6.76888773820283381802... i`, so its ball excludes 0;
  - it matches Proposition A.1(c).
- **Cross-check (not part of the certificate).** `I(y*, v*)` is recomputed in floating point with `mpmath`, by numerical integration over small circles centred at the roots of `R`. This does not use `sympy` or Arb. The two values agree to within `1e-39`.

Newton refinement of the starting values is not part of the proof; it only produces the centre of the box. The last block prints the enclosure of `y*` to 30 digits; the printed balls are widened to fit that many digits, and the actual radii are below `1e-56`.

Expected last line: `ALL CHECKS PASSED`, exit status 0. See `expected_output/certify.txt`, which is the standard output; a line with the versions and running time is written to standard error.

## `python3 identities.py` (a few seconds)

This makes 55 exact checks, one printed line each. They use integer arithmetic, or polynomial identities in `sympy`.

- **Part I (Section 3).**
  - The degrees in (3.1) and (3.2), and `c2 + c3` as the coefficient of `s^10`.
  - The identity (3.5) `u_0 u_1^2 u_2^16 u_3^21 u_4^30 = Psi A^35`, and the degree of `A` in (3.4).
  - The determinant form (3.6) of `Delta_v`, and `A Delta_v/(u3 u4) = P Delta_v/(c3 c4 s R^2)`.
  - The weights in Lemma 3.5(a).
  - The contracted Euler identity of Lemma 3.5(b), for all ten triples, as an exact polynomial identity in the coordinates and two tangent vectors.
  - The determinant and normalisation in Lemma 3.9.
- **Part II (`alpha`, `kappa`, `beta`).**
  - Lemma 3.1: `|t alpha|`, the set `T`, `dim M_alpha = 24`, and the pair sums.
  - The exponents `210, 35, 105, 35, 35` of `s, t, R, P, Q` stated before (3.4).
  - Lemma 4.1 and Theorem 4.3.
  - All 24 rows of Table 2.
  - Remark 1.2(e).
- **Part III (`gamma`, Section 7.1).**
  - Every membership in Proposition 7.1.
  - The multiset identity `gamma*epsilon ~ 2beta*varsigma1*varsigma2*sigma3`, and (B.1).
  - The hypotheses of Theorem 2.5: `p = 5`, `d = 7`, `i = 5`, `gcd(5, 7) = 1`.
  - The partial juxtapositions used in Proposition 7.2.
  - `varsigma1*(36,34,35,35) ~ sigma_{2,18}*sigma_{2,1}`, `varsigma2*(52,18,35,35) ~ sigma_{2,26}*sigma_{2,17}` and `sigma3 ~ sigma_{5,10}` at level 70.
- **Part IV (`delta`, Section 7.2).**
  - Every membership in Proposition 7.3.
  - `47 delta` and the identity `3gamma*tau ~ 47delta*epsilon'`.
  - `gcd(47, 210) = 1`, and `tau ~ sigma_{3,24}` at level 210.
  - The 48 distinct conjugates of `delta` (Proposition 7.4).
- **Part V (Remark 5.9(a), Lemma 5.6, Remark 7.5).**
  - The standard characters of level 35, and which of them meet `{±1, ±6}`; the parity `varpi` is odd on `beta`.
  - At the levels 35, 70 and 210, the parity `nu_5` of the number of entries of order 5 vanishes on every standard character and on every pair `(y, -y)`, and is odd on `beta`, `gamma` and `delta` respectively.

Expected last line: `all 55 checks passed`, exit status 0. See `expected_output/identities.txt`.

## `python3 gap.py` (about 20 seconds)

This makes 134 exact checks, one printed line each (124 without `python-flint`). Every claim is proved by exact integer arithmetic in pure Python; FLINT is used only for an independent computation of the same indices. In the notation of Section 5, `R_m` is the free abelian group on the symbols `(y)`, `y` in `Z/m - {0}`, and

```
K_m = { x : sum_y c_y <ty> = (m/2) len(x) for all units t }
B_m = { x in K_m : len(x) even }                       (Lemma 5.1)
S_m = subgroup generated by the pairs (y) + (-y) and the u(sigma_{p,i})
rho_12 = (1,6,8,9), rho_15 = (1,6,10,13), rho_20 = (1,4,17,18),
rho_21 = (1,4,18,19), rho_28 = (1,9,21,25), rho_35 = beta = (1,2,16,21,30,17,22,31)
rho_q^(m) = (m/q) rho_q  (entries multiplied by m/q)
```

- **Part A (Lemmas 5.6 and 5.7).** For `q = 12, 15, 20, 21, 28, 35`:
  - `rho_q` is a Hodge character, and `q` is a product of two elements of `P = {4, 3, 5, 7, 11, ...}`;
  - the order parity `nu_l` (`l = 3, 5, 5, 7, 4, 5`) vanishes on all generators of `S_q`, and `rho_q` has exactly one entry of order `l`;
  - `B_q = S_q + Z u(rho_q)`, and `2 u(rho_q)` lies in `S_q`, so `B_q/S_q = Z/2`.
  - Also: `(1,4,16,9,15,18)` is a Hodge character of level 21 with `nu_7 = 1`, and the entries of `(1,9,25,12,20,24)` add up to 91, which is not divisible by 28.
- **Part B (Proposition 5.8, Remark 5.9).** For `m = 35k`, `1 <= k <= 10`, with `Q(m)` the set of `q` in `{12, 15, 20, 21, 28, 35}` dividing `m`:
  - `Q(m)` is the set of divisors `d > 1` of `m` with `mu'(d) = 1`, and `|Q(m)| = 2^(r-1) - 1`;
  - **(a)** `B_m = S_m + sum_q Z u(rho_q^(m))`;
  - **(b)** `2 u(rho_q^(m))` lies in `S_m`, and `B_m/S_m = (Z/2)^|Q(m)|`, with explicit parity homomorphisms `phi_q` (each set `W_q` is printed, as a union of order classes or as a list of residues);
  - `u(rho_q^(m)) + t u(rho_q^(m))` lies in `S_m` for every unit `t`, so `x + t x` lies in `S_m` for every `x` in `B_m` (Remark 5.9(b));
  - **(c)** `[B_m : S_m + sum_{q != 35} Z u(rho_q^(m))] = 2`;
  - the order parities `nu_l` vanish on `S_m`, and separate the classes of the `rho_q^(m)` except for `m = 280, 315`;
  - FLINT cross-check: `[K_m : S_m] = 2^|Q(m)| [K_m : B_m]`, and similarly for the other subgroups.
- **Part C (Remark 7.5).** `u(gamma) - u(2 beta)` lies in `S_70`, `u(delta) - u(6 beta)` lies in `S_210`, and the entries of order 5 of `beta`, `gamma`, `delta` are `21`, `42`, `168`.
- **Part D (Remark 5.10).** For the ten levels:
  - no Hodge character with four entries lies in the class of `beta`, that is, outside `S_m + sum_{q != 35} Z u(rho_q^(m))`;
  - for odd `m`, the same holds for six entries;
  - for even `m`, the six-entry character `(m/70) gamma` does lie in the class of `beta`.

**How (a) is proved.** Every generator is checked to lie in `B_m`. The `|(Z/m)^x|/2` linear forms defining `K_m` have full rank modulo the prime `2^61 - 1`, so `K_m` has rank `r_m = m - 1 - |(Z/m)^x|/2`. An echelon basis of `L_m + Z (m/2)` (`m` even), resp. `L_m` (`m` odd), is computed with unimodular row operations; it has `r_m` rows. For every prime dividing the product of its pivots, its rank modulo that prime is `r_m`. The only such prime occurs for `m = 315`, where the pivot product is 6. So the subgroup is saturated of rank `r_m`, it equals `K_m`, and hence `L_m = B_m`.

**How Part D is proved.** The homomorphism `phi_35` counts the entries whose order lies in a fixed set, so it is invariant under `(Z/m)^x`. Every orbit of `(Z/m)^x` on Hodge characters contains a character with an entry dividing `m`. All Hodge characters with four, resp. six, entries that have such an entry are listed exactly, by a meet-in-the-middle computation on an exact integer key. `phi_35` vanishes on all of them.

Expected last line: `all 134 checks passed`, exit status 0. See `expected_output/gap.txt`, which is the standard output; a line with the running time is written to standard error.

## `python3 level420.py` (about 1 minute)

This makes 31 exact checks, one printed line each (30 without `python-flint`). Here `m = 420`, `n_N(x)` is the number of entries of order `N`, `nu_15 = n_15 mod 2`, `n'_420(x)` is the number of unit entries `u` with `u = +-2 mod 5`, and

```
L''_420 = S_420 + sum_{q in {12,15,20,21,28,35}} Z u(rho_q^(420))
w_1 = (56,69,176,191,193,199,204,209,210,236,243,253,257,267,296,301)
w_2 = (1,17,28,132,141,168,191,199,236,253,257,267,296,371,387,416)
R   = {2,3,5,6,7,10,14,21,30,42,70,210}   (orders)
```

Characters mod 420 are evaluated exactly in `Z[zeta_12]`.

- **Part A (Lemma 8.2, Proposition 8.3).**
  - `Q(420) = {12,15,20,21,28,35,420}`;
  - the 700 generators of `S_420` lie in `B_420`, and `nu_15` vanishes on all of them and on the six `rho_q^(420)`;
  - `w_1`, `w_2` are Hodge characters with 16 entries and exactly one entry of order 15;
  - `B_420 = L''_420 + Z u(w_1)` (same method as Proposition 5.8(a)), and `2 u(w_1)` lies in `S_420`; so `L''_420` is the kernel of `nu_15` on `B_420` and has index 2;
  - FLINT cross-check: `[K_420 : S_420] = 256`, `[K_420 : L''_420] = 4`.
- **Part B (Lemma 8.4 and (8.2)-(8.12)).**
  - For all 48 odd characters `chi` mod 420 and all `y != 0`, the sum `sum_t conj(chi)(t) (2<ty> - 420)` agrees with the formula in the proof of Lemma 8.4, and `B_{1,chi} != 0`;
  - each of the eleven relations (8.2)-(8.12) has the coefficients printed in the paper (times the stated factor `kappa`), and vanishes on all 707 generators of `B_420` from Part A. In (8.11) the coefficient `1 - conj(zeta_3)` is printed as `2 + zeta_3`.
- **Part C (Lemmas 8.5-8.8, Theorem 8.1(iii)).**
  - The parity relations of Lemma 8.5, the congruence `n'_420 = n_15 mod 2` of Lemma 8.6 and the six-class congruence of Lemma 8.7, each on all generators of `B_420`;
  - the norm facts used in Lemma 8.8(b);
  - cross-check of Lemmas 8.6 and 8.7 by enumeration: the multisets of 2 and 4 units containing 1 that satisfy `s_420(chi) = 0` for the seven odd primitive characters of conductor 420 number 4 and 808, and `n'_420` is even on all of them (none of size 1 or 3);
  - `[w_1]`, `[w_2]` have 96 elements each and are disjoint;
  - `w_1` and `w_2` have the count pattern `(n_420, n_15, n_35, n_60, n_105, n_140, n_R) = (6, 1, 1, 1, 3, 3, 1)`, so the bound 16 of Theorem 8.1(ii) is attained.
- **Part D (cross-check).** `nu_15` vanishes on all 5235 Hodge characters of level 420 with four entries having an entry dividing 420.

Expected last line: `all 31 checks passed`, exit status 0. See `expected_output/level420.txt`; running times are written to standard error.

## Notation

The characters are given explicitly in the scripts, as in the paper:

```
alpha  = (1, 2, 16, 21, 30)                in (Z/35)^5
kappa  = (17, 22, 31)
beta   = alpha*kappa = (1, 2, 16, 21, 30, 17, 22, 31)
gamma  = (1, 20, 24, 42, 61, 62)           in (Z/70)^6
delta  = (2, 9, 129, 142, 168, 180)        in (Z/210)^6
varsigma1 = (1, 18, 53, 68)
varsigma2 = (17, 26, 36, 61)
sigma3 = 2*sigma_{5,5} = (10, 24, 38, 52, 66, 20)
tau    = (24, 94, 138, 164)
epsilon, epsilon' as in Section 7
```
