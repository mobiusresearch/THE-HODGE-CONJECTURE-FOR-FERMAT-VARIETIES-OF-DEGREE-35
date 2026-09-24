# The Hodge conjecture for Fermat varieties of degree 35

Trevin Peterson

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22931340.svg)](https://doi.org/10.5281/zenodo.22931340)

**Paper:** [`paper/main.pdf`](paper/main.pdf) (LaTeX source: [`paper/main.tex`](paper/main.tex))

## Abstract

Let X^n_m ⊂ P^{n+1} be the Fermat variety of dimension n and degree m. We give a proof, independent of Kang's theorem on the generalized Hodge conjecture for Fermat threefolds and fourfolds, of the Hodge conjecture for X^n_m for every n when m = 2^a 3^b 5^c 7^d is not divisible by 420. For cd = 0 this is due to Aoki. Martelotte, Miranda, Movasati and Villaflor, in work posted shortly before the first version of this paper, independently proved the Hodge conjecture for all Fermat varieties of degree less than 65 other than 44, 51 and 52, among them the degree 35, by a different method which rests on Kang's theorem; combined with Aoki's reduction, their result also covers the other degrees above. By Aoki's work on the gap group, what is missing beyond Aoki's results is the algebraicity of the Hodge classes attached to one exceptional character of X^6_35. We obtain it from a coniveau statement: the 24-dimensional rational sub-Hodge structure of level one of H^3(X^3_35, Q) attached to the character (1,2,16,21,30) of μ_35^5 is supported on a divisor. This is a case of the generalized Hodge conjecture which is included in the statement of Kang's theorem; our proof of it is explicit and does not use that theorem. It uses an explicit one-parameter family of curves on X^3_35, a residue formula for the derivative of the Abel–Jacobi map along the family, and a computer-verified interval-arithmetic certificate at one point of the family. For the ten degrees m = 35k, 1 ≤ k ≤ 10, the reduction to this one class is verified by an exact finite computation which does not use Aoki's structure theorem. As consequences we obtain the Hodge conjecture for products of Fermat varieties of one such degree and for abelian varieties of Fermat type of these degrees, among them the Jacobian of the curve y^2 = x^35 − 1. Finally, for m = 420 we show that the Hodge characters in the one remaining class of the gap group have dimension at least 14, a bound which is attained, and deduce the Hodge conjecture for X^n_420 for n ≤ 12 (for n ≤ 4 this is also included in the statement of Kang's theorem).

## Changes in version 1.2

- Added reference to Martelotte–Miranda–Movasati–Villaflor [arXiv:2609.27301] and adjusted the description of the contribution accordingly.
- Added references to S.-J. Kang, *Refined motivic dimension of some Fermat varieties*, Bull. Aust. Math. Soc. 93 (2016), whose Corollary 3.2 states the generalized Hodge conjecture for Fermat threefolds and fourfolds, and to its companion paper Kang, Canad. Math. Bull. 58 (2015). The abstract, the introduction and Remarks 1.2 and 8.9 now say which statements are also included in the statements of these works; the proofs given here do not use them.
- The mathematical content, the proofs and the verification scripts are unchanged.

## Changes in version 1.1

- New Theorem 1.1(E) and Section 8: at the degree m = 420, every Hodge character in the class of the gap group not covered by the representatives of version 1.0 has at least 16 entries (dimension at least 14); two explicit characters of dimension 14 show that this bound is attained; consequently the Hodge conjecture holds for X^n_420 for every n ≤ 12, and for products of Fermat varieties of degree 420 of small total dimension. The proof is by elementary identities for character sums; the case n ≥ 14 at m = 420 and the degrees 420k, k ≥ 2, remain open.
- New verification script `verification/level420.py` with expected output.
- The DOI is printed on the first page. The results and proofs of version 1.0 are unchanged.

## Status

This is a preprint. It has not yet been peer reviewed. Comments and corrections are welcome through GitHub issues.

## Verification

The computational parts of the proof are in [`verification/`](verification/). See [`verification/README.md`](verification/README.md) for what each script checks and how each check maps to the paper.

```
pip install python-flint sympy mpmath
cd verification
python3 certify.py      # interval-arithmetic certificate (Proposition A.1); ends with: ALL CHECKS PASSED
python3 identities.py   # 55 exact checks; ends with: all 55 checks passed
python3 gap.py          # 134 exact checks on the gap group; ends with: all 134 checks passed
python3 level420.py     # 31 exact checks for the level 420; ends with: all 31 checks passed
```

The expected outputs are in `verification/expected_output/`.

## Building the paper

`paper/main.tex` is self-contained and builds with any standard LaTeX distribution (for example `latexmk -pdf main.tex` or `tectonic main.tex`).

## Licence

The paper (`paper/`) is licensed under [CC BY 4.0](LICENSE-paper). The verification code (`verification/`) is licensed under the [MIT licence](LICENSE).

## Citation

Trevin Peterson, *The Hodge conjecture for Fermat varieties of degree 35*, preprint (2026), doi:[10.5281/zenodo.22931340](https://doi.org/10.5281/zenodo.22931340).

This DOI always resolves to the latest version (currently 1.2); version 1.1 is doi:[10.5281/zenodo.22933605](https://doi.org/10.5281/zenodo.22933605) and version 1.0 is doi:[10.5281/zenodo.22931341](https://doi.org/10.5281/zenodo.22931341). See also [`CITATION.cff`](CITATION.cff).
