#!/usr/bin/env python3
"""P7: two numbers of the paper that no other script tests.

(a) Section 7: for k = 2 the ratio of the L^2(S^4) norm squared of the harmonic polynomial Y of block l
    to 2 pi^2 = vol(S^3) times the Fischer norm squared of its Cauchy datum on t = 0 is 4/21 (l = 0) and
    4/105 (l = 1). It is computed here directly from monomials (exact integrals over the sphere).
(b) Section 8: the sectors of circle charge 1 (polynomials of bidegree (p, q), p - q = 1, in z, zbar) of the
    degree-k harmonic polynomials on R^5 = C^2 x R have dimensions 0, 2, 2, 6, 6, 12, 12 for k = 0..6.

Run from the repository root:  OMP_NUM_THREADS=1 nice -n 19 python3 code/p7_extra_checks.py
"""
import itertools
import sys

from sympy import Poly, Rational, factorial, gamma, gegenbauer, pi, simplify, symbols, sqrt, S
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

ok = True


def check(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print(("held   " if cond else "FAILED ") + name + (" " + detail if detail else ""))


# ---------------------------------------------------------------- (a) norm ratios, k = 2
y1, y2, y3, y4, t = symbols("y1 y2 y3 y4 t", real=True)
V5 = (y1, y2, y3, y4, t)
V4 = (y1, y2, y3, y4)


def sphere_integral(poly, vars_):
    """Integral of a polynomial over the unit sphere S^(n-1) in R^n, n = len(vars_) (exact)."""
    n = len(vars_)
    total = S(0)
    for mon, c in Poly(poly, *vars_).terms():
        if any(e % 2 for e in mon):
            continue
        b = [Rational(e + 1, 2) for e in mon]
        val = 2
        for bi in b:
            val *= gamma(bi)
        val /= gamma(sum(b))
        total += c * val
    return simplify(total)


def fischer(poly, vars_):
    """Fischer inner product <p, p> = sum over monomials of c^2 * alpha!  (real variables)."""
    tot = S(0)
    for mon, c in Poly(poly, *vars_).terms():
        f = 1
        for e in mon:
            f *= factorial(e)
        tot += c**2 * f
    return tot


def ratio(Y, datum):
    """L^2(S^4) norm^2 of Y over (vol(S^3) times the Fischer norm^2 of its Cauchy datum), vol(S^3) = 2 pi^2.
    The quotient does not depend on the scale of h."""
    return simplify(sphere_integral(Y**2, V5) / (2 * pi**2 * fischer(datum, V4)))


rho2 = y1**2 + y2**2 + y3**2 + y4**2
# l = 0 (j = 0), k = 2: Y = h * rho^2 * C^(3/2)_2(t/rho) with h = 1 is a harmonic polynomial
# |x|^2 C^(3/2)_2(t/|x|) with |x|^2 = rho2 + t^2 and C^(3/2)_2(x) = (15 x^2 - 3)/2
Y0 = (15 * t**2 - 3 * (rho2 + t**2)) / 2
lap = lambda f: sum(f.diff(v, 2) for v in V5)
check("Y (l = 0, k = 2) is harmonic on R^5", simplify(lap(Y0)) == 0)
r0 = ratio(Y0, Y0.subs(t, 0))
# l = 1 (j = 2), k = 2: Y = h = y1^2 + y2^2 - y3^2 - y4^2 (the third KS coordinate)
h1 = y1**2 + y2**2 - y3**2 - y4**2
check("Y (l = 1, k = 2) is harmonic on R^5", simplify(lap(h1)) == 0)
r1 = ratio(h1, h1)
print("ratios computed:", r0, r1)
check("l = 0, k = 2: ratio 4/21", r0 == Rational(4, 21), f"(got {r0})")
check("l = 1, k = 2: ratio 4/105", r1 == Rational(4, 105), f"(got {r1})")


# ---------------------------------------------------------------- (b) charge-1 sector sizes
def charge_dim(k, charge):
    """Dimension of harmonic polynomials of degree k on C^2 x R with p - q = charge (exact rank)."""
    # variables: z1, z2 (weight +1), w1, w2 = conj z (weight -1), t (weight 0)
    mons = [m for m in itertools.product(range(k + 1), repeat=5)
            if sum(m) == k and (m[0] + m[1]) - (m[2] + m[3]) == charge]
    if not mons:
        return 0
    # Laplacian = 4 (d_z1 d_w1 + d_z2 d_w2) + d_t^2, acting monomial -> monomial
    targets = {}
    rows = []

    def add(img, m, coef, col):
        if coef:
            idx = targets.setdefault(img, len(targets))
            rows.append((idx, col, coef))

    for col, m in enumerate(mons):
        a1, a2, b1, b2, c = m
        if a1 and b1:
            add((a1 - 1, a2, b1 - 1, b2, c), m, 4 * a1 * b1, col)
        if a2 and b2:
            add((a1, a2 - 1, b1, b2 - 1, c), m, 4 * a2 * b2, col)
        if c >= 2:
            add((a1, a2, b1, b2, c - 2), m, c * (c - 1), col)
    if not targets:
        return len(mons)
    M = [[0] * len(mons) for _ in range(len(targets))]
    for i, j, v in rows:
        M[i][j] += v
    rank = DomainMatrix([[QQ(x) for x in r] for r in M], (len(targets), len(mons)), QQ).rank()
    return len(mons) - rank


expected = [0, 2, 2, 6, 6, 12, 12]
got = [charge_dim(k, 1) for k in range(7)]
formula = [sum(j + 1 for j in range(1, k + 1, 2)) for k in range(7)]
check("charge-1 sizes for k = 0..6 equal 0, 2, 2, 6, 6, 12, 12", got == expected, f"(got {got})")
check("they equal sum over odd j <= k of (j + 1)", formula == expected)
check("charge 0 recovers 1, 1, 4, 4, 9, 9, 16", [charge_dim(k, 0) for k in range(7)] == [1, 1, 4, 4, 9, 9, 16])

print("P7: PASS" if ok else "P7: FAIL")
sys.exit(0 if ok else 1)
