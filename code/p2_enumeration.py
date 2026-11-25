#!/usr/bin/env python3
"""P2: dimensions of the circle-invariant harmonic polynomials of degree k on R^5.

Independent of the S^4 picture. R^5 = C^2 x R with coordinates z1, z2, t; the circle
acts by (z1, z2) -> e^{i a}(z1, z2). A monomial z^a zbar^b t^c is invariant iff
|a| = |b|. We build the invariant monomials of degree k, apply the Laplacian
4 (d_z1 d_zbar1 + d_z2 d_zbar2) + d_t^2 as an exact integer matrix, and take its
kernel dimension by exact rank (fractions). Expected: floor(k/2 + 1)^2 = 1, 1, 4, 4, 9, 9, 16, 16, 25.
Prints PASS or FAIL.
"""
import sys
from fractions import Fraction
from itertools import product


def monomials(k):
    out = []
    for c in range(k + 1):
        rest = k - c
        if rest % 2:
            continue
        h = rest // 2  # |a| = |b| = h
        for a1 in range(h + 1):
            for b1 in range(h + 1):
                out.append((a1, h - a1, b1, h - b1, c))  # exponents of z1, z2, zb1, zb2, t
    return out


def laplacian(m):
    z1, z2, w1, w2, t = m
    res = {}
    def add(mon, coef):
        if coef:
            res[mon] = res.get(mon, 0) + coef
    if z1 and w1:
        add((z1 - 1, z2, w1 - 1, w2, t), 4 * z1 * w1)
    if z2 and w2:
        add((z1, z2 - 1, w1, w2 - 1, t), 4 * z2 * w2)
    if t >= 2:
        add((z1, z2, w1, w2, t - 2), t * (t - 1))
    return res


def rank(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rk, ncol = 0, len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def harmonic_invariant_dim(k):
    mons = monomials(k)
    imgs = [laplacian(m) for m in mons]
    targets = sorted({t for im in imgs for t in im})
    if not targets:
        return len(mons)
    mat = [[im.get(t, 0) for im in imgs] for t in targets]
    return len(mons) - rank(mat)


def main():
    got = [harmonic_invariant_dim(k) for k in range(9)]
    want = [1, 1, 4, 4, 9, 9, 16, 16, 25]
    print("circle-invariant harmonics, degrees 0..8:", got)
    print("expected                               :", want)
    ok = got == want
    print("P2:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
