#!/usr/bin/env python3
"""P6: scalar curvature of g_mu = (1/4) sech^2(mu x) (dx^2 + g_{S^2}) on R x S^2.

Symbolic (sympy), from the Christoffel symbols and the Ricci tensor of the metric
in coordinates (x, theta, phi). Claim: Scal = 4 (2 (1 - mu^2) cosh^2(mu x) + 6 mu^2),
constant exactly when mu = 1 (then g_1 is the round 3-sphere of radius 1/2, Scal = 24).
Prints PASS or FAIL.
"""
import sys
import sympy as sp

x, th, ph = sp.symbols("x theta phi", real=True)
mu = sp.symbols("mu", positive=True)
X = [x, th, ph]
f = sp.sech(mu * x)**2 / 4
g = sp.diag(f, f, f * sp.sin(th)**2)
ginv = g.inv()


def christoffel(a, b, c):
    return sum(ginv[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
               for d in range(3)) / 2


G = [[[sp.simplify(christoffel(a, b, c)) for c in range(3)] for b in range(3)] for a in range(3)]


def ricci(b, c):
    r = 0
    for a in range(3):
        r += sp.diff(G[a][b][c], X[a]) - sp.diff(G[a][b][a], X[c])
        for d in range(3):
            r += G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a]
    return r


scal = sp.simplify(sum(ginv[b, c] * ricci(b, c) for b in range(3) for c in range(3)))
claim = 4 * (2 * (1 - mu**2) * sp.cosh(mu * x)**2 + 6 * mu**2)
diff = sp.simplify((scal - claim).rewrite(sp.exp))
print("Scal =", sp.simplify(scal.rewrite(sp.exp)))
ok1 = diff == 0
ok2 = sp.simplify(sp.diff(scal.subs(mu, 1), x)) == 0 and sp.simplify(scal.subs(mu, 1) - 24) == 0
nonconst = all(sp.simplify(sp.diff(scal.subs(mu, m), x)) != 0
               for m in (sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(2, 3), 2))
print("formula matches:", ok1, "| constant 24 at mu = 1:", ok2, "| nonconstant for mu != 1:", nonconst)
ok = ok1 and ok2 and nonconst
print("P6:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
