#!/usr/bin/env python3
"""Symbolic checks for the lemmas of Theorem A (see the paper, paper.md).

Light: a few seconds to a minute. Needs sympy only. Every check prints
"held" or "FAILED" and the script exits nonzero if any check fails.

Conventions (as in the paper):
  u = (u1, u2, u3, u4) in R^4, z1 = u1 + i u2, z2 = u3 + i u4.
  KS (Hopf) map pi(u) = x with x1 + i x2 = 2 z1 conj(z2), x3 = |z1|^2 - |z2|^2,
  so |x| = |u|^2 =: r. Circle action: (z1, z2) -> e^{i alpha} (z1, z2).
  a^2 = R, Omega(u) = 2 / (1 + |u|^2 / R)  (stereographic factor of S^4_a).
  3D problem: -Lap_3 psi = Z w_R psi, w_R(r) = 2 / (r (1 + r/R)^2).
  DO 1971: [-(1/2) Lap + U_mu] psi = 0, U_mu = -2 v / (r^2 R^2 [(r/R)^mu + (R/r)^mu]^2).
"""

import sys
import sympy as sp

FAILS = []


def report(name, ok, extra=""):
    print(f"{'held  ' if ok else 'FAILED'}  {name}" + (f"  ({extra})" if extra else ""))
    if not ok:
        FAILS.append(name)


u1, u2, u3, u4 = sp.symbols("u1 u2 u3 u4", real=True)
U = [u1, u2, u3, u4]
R = sp.symbols("R", positive=True)
s2 = u1**2 + u2**2 + u3**2 + u4**2          # |u|^2 = r

X1 = 2 * (u1 * u3 + u2 * u4)
X2 = 2 * (u2 * u3 - u1 * u4)
X3 = u1**2 + u2**2 - u3**2 - u4**2
PI = [X1, X2, X3]


def lap4(f):
    return sum(sp.diff(f, v, 2) for v in U)


def grad4(f):
    return [sp.diff(f, v) for v in U]


# ---------------------------------------------------------------- Lemma 1
print("Lemma 1: KS lift")
report("|pi(u)|^2 = |u|^4", sp.expand(X1**2 + X2**2 + X3**2 - s2**2) == 0)
report("components of pi are harmonic", all(sp.expand(lap4(x)) == 0 for x in PI))
G = sp.Matrix(3, 3, lambda i, j: sp.expand(sum(a * b for a, b in zip(grad4(PI[i]), grad4(PI[j])))))
report("horizontal conformality: grad x_i . grad x_j = 4|u|^2 delta_ij",
       sp.simplify(G - 4 * s2 * sp.eye(3)) == sp.zeros(3))
# circle invariance of pi: generator of (z1,z2) -> e^{i a}(z1,z2)
Xgen = lambda f: -u2 * sp.diff(f, u1) + u1 * sp.diff(f, u2) - u4 * sp.diff(f, u3) + u3 * sp.diff(f, u4)
report("pi is invariant under the KS circle", all(sp.expand(Xgen(x)) == 0 for x in PI))
report("circle generator is orthogonal to every grad x_i",
       all(sp.expand(sum(g * c for g, c in zip(grad4(x), [-u2, u1, -u4, u3]))) == 0 for x in PI))
# Laplacian identity on explicit nontrivial test functions
y1, y2, y3 = sp.symbols("y1 y2 y3", real=True)
tests = [sp.exp(y1) * sp.sin(y2) * y3**2, 1 / (1 + y1**2 + y2**2 + y3**2), y1 * y2 * y3 + y3**3,
         sp.exp(-sp.sqrt(y1**2 + y2**2 + y3**2))]
ok = True
for f in tests:
    lap3 = sum(sp.diff(f, v, 2) for v in (y1, y2, y3))
    sub = dict(zip((y1, y2, y3), PI))
    lhs = lap4(f.subs(sub))
    rhs = 4 * s2 * lap3.subs(sub)
    pt = {u1: sp.Rational(3, 7), u2: -sp.Rational(2, 5), u3: sp.Rational(5, 11), u4: sp.Rational(1, 3)}
    ok &= abs(sp.N((lhs - rhs).subs(pt), 40)) < 1e-30
    gl = sum(g**2 for g in grad4(f.subs(sub)))
    gr = 4 * s2 * sum(sp.diff(f, v)**2 for v in (y1, y2, y3)).subs(sub)
    ok &= abs(sp.N((gl - gr).subs(pt), 40)) < 1e-30
report("Lap_4(f o pi) = 4|u|^2 (Lap_3 f) o pi and |grad_4(f o pi)|^2 = 4|u|^2 |grad_3 f|^2 o pi (4 test f)", ok)
J = sp.Matrix([grad4(x) for x in PI])
jac = sp.sqrt(sp.factor((J * J.T).det()))
report("coarea Jacobian sqrt(det(D pi D pi^T)) = 8|u|^3", sp.simplify(jac**2 - 64 * s2**3) == 0)
# radial check of the pushforward measure: int_{R^4} f(|u|^2) du = (pi/4) int_{R^3} f(|x|)/|x| dx
rho, r = sp.symbols("rho r", positive=True)
f_test = sp.exp(-r) * (1 + r**2)
lhs = 2 * sp.pi**2 * sp.integrate(f_test.subs(r, rho**2) * rho**3, (rho, 0, sp.oo))
rhs = sp.pi / 4 * 4 * sp.pi * sp.integrate(f_test / r * r**2, (r, 0, sp.oo))
report("pushforward of du is (pi/4) dx/|x| (radial test)", sp.simplify(lhs - rhs) == 0)
wR = 2 / (r * (1 + r / R)**2)
Om = 2 / (1 + s2 / R)
report("lifted weight 4|u|^2 w_R(|u|^2) = 8/(1+|u|^2/R)^2 = 2 Omega^2",
       sp.simplify(4 * s2 * wR.subs(r, s2) - 2 * Om**2) == 0)

# ---------------------------------------------------------------- Lemma 2
print("Lemma 2: conformal weight and covariance")
report("Yamabe equation -Lap_4 Omega = (2/R) Omega^3", sp.simplify(-lap4(Om) - 2 / R * Om**3) == 0)
Phi = sp.Function("Phi")(*U)
lap_g = Om**-4 * sum(sp.diff(Om**2 * sp.diff(Phi, v), v) for v in U)
cov = -lap4(Om * Phi) - Om**3 * (-lap_g + 2 / R * Phi)
report("covariance -Lap_4(Omega Phi) = Omega^3 (-Lap_g Phi + (2/R) Phi), generic Phi", sp.simplify(cov) == 0)
# scalar curvature of g = Omega^2 delta in dim 4: Scal = -6 Omega^-3 Lap Omega
scal = sp.simplify(-6 * Om**-3 * lap4(Om))
report("Scal(g) = 12/R, so L = -Lap_g + Scal/6 = -Lap_g + 2/a^2", sp.simplify(scal - 12 / R) == 0)
# Laplace-Beltrami on S^4 (unit) acting on sin^j(th) C^{j+3/2}_{k-j}(cos th) h_j:
th = sp.symbols("theta")
ok = True
for k in range(0, 7):
    for j in range(0, k + 1):
        g = sp.sin(th)**j * sp.gegenbauer(k - j, sp.Rational(2 * j + 3, 2), sp.cos(th))
        lb = sp.diff(sp.sin(th)**3 * sp.diff(g, th), th) / sp.sin(th)**3 - j * (j + 2) * g / sp.sin(th)**2
        resid = lb + k * (k + 3) * g
        ok &= all(abs(sp.N(resid.subs(th, tv), 40)) < 1e-30 for tv in [sp.Rational(1, 3), sp.Rational(6, 5), sp.Rational(5, 2)])
report("sin^j C^{j+3/2}_{k-j}(cos) h_j has -Lap_{S^4} eigenvalue k(k+3), 0 <= j <= k <= 6", ok)
kk = sp.symbols("k", integer=True, nonnegative=True)
report("k(k+3) + 2 = (k+1)(k+2)", sp.expand(kk * (kk + 3) + 2 - (kk + 1) * (kk + 2)) == 0)
report("dim H_k(R^5) = sum_{j<=k} (j+1)^2 (k <= 20)",
       all(sp.binomial(k + 4, 4) - sp.binomial(k + 2, 4) == sum((j + 1)**2 for j in range(k + 1)) for k in range(21)))

# ---------------------------------------------------------------- Lemma 3
print("Lemma 3: spectral equivalence (explicit eigenfunctions and asymptotics)")
report("Hardy-type bound r^2 w_R(r) <= R/2: R/2 - r^2 w_R = R (r - R)^2 / (2 (r + R)^2)",
       sp.simplify((R / 2 - r**2 * wR) - R * (r - R)**2 / (2 * (r + R)**2)) == 0)
# DO Eq. (26) at mu = 1/2 equals Omega * (S^4 harmonic) with cos(theta) = (R - r)/(R + r)
ok = True
for l in range(0, 4):
    for nr in range(0, 4):
        psi_do = (r / R)**l * (r / R + 1)**(-(2 * l + 1)) * sp.gegenbauer(nr, sp.Rational(4 * l + 3, 2), (R - r) / (R + r))
        ct = (R - r) / (R + r)
        st = 2 * sp.sqrt(r * R) / (R + r)
        lift = 2 / (1 + r / R) * st**(2 * l) * sp.gegenbauer(nr, sp.Rational(4 * l + 3, 2), ct)
        ratio = sp.simplify(psi_do / lift)
        ok &= ratio.free_symbols <= {R} and sp.simplify(sp.diff(ratio, r)) == 0
report("DO Eq. (26), mu = 1/2, equals Omega times sin^{2l} C^{2l+3/2}_{n_r}(cos theta), l, n_r <= 3", ok)
# radial ODE of the 3D problem is solved by those functions with Z = N(N+1)/(2R)
ok = True
for l in range(0, 4):
    for nr in range(0, 4):
        N = nr + 2 * l + 1
        psi = (r / R)**l * (r / R + 1)**(-(2 * l + 1)) * sp.gegenbauer(nr, sp.Rational(4 * l + 3, 2), (R - r) / (R + r))
        Z = sp.Rational(N * (N + 1), 2) / R
        res = -(sp.diff(psi, r, 2) + 2 / r * sp.diff(psi, r) - l * (l + 1) / r**2 * psi) - Z * wR * psi
        ok &= sp.simplify(res) == 0
report("those functions solve -Lap psi = Z w_R psi with Z = N(N+1)/(2R), N = n_r + 2l + 1 = n + l", ok)
ok = True
for l in range(0, 4):
    for nr in range(0, 3):
        psi = (r / R)**l * (r / R + 1)**(-(2 * l + 1)) * sp.gegenbauer(nr, sp.Rational(4 * l + 3, 2), (R - r) / (R + r))
        a0 = sp.limit(psi / r**l, r, 0)
        ainf = sp.limit(psi * r**(l + 1), r, sp.oo)
        ok &= (a0 != 0 and a0.is_finite) and (ainf != 0 and ainf.is_finite)
report("asymptotics psi ~ c0 r^l at 0 and psi ~ c1 r^{-l-1} at infinity, c0, c1 nonzero", ok)

# ---------------------------------------------------------------- Lemma 5
print("Lemma 5: count")
z1, z2, w1, w2, t = sp.symbols("z1 z2 w1 w2 t")


def invariant_harmonic_dim(k, charge=0):
    """dim of degree-k harmonic polynomials on R^5 = C^2 x R with circle charge `charge`.
    Monomials z^a w^b t^c (w = conj z), charge = |a| - |b|. Lap_5 = 4 sum d_z d_w + d_t^2."""
    mons = []
    for c in range(k + 1):
        rest = k - c
        if (rest - charge) % 2:
            continue
        p = (rest + charge) // 2
        q = (rest - charge) // 2
        if p < 0 or q < 0:
            continue
        for a1 in range(p + 1):
            for b1 in range(q + 1):
                mons.append(z1**a1 * z2**(p - a1) * w1**b1 * w2**(q - b1) * t**c)
    if not mons:
        return 0
    imgs = [sp.expand(4 * (sp.diff(m, z1, w1) + sp.diff(m, z2, w2)) + sp.diff(m, t, 2)) for m in mons]
    targets = sorted({mm for im in imgs for mm in sp.Poly(im, z1, z2, w1, w2, t).monoms()}) if any(imgs) else []
    if not targets:
        return len(mons)
    M = sp.Matrix([[sp.Poly(im, z1, z2, w1, w2, t).coeff_monomial(tm) if im != 0 else 0 for im in imgs]
                   for tm in targets])
    return len(mons) - M.rank()


counts = [invariant_harmonic_dim(k) for k in range(9)]
expect = [(k // 2 + 1)**2 for k in range(9)]
report("independent enumeration of circle-invariant harmonics on R^5, k = 0..8",
       counts == expect, f"got {counts}")
charge1 = [invariant_harmonic_dim(k, 1) for k in range(7)]
print(f"        (charge-1 sector sizes, for the spin remark: {charge1})")
# H_{p,q}(C^2) dims and the Hopf pullback of degree-l solid harmonics
ok = True
for p in range(5):
    for q in range(5):
        mons = [z1**a * z2**(p - a) * w1**b * w2**(q - b) for a in range(p + 1) for b in range(q + 1)]
        imgs = [sp.expand(sp.diff(m, z1, w1) + sp.diff(m, z2, w2)) for m in mons]
        if p == 0 or q == 0:
            d = len(mons)
        else:
            targ = sorted({mm for im in imgs for mm in sp.Poly(im, z1, z2, w1, w2).monoms()})
            M = sp.Matrix([[sp.Poly(im, z1, z2, w1, w2).coeff_monomial(tm) for im in imgs] for tm in targ])
            d = len(mons) - M.rank()
        ok &= d == p + q + 1
report("dim H_{p,q}(C^2) = p + q + 1 for p, q <= 4", ok)
ok = True
for l in range(0, 5):
    # Pullback of the full space: Lap_4 (P o pi) = 0 for every harmonic P of degree l (already by Lemma 1),
    # so check one generic harmonic polynomial per l:
    P = sp.expand((y1 + sp.I * y2)**l)          # harmonic, degree l
    ok &= sp.expand(lap4(P.subs({y1: X1, y2: X2, y3: X3}))) == 0
report("Hopf pullback of degree-l harmonics is harmonic of degree 2l on R^4 (l <= 4)", ok)
report("sum_{l <= floor(k/2)} (2l+1) = (floor(k/2)+1)^2 = floor((N+1)/2)^2 (k <= 30)",
       all(sum(2 * l + 1 for l in range(k // 2 + 1)) == ((k + 2) // 2)**2 for k in range(31)))

# ---------------------------------------------------------------- Lemma 6
print("Lemma 6: Kelvin = pole swap")
iota = [R * v / s2 for v in U]
sub_i = dict(zip(U, iota))
report("pi(R u/|u|^2) = R^2 pi(u)/|u|^4", all(sp.simplify(x.subs(sub_i, simultaneous=True) - R**2 * x / s2**2) == 0 for x in PI))
report("Omega(iota u) = (|u|^2/R) Omega(u)", sp.simplify(Om.subs(sub_i, simultaneous=True) - s2 / R * Om) == 0)
a = sp.sqrt(R)
xi = [2 * R * v / (R + s2) for v in U] + [a * (s2 - R) / (R + s2)]
xi_i = [sp.simplify(c.subs(sub_i, simultaneous=True)) for c in xi]
report("inverse stereographic image: iota is the reflection t -> -t of S^4_a",
       all(sp.simplify(xi_i[i] - xi[i]) == 0 for i in range(4)) and sp.simplify(xi_i[4] + xi[4]) == 0)
report("sum xi^2 = R (the map lands on the sphere of radius sqrt(R))", sp.simplify(sum(c**2 for c in xi) - R) == 0)
report("Tietz weight is Kelvin invariant: (R/r)^4 w_R(R^2/r) = w_R(r)",
       sp.simplify((R / r)**4 * wR.subs(r, R**2 / r) - wR) == 0)
report("4D Kelvin weight (a/|u|)^2 equals 3D Kelvin weight R/r (constant 1)", sp.simplify(a**2 / s2 - R / s2) == 0)

# ---------------------------------------------------------------- Lemma 4
print("Lemma 4: spindle and DO Eq. (25)")
mu, v, tau, sig = sp.symbols("mu v tau sigma", positive=True)
U_mu = -2 * v / (r**2 * R**2 * ((r / R)**mu + (R / r)**mu)**2)
# lifted weight W = 4 r (-2 U_mu) at r = |u|^2 = R e^{2 tau}
W = (4 * r * (-2 * U_mu)).subs(r, R * sp.exp(2 * tau))
Om_mu2 = 4 * mu**2 * R / (R * sp.exp(2 * tau)) / sp.cosh(2 * mu * tau)**2
report("lifted DO weight W_mu = (v / (mu^2 R^3)) Omega_mu^2, Omega_mu^2 = 4 mu^2 R sech^2(2 mu tau)/|u|^2",
       sp.simplify((W - v / (mu**2 * R**3) * Om_mu2).rewrite(sp.exp)) == 0)
report("mu = 1/2: Omega_mu = round Omega", sp.simplify((Om_mu2.subs(mu, sp.Rational(1, 2)) -
                                                          (2 / (1 + sp.exp(2 * tau)))**2).rewrite(sp.exp)) == 0)
# metric: g = Omega_mu^2 |u|^2 (dtau^2 + g_S3) = 4 mu^2 R sech^2(2 mu tau)(dtau^2 + g_S3);
# with cos(theta) = -tanh(2 mu tau): d theta = 2 mu sech(2 mu tau) d tau, sin(theta) = sech(2 mu tau)
dth = sp.diff(sp.acos(-sp.tanh(2 * mu * tau)), tau)
report("g_mu = R (d theta^2 + (2 mu)^2 sin^2 theta g_S3)",
       sp.simplify(dth**2 - 4 * mu**2 / sp.cosh(2 * mu * tau)**2) == 0)
# DO Eq. (26) solves Eq. (24) with v = v_N of Eq. (25), for several mu (numeric residual, 40 digits)
ok = True
worst = 0
for mu_v in [sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(2, 3), sp.Integer(1), sp.Integer(2)]:
    for l in range(0, 3):
        for nr in range(0, 3):
            Nq = nr + l + 1 + (1 / mu_v - 1) * l
            vN = mu_v**2 * (Nq + 1 / (2 * mu_v)) * (Nq + 1 / (2 * mu_v) - 1)  # R = 1
            lam = (2 * l + 1) / (2 * mu_v) + sp.Rational(1, 2)
            psi = r**l * (r**(2 * mu_v) + 1)**(-(2 * l + 1) / (2 * mu_v)) * \
                sp.gegenbauer(nr, lam, (1 - r**(2 * mu_v)) / (1 + r**(2 * mu_v)))
            Uq = -2 * vN / (r**2 * (r**mu_v + r**(-mu_v))**2)
            res = -sp.Rational(1, 2) * (sp.diff(psi, r, 2) + 2 / r * sp.diff(psi, r) - l * (l + 1) / r**2 * psi) + Uq * psi
            for rv in [sp.Rational(1, 3), sp.Rational(7, 5), 4]:
                val = abs(sp.N(res.subs(r, rv), 40)) / (abs(sp.N(Uq.subs(r, rv) * psi.subs(r, rv), 40)) + 1e-40)
                worst = max(worst, float(val))
ok = worst < 1e-25
report("DO Eq. (26) solves Eq. (24) at v = v_N of Eq. (25), mu in {1/3,1/2,2/3,1,2}, l, n_r <= 2",
       ok, f"worst relative residual {worst:.1e}")
# spindle separated equation: Poschl-Teller eigenvalue s(s+1) with s = n_r + (j+1)/(2 mu), j = 2l,
# equals lambda a^2 = v_N / (mu^2 R^2) from Eq. (25)
nr_s, l_s = sp.symbols("n_r l", nonnegative=True)
s_pt = nr_s + (2 * l_s + 1) / (2 * mu)
Nq = nr_s + l_s + 1 + (1 / mu - 1) * l_s
vN = R**2 * mu**2 * (Nq + 1 / (2 * mu)) * (Nq + 1 / (2 * mu) - 1)
report("Poschl-Teller s(s+1), s = n_r + (2l+1)/(2 mu), equals lambda a^2 = v_N/(mu^2 R^2) (DO Eq. 25)",
       sp.simplify(s_pt * (s_pt + 1) - vN / (mu**2 * R**2)) == 0)
report("mu = 1/2: lambda = 2Z with Z = 2 v_N / R^3 = N(N+1)/(2R)",
       sp.simplify((2 * vN / R**3).subs(mu, sp.Rational(1, 2)) -
                   (nr_s + 2 * l_s + 1) * (nr_s + 2 * l_s + 2) / (2 * R)) == 0)
# Poschl-Teller ground and excited states: sech^kappa C^{kappa+1/2}_n(tanh) at depth s(s+1), s = kappa + n
ok = True
for kap in [sp.Rational(1, 2), sp.Rational(3, 2), 2, sp.Rational(5, 3)]:
    for n in range(3):
        s_ = kap + n
        wf = sp.cosh(sig)**(-kap) * sp.gegenbauer(n, kap + sp.Rational(1, 2), sp.tanh(sig))
        res = -sp.diff(wf, sig, 2) + kap**2 * wf - s_ * (s_ + 1) / sp.cosh(sig)**2 * wf
        for sv in [sp.Rational(-3, 2), sp.Rational(1, 3), 2]:
            ok &= abs(sp.N(res.subs(sig, sv), 40)) < 1e-30
report("Poschl-Teller: sech^kappa C^{kappa+1/2}_n(tanh sigma) solves -w'' + kappa^2 w = s(s+1) sech^2 w, s = kappa+n", ok)

print()
print("ALL HELD" if not FAILS else f"FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
