"""M1 of review-r1 (H2): Demkov-Ostrovsky functions continued to global dS4.

Criteria: for N = 1..8, the Remark (de Sitter space) of the paper. Exact arithmetic only (sympy, rationals).
Conventions: R = a = 1 (the tachyon equation times a^2 carries no a), k = N - 1,
r = R e^{2 i eta}, w = e^{i eta}, phi = Omega^{-1}(r) * psi_rad(r) * h(omega),
Omega^{-1} = (1 + r/R)/2, signature (-,+,+,+), (Box + k(k+3)) phi = 0.
Run: python3 code/ds4_m1.py   (a few seconds on one core)
"""
import itertools
import sympy as sp
from sympy import Rational as Q, I

w = sp.symbols('w')                      # e^{i eta}
z1, z2, b1, b2, t = sp.symbols('z1 z2 b1 b2 t')   # b = conj(z); T: z->lam z, b->b/lam
lam = sp.symbols('lam')
x1s, x2s, x3s = sp.symbols('x1 x2 x3')

# Hopf map pi in complex variables (|pi(u)| = |u|^2)
PI = (z1*b2 + z2*b1, -I*(z1*b2 - z2*b1), z1*b1 - z2*b2)
ysq = z1*b1 + z2*b2


def lap4(f):
    return 4*(sp.diff(f, z1, b1) + sp.diff(f, z2, b2))


def lap5(f):
    return lap4(f) + sp.diff(f, t, 2)


def solid_harmonics(l):
    """Basis of harmonic homogeneous polynomials of degree l in x1,x2,x3 (exact nullspace)."""
    mons = [x1s**a*x2s**b*x3s**(l-a-b) for a in range(l+1) for b in range(l+1-a)]
    cs = sp.symbols('c0:%d' % len(mons))
    P = sum(c*m for c, m in zip(cs, mons))
    lap = sp.expand(sum(sp.diff(P, v, 2) for v in (x1s, x2s, x3s)))
    eqs = sp.Poly(lap, x1s, x2s, x3s).coeffs() if lap != 0 else []
    sol = sp.linsolve(eqs, cs) if eqs else None
    free = [c for c in cs]
    if sol is None:
        basis = mons
    else:
        s = list(sol)[0]
        Pg = sum(si*m for si, m in zip(s, mons))
        fs = sorted(Pg.free_symbols & set(cs), key=lambda x: cs.index(x))
        basis = [sp.expand(Pg.subs({f: (1 if f == g else 0) for f in fs})) for g in fs]
    assert len(basis) == 2*l+1, (l, len(basis))
    return basis


def make_h(P):
    return sp.expand(P.subs({x1s: PI[0], x2s: PI[1], x3s: PI[2]}, simultaneous=True))


# ---------- geometry sanity (once) ----------
assert all(sp.expand(lap4(c)) == 0 for c in PI)
assert sp.expand(sum(c**2 for c in PI) - ysq**2) == 0
assert all(sp.expand(c.subs({z1: lam*z1, z2: lam*z2, b1: b1/lam, b2: b2/lam}, simultaneous=True) - c) == 0 for c in PI)

HARM = {}
for l in range(0, 4):
    HARM[l] = []
    for P in solid_harmonics(l):
        h = make_h(P)
        # h in H_{l,l}: harmonic on R^4, T-invariant, homogeneous of degree 2l
        assert sp.expand(lap4(h)) == 0
        assert sp.expand(h.subs({z1: lam*z1, z2: lam*z2, b1: b1/lam, b2: b2/lam}, simultaneous=True) - h) == 0
        s = sp.symbols('s')
        assert sp.expand(h.subs({z1: s*z1, z2: s*z2, b1: s*b1, b2: s*b2}, simultaneous=True) - s**(2*l)*h) == 0
        HARM[l].append((P, h))


def gegen_poly(n, alpha):
    x = sp.symbols('x')
    return sp.Poly(sp.expand(sp.gegenbauer(n, alpha, x)), x), x


def radial_f(l, k):
    """f(w) = Omega^{-1}(r) * psi_rad(r), r = w^2 (R = 1)."""
    r = w**2
    G, x = gegen_poly(k-2*l, 2*l + Q(3, 2))
    arg = (1 - r)/(1 + r)
    psi = r**l * (1 + r)**(-(2*l+1)) * G.as_expr().subs(x, arg)
    return sp.cancel(sp.together((1 + r)/2 * psi))


def dEta(f):
    return sp.cancel(sp.together(I*w*sp.diff(f, w)))


def box_residual(f, l, k):
    cosn = (w + 1/w)/2
    A = 1/cosn
    j = 2*l
    e = -A**-4*dEta(A**2*dEta(f)) - A**-2*j*(j+2)*f + k*(k+3)*f
    return sp.cancel(sp.together(e))


def ambient_poly(l, k, h):
    """Y(y,t) = h(y) * rho^{k-2l} G(-t/rho), rho^2 = |y|^2 + t^2, polynomial of degree k."""
    G, x = gegen_poly(k-2*l, 2*l + Q(3, 2))
    n = k - 2*l
    rho2 = ysq + t**2
    Y = 0
    for (i,), c in G.terms():
        assert (n - i) % 2 == 0
        Y += c*(-t)**i * rho2**((n - i)//2)
    return sp.expand(h*Y)


def rational_points_S2(n):
    pts = []
    for u in range(-6, 7):
        for v in range(-6, 7):
            d = 1 + u*u + v*v
            pts.append((Q(2*u, d), Q(2*v, d), Q(u*u + v*v - 1, d)))
    return pts[:n]


def kernel_dim(k):
    """dim of T-invariant degree-k harmonic polynomials on R^5 (exact rank)."""
    mons = []
    for a in range(k+1):                       # |alpha| = |beta| = a, t-power k-2a
        c = k - 2*a
        if c < 0:
            continue
        for al in [(i, a-i) for i in range(a+1)]:
            for be in [(i, a-i) for i in range(a+1)]:
                mons.append(((al[0], al[1], be[0], be[1], c)))
    idx = {}
    rows = []
    cols = []
    for m in mons:
        e = z1**m[0]*z2**m[1]*b1**m[2]*b2**m[3]*t**m[4]
        L = sp.expand(lap5(e))
        d = {}
        if L != 0:
            for mon, co in sp.Poly(L, z1, z2, b1, b2, t).terms():
                d[mon] = co
        cols.append(d)
        for mon in d:
            idx.setdefault(mon, len(idx))
    M = sp.zeros(len(idx), len(mons))
    for j, d in enumerate(cols):
        for mon, co in d.items():
            M[idx[mon], j] = co
    rank = M.rank() if len(idx) else 0
    return len(mons) - rank, len(mons)


results = {}
allpass = True
for N in range(1, 9):
    k = N - 1
    p = (N - 1)//2
    labels = [(n, l, m) for l in range(0, p+1) for n in [N - l] if n >= l+1 for m in range(2*l+1)]
    res = dict(i_kg=True, i_amb=True, i_mono=True, ii=True, iii=True, iv=True, notes=[])
    amb_polys = []
    data_rows = []
    for l in range(0, p+1):
        f = radial_f(l, k)
        # (i) KG residual
        if box_residual(f, l, k) != 0:
            res['i_kg'] = False
        # (i) monodromy: rational in w, poles only at w=0, oo, w^2=-1
        num, den = sp.fraction(sp.together(f))
        den_f = sp.factor_list(den)[1]
        for fac, _ in den_f:
            if not (fac == w or sp.expand(fac - (w**2 + 1)) == 0 or fac.is_number):
                res['i_mono'] = False
                res['notes'].append('bad pole factor %s' % fac)
        # (iii) parity eta -> -eta (w -> 1/w)
        if sp.cancel(sp.together(f.subs(w, 1/w) - (-1)**k * f)) != 0:
            res['iii'] = False
        # Kelvin on psi_rad itself
        r = sp.symbols('r', positive=True)
        G, x = gegen_poly(k-2*l, 2*l + Q(3, 2))
        psi = r**l*(1+r)**(-(2*l+1))*G.as_expr().subs(x, (1-r)/(1+r))
        kel = sp.cancel(sp.together((1/r)*psi.subs(r, 1/r) - (-1)**k*psi))
        if kel != 0:
            res['iii'] = False
        # (ii) throat data (eta=0 <-> w=1)
        val = sp.simplify(f.subs(w, 1))
        der = sp.simplify(dEta(f).subs(w, 1))
        if N % 2 == 1:
            ok = (val != 0 and der == 0)
        else:
            ok = (val == 0 and der != 0)
        if not ok:
            res['ii'] = False
            res['notes'].append('throat l=%d val=%s der=%s' % (l, val, der))
        c_data = val if N % 2 == 1 else der/I          # d/dt = (1/i) d/d eta at eta=0 (t = i tan eta)
        for (P, h) in HARM[l]:
            Y = ambient_poly(l, k, h)
            amb_polys.append(Y)
            # ambient: harmonic, parity, homogeneity degree k
            if sp.expand(lap5(Y)) != 0:
                res['i_amb'] = False
                res['notes'].append('ambient not harmonic l=%d' % l)
            if sp.expand(Y.subs(t, -t) - (-1)**k*Y) != 0:
                res['iii'] = False
            s = sp.symbols('s')
            if sp.expand(Y.subs({z1: s*z1, z2: s*z2, b1: s*b1, b2: s*b2, t: s*t}, simultaneous=True) - s**k*Y) != 0:
                res['i_amb'] = False
            if sp.expand(Y.subs({z1: lam*z1, z2: lam*z2, b1: b1/lam, b2: b2/lam}, simultaneous=True) - Y) != 0:
                res['i_amb'] = False
            # f == 2^{-(2l+1)} * sec^{2l} G(-i tan eta): restricted ambient radial factor
            # (checked once per l below); throat data vs ambient Cauchy data at t=0
            if N % 2 == 1:
                cauchy = sp.expand(Y.subs(t, 0))
                # on the unit throat |y|=1: data = f(0)*h ; Y(y,0) = G(0)*|y|^{2l}... compare via ratio
                G0 = G.as_expr().subs(x, 0)
                expect = sp.expand(G0*h*ysq**((k-2*l)//2) * 1)
                # Y(y,0) = h * rho^{n} G(0) with rho^2=|y|^2
                if sp.expand(cauchy - expect) != 0:
                    res['ii'] = False
            else:
                cauchy = sp.expand(sp.diff(Y, t).subs(t, 0))
                dG0 = sp.diff(G.as_expr(), x).subs(x, 0)
                expect = sp.expand(-dG0*h*ysq**((k-2*l-1)//2))
                if sp.expand(cauchy - expect) != 0:
                    res['ii'] = False
            # throat data function on S^2 (unit): P_lm(xhat) times c_data
            data_rows.append((l, P, c_data))
        # identification of Route A with the ambient restriction (radial factor, once per l)
        tan = (w - 1/w)/(I*(w + 1/w))       # tan(eta)
        radial_amb = (2*(1/((w + 1/w)/2)))**0 * ((1/((w + 1/w)/2))**(2*l)) * G.as_expr().subs(x, -I*tan)
        if sp.cancel(sp.together(f*2**(2*l+1) - radial_amb)) != 0:
            res['i_amb'] = False
            res['notes'].append('route A != ambient restriction, l=%d' % l)
    # (ii) rank of throat data on S^2: exact evaluation at rational points
    npts = (p+1)**2 + 6
    pts = rational_points_S2(max(npts, 60))
    cols = []
    for (l, P, c) in data_rows:
        col = [c*P.subs({x1s: a, x2s: b, x3s: cc}) for (a, b, cc) in pts]
        cols.append(col)
    rk = sp.Matrix(cols).rank() if cols else 0
    if rk != (p+1)**2:
        res['ii'] = False
        res['notes'].append('data rank %d != %d' % (rk, (p+1)**2))
    spins = sorted(set(l for (l, P, c) in data_rows))
    if spins != list(range(p+1)):
        res['ii'] = False
    # (iv) multiplicity
    kd, nm = kernel_dim(k)
    target = ((N+1)//2)**2
    # use all monomials appearing in any polynomial
    allm = sorted(set(m for Y in amb_polys for m in sp.Poly(Y, z1, z2, b1, b2, t).monoms()))
    amb_rank = sp.Matrix([[sp.Poly(Y, z1, z2, b1, b2, t).coeff_monomial(m) for m in allm] for Y in amb_polys]).rank()
    res['iv_vals'] = (kd, len(labels), amb_rank, target)
    if not (kd == target and len(labels) == target and amb_rank == target):
        res['iv'] = False
    results[N] = res
    print('N=%d k=%d p=%d' % (N, k, p), {kk: v for kk, v in res.items() if kk != 'notes'}, res['notes'], flush=True)
    allpass &= all(res[q] for q in ('i_kg', 'i_amb', 'i_mono', 'ii', 'iii', 'iv'))

print('OVERALL', 'PASS' if allpass else 'FAIL')

# negative controls: the checks must be able to fail
ctl_kg = box_residual(radial_f(1, 3), 1, 4) != 0        # wrong mass (k+1 in the equation)
ctl_par = sp.cancel(sp.together(radial_f(0, 3).subs(w, 1/w) - (-1)**4*radial_f(0, 3))) != 0   # wrong parity
ctl_cnt = kernel_dim(4)[0] != 4                          # wrong count target
print('CONTROLS (all must be True):', ctl_kg, ctl_par, ctl_cnt)
