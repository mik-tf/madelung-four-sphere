"""Checks for Theorem D: zero-energy Green's function of
-Lap - Z w_R on R^3, w_R = 2/(r (1 + r/R)^2), and the addition formula for each
Madelung group).

Every check prints 'held' or 'FAILED' with the measured error.
Part 1 is symbolic (sympy), part 2 numerical (mpmath, numpy, scipy).
The independent reference for the Green's function is the partial-wave sum built
from the radial Legendre solutions (no S^4, no KS map), as in payoff_green_check.py,
now also for Z R < -1/8 (complex s) and Z = 0.
Run: OMP_NUM_THREADS=1 nice -n 19 python3 code/theorem_d_check.py
"""
import numpy as np, mpmath as mp, sympy as sp
from scipy.special import lpmv, eval_legendre, eval_gegenbauer
from scipy.integrate import quad

mp.mp.dps = 25
rng = np.random.default_rng(11)
R = 1.7
ALL = []


def report(name, err, tol):
    ok = bool(err <= tol)
    ALL.append(ok)
    print(f"[{'held' if ok else 'FAILED'}] {name}: {float(err):.2e} (tol {tol:.0e})")


# ------------------------------------------------------------------ part 1: symbolic
def symbolic():
    t, w, s = sp.symbols('t w s')
    # (S1) g(t) = K F(1-s, s+2; 2; (1+t)/2) solves the zonal S^4 equation
    #      (1-t^2) g'' - 4 t g' + (s(s+1) - 2) g = 0: it is the hypergeometric equation.
    F = sp.Function('F')
    hyp = w * (1 - w) * F(w).diff(w, 2) + (2 - 4 * w) * F(w).diff(w) - (1 - s) * (s + 2) * F(w)
    G = F((1 + t) / 2)
    zon = (1 - t ** 2) * G.diff(t, 2) - 4 * t * G.diff(t) + (s * (s + 1) - 2) * G
    zon_w = zon.subs(t, 2 * w - 1)
    report("S1 zonal S^4 equation = hypergeometric equation (c=2, a=1-s, b=s+2)",
           1 if sp.simplify(sp.expand(zon_w.doit() - hyp)) != 0 else 0, 0)
    # (S1b) the same equation is the derivative of Legendre's equation: h'' Legendre => h' zonal
    h = sp.Function('h')
    leg = (1 - t ** 2) * h(t).diff(t, 2) - 2 * t * h(t).diff(t) + s * (s + 1) * h(t)
    gg = h(t).diff(t)
    zon2 = (1 - t ** 2) * gg.diff(t, 2) - 4 * t * gg.diff(t) + (s * (s + 1) - 2) * gg
    report("S1b d/dt(Legendre eq) = zonal S^4 eq for h'", 1 if sp.simplify(zon2 - leg.diff(t)) != 0 else 0, 0)
    # (S2) Z = 0: g0 = 1/(8 pi^2 (1-t)) solves it with s(s+1) = 0
    g0 = 1 / (8 * sp.pi ** 2 * (1 - t))
    report("S2 g0 solves the s(s+1)=0 equation",
           1 if sp.simplify((1 - t ** 2) * g0.diff(t, 2) - 4 * t * g0.diff(t) - 2 * g0) != 0 else 0, 0)
    # (S3) Hermitian product: 2 |<u,v>|^2 = r r' + x.y (Hopf form of the KS map)
    a1, a2, a3, a4, b1, b2, b3, b4 = sp.symbols('a1:5 b1:5', real=True)
    z1, z2, w1, w2 = a1 + sp.I * a2, a3 + sp.I * a4, b1 + sp.I * b2, b3 + sp.I * b4
    def ks(p, q):
        c = 2 * p * sp.conjugate(q)
        return sp.Matrix([sp.re(c), sp.im(c), sp.expand(p * sp.conjugate(p) - q * sp.conjugate(q))])
    x, y = ks(z1, z2), ks(w1, w2)
    r = sp.expand(z1 * sp.conjugate(z1) + z2 * sp.conjugate(z2))
    rp = sp.expand(w1 * sp.conjugate(w1) + w2 * sp.conjugate(w2))
    hp = z1 * sp.conjugate(w1) + z2 * sp.conjugate(w2)
    lhs = sp.expand(2 * hp * sp.conjugate(hp))
    rhs = sp.expand(r * rp + (x.T * y)[0])
    report("S3 2|<u,v>|^2 = r r' + x.y", 1 if sp.simplify(lhs - rhs) != 0 else 0, 0)
    report("S3b |x| = |u|^2", 1 if sp.simplify(sp.expand((x.T * x)[0]) - r ** 2) != 0 else 0, 0)
    # (S4) stereographic map of S^4_a, a^2 = R: unit vector (sin th u/|u|, cos th) with
    #      cos th = (r-R)/(r+R), sin th = 2 sqrt(rR)/(r+R); chordal identity
    #      |xi(u) - xi(v)|^2 = Om(u) Om(v) |u - v|^2.
    p, q, d, Rs = sp.symbols('p q d R', positive=True)  # p=|u|^2, q=|v|^2, d=u.v
    a = sp.sqrt(Rs)
    def xi_dot(p1, q1, d1):  # xi(u).xi(v) from |u|^2, |v|^2, u.v
        return (4 * Rs ** 2 * d1) / ((Rs + p1) * (Rs + q1)) + a ** 2 * (p1 - Rs) * (q1 - Rs) / ((Rs + p1) * (Rs + q1))
    chord = xi_dot(p, p, p) + xi_dot(q, q, q) - 2 * xi_dot(p, q, d)
    Om = lambda m: 2 / (1 + m / Rs)
    report("S4 |xi(u)|^2 = a^2", 1 if sp.simplify(xi_dot(p, p, p) - Rs) != 0 else 0, 0)
    report("S4b chordal |xi(u)-xi(v)|^2 = Om Om |u-v|^2",
           1 if sp.simplify(chord - Om(p) * Om(q) * (p + q - 2 * d)) != 0 else 0, 0)
    cth = (p - Rs) / (p + Rs); sth = 2 * sp.sqrt(p * Rs) / (p + Rs)
    report("S4c sin^2 + cos^2 = 1 and 2 a |u|/(R + |u|^2) = sin th",
           1 if sp.simplify(cth ** 2 + sth ** 2 - 1) != 0 or sp.simplify(2 * a * sp.sqrt(p) / (Rs + p) - sth) != 0 else 0, 0)
    # (S5) Hostler variables: rho_pm = r + r' -/+ ... ; rho_+ rho_- = 2 (r r' + x.y) and
    #      r + r' - sqrt(rho_+ rho_-) = (sqrt rho_+ - sqrt rho_-)^2 / 2
    r_, rp_, D = sp.symbols('r rp D', positive=True)
    xy = (r_ ** 2 + rp_ ** 2 - D ** 2) / 2
    rhp, rhm = r_ + rp_ + D, r_ + rp_ - D
    report("S5 rho_+ rho_- = 2 (r r' + x.y)", 1 if sp.simplify(rhp * rhm - 2 * (r_ * rp_ + xy)) != 0 else 0, 0)
    report("S5b (sqrt rho+ -/+ sqrt rho-)^2/2 = r + r' -/+ sqrt(rho+ rho-)",
           1 if sp.simplify(sp.expand((sp.sqrt(rhp) - sp.sqrt(rhm)) ** 2 / 2) - (r_ + rp_ - sp.sqrt(rhp) * sp.sqrt(rhm))) != 0 else 0, 0)
    # (S6) Hopf average identity, exact, l = 0..8: (1/2pi) int U_{2l}(c cos a) da = P_l(2c^2 - 1),
    #      odd degrees average to 0.
    c, al = sp.symbols('c alpha', real=True)
    bad = 0
    for n in range(0, 17):
        avg = sp.integrate(sp.expand(sp.chebyshevu(n, c * sp.cos(al))), (al, 0, 2 * sp.pi)) / (2 * sp.pi)
        target = sp.legendre(n // 2, 2 * c ** 2 - 1) if n % 2 == 0 else 0
        bad += sp.simplify(avg - target) != 0
    report("S6 Hopf average identity, degrees 0..16 (l = 0..8), exact", bad, 0)
    # (S6b) the generating-function step: average of 1/(1 - 2 h c cos a + h^2)
    h_ = sp.symbols('h', positive=True)
    ser_l = sp.series(1 / sp.sqrt(1 - 2 * h_ ** 2 * (2 * c ** 2 - 1) + h_ ** 4), h_, 0, 14).removeO()
    ser_r = sum(sp.legendre(l, 2 * c ** 2 - 1) * h_ ** (2 * l) for l in range(7))
    report("S6b (1 - 2H x + H^2)^(-1/2), H = h^2, x = 2c^2-1: Legendre series to h^13",
           1 if sp.simplify(sp.expand(ser_l - ser_r)) != 0 else 0, 0)


# ------------------------------------------------------------------ part 2: numerical
def s_of(ZR):
    return (-1 + mp.sqrt(1 + 8 * mp.mpf(ZR))) / 2  # complex (Re s = -1/2) when ZR < -1/8


def g(t, s):
    """unit S^4 kernel, hypergeometric form; s = 0 handled exactly."""
    if s == 0:
        return 1 / (8 * mp.pi ** 2 * (1 - t))
    return s * (s + 1) * mp.hyp2f1(1 - s, s + 2, 2, (1 + t) / 2) / (16 * mp.pi * mp.sin(mp.pi * s))


def g_legendre(t, s):  # the form of the statement: P_s'(-t) / (8 pi sin pi s)
    return mp.diff(lambda z: mp.legenp(s, 0, z, type=2), -t) / (8 * mp.pi * mp.sin(mp.pi * s))


def gprime(t, s):
    return s * (s + 1) * (1 - s) * (s + 2) * mp.hyp2f1(2 - s, s + 3, 3, (1 + t) / 2) / (64 * mp.pi * mp.sin(mp.pi * s))


ZRS = (0.3, 1.7, 4.2, -0.05, -0.125, -0.3, -2.0)


def check_unit_kernel():
    e = 0
    for ZR in ZRS:
        s = s_of(ZR)
        for t in (-0.9, -0.2, 0.4, 0.95):
            e = max(e, abs(g(t, s) - g_legendre(t, s)) / abs(g(t, s)))
    report("N1 hypergeometric form = P_s'(-t)/(8 pi sin pi s) (numerical derivative)", e, 1e-12)
    e = 0
    for ZR in ZRS:
        s = s_of(ZR)
        e = max(e, abs(mp.im(g(0.3, s))))
    report("N1b g real for all ZR (incl. complex s, ZR < -1/8): max |Im g|", e, 1e-20)
    e = 0
    for ZR in ZRS:
        s = s_of(ZR)
        t = 1 - mp.mpf(10) ** -12
        e = max(e, abs(g(t, s) * 8 * mp.pi ** 2 * (1 - t) - 1), abs(gprime(t, s) * 8 * mp.pi ** 2 * (1 - t) ** 2 - 1))
        d = abs(mp.diff(lambda z: g(z, s), 0.37) - gprime(0.37, s)) / abs(gprime(0.37, s))
        e = max(e, d)
    report("N2 g ~ 1/(8 pi^2 (1-t)), g' ~ 1/(8 pi^2 (1-t)^2) at 1-t = 1e-12; g' formula", e, 1e-9)
    # Funk-Hecke eigenvalues: 2 pi^2 int g P_k (1-t^2) dt = 1/((k+1)(k+2) - s(s+1))
    e = 0
    for ZR in ZRS + (0.0,):
        s = s_of(ZR) if ZR != 0 else mp.mpf(0)
        for k in range(7):
            Pk = lambda t: mp.gegenbauer(k, 1.5, t) / mp.gegenbauer(k, 1.5, 1)
            # within 1e-18 of t = 1 use the limit (1+t)(1-t) g -> 2/(8 pi^2) (Step 2); hyp2f1 overflows at w = 1
            fh = lambda t: (1 + t) / (8 * mp.pi ** 2) if 1 - t < mp.mpf(10) ** -18 else g(t, s) * Pk(t) * (1 - t ** 2)
            I = 2 * mp.pi ** 2 * mp.quad(fh, [-1, 0, 0.9, 0.999, 1])
            ref = 1 / ((k + 1) * (k + 2) - 2 * mp.mpf(ZR))
            e = max(e, abs(I - ref) / abs(ref))
    report("N3 Funk-Hecke eigenvalues on H_k, k = 0..6, 8 couplings", e, 1e-12)


# --- geometry
def sph(x):
    r = np.linalg.norm(x); return r


def lift(x):
    r = np.linalg.norm(x)
    z1 = np.sqrt((r + x[2]) / 2)
    z2 = np.conj((x[0] + 1j * x[1]) / 2) / z1
    return z1, z2


def S4(z1, z2):
    n = abs(z1) ** 2 + abs(z2) ** 2
    sc = np.sqrt(R) * 2 / (R + n)
    return np.array([sc * z1.real, sc * z1.imag, sc * z2.real, sc * z2.imag, (n - R) / (n + R)])


def Om(r): return 2 / (1 + r / R)


def cosg_direct(x, y, a):
    z1, z2 = lift(x); w1, w2 = lift(y); e = np.exp(1j * float(a))
    return float(S4(z1, z2) @ S4(e * w1, e * w2))


def AB(x, y):
    r, rp = sph(x), sph(y)
    ct, ctp = (r - R) / (r + R), (rp - R) / (rp + R)
    st, stp = np.sqrt(1 - ct ** 2), np.sqrt(1 - ctp ** 2)
    cosTh = x @ y / (r * rp)
    return ct * ctp, st * stp * np.sqrt((1 + cosTh) / 2)


def rand_pair(sep=False):
    while True:
        x = rng.normal(size=3); y = rng.normal(size=3)
        x *= np.exp(rng.uniform(-1, 1)) * R / np.linalg.norm(x)
        y *= np.exp(rng.uniform(-1, 1)) * R / np.linalg.norm(y)
        if not sep or abs(np.log(np.linalg.norm(x) / np.linalg.norm(y))) > 0.4:
            return x, y


def check_geometry():
    e = 0
    for _ in range(40):
        x, y = rand_pair()
        z1, z2 = lift(x); w1, w2 = lift(y)
        r, rp = sph(x), sph(y)
        hp = z1 * np.conj(w1) + z2 * np.conj(w2)
        A, B = AB(x, y)
        phi = np.angle(hp)
        for a in np.linspace(0, 2 * np.pi, 13):
            e = max(e, abs(cosg_direct(x, y, a) - (A + B * np.cos(a - phi))))
        D = np.linalg.norm(x - y); rho_p, rho_m = r + rp + D, r + rp - D
        cp, cm = A + B, A - B
        e = max(e, abs(1 - cp - Om(r) * Om(rp) * (np.sqrt(rho_p) - np.sqrt(rho_m)) ** 2 / (4 * R)),
                abs(1 - cm - Om(r) * Om(rp) * (np.sqrt(rho_p) + np.sqrt(rho_m)) ** 2 / (4 * R)))
    report("N4 cos gamma_a = A + B cos(a - phi) exactly; 1 - c_pm in Hostler variables", e, 1e-12)
    # payoff identity: fibre average of 1/|u - e^{ia} v|^2 = 1/|x - y|
    e = 0
    for _ in range(10):
        x, y = rand_pair()
        z1, z2 = lift(x); w1, w2 = lift(y)
        u = np.array([z1, z2]); v = np.array([w1, w2])
        f = lambda a: 1 / np.sum(abs(u - np.exp(1j * a) * v) ** 2)
        avg = quad(f, 0, 2 * np.pi, epsabs=0, epsrel=1e-13, limit=200)[0] / (2 * np.pi)
        e = max(e, abs(avg * np.linalg.norm(x - y) - 1))
    report("N5 KS: fibre average of |u - e^{ia} v|^-2 = |x - y|^-1", e, 1e-11)


# --- Green's function
def G_radial(x, y, ZR, lmax=120):
    s = s_of(ZR)
    r, rp = sph(x), sph(y); cosTh = x @ y / (r * rp)
    sg = lambda rr: mp.log(rr / R) / 2
    tot = mp.mpf(0)
    for l in range(lmax):
        k = 2 * l + 1
        y1 = lambda t: mp.legenp(s, -k, -mp.tanh(t), type=2)
        y2 = lambda t: mp.legenp(s, -k, mp.tanh(t), type=2)
        t0 = mp.mpf(0.3)
        W = (mp.diff(y1, t0) * y2(t0) - y1(t0) * mp.diff(y2, t0)) / 2
        lo, hi = (r, rp) if r < rp else (rp, r)
        gl = y1(sg(lo)) * y2(sg(hi)) / mp.sqrt(lo * hi) / W
        term = (2 * l + 1) / (4 * mp.pi) * eval_legendre(l, cosTh) * gl
        tot += term
        if l > 8 and abs(term) < mp.mpf(10) ** -18 * abs(tot): break
    return tot


def G_closed(x, y, ZR, direct=False):
    s = s_of(ZR) if ZR != 0 else mp.mpf(0)
    r, rp = sph(x), sph(y)
    A, B = AB(x, y)
    if direct:
        avg = mp.quad(lambda a: g(cosg_direct(x, y, a), s), np.linspace(0, 2 * np.pi, 9)) / (2 * mp.pi)
    else:
        avg = mp.quad(lambda a: g(A + B * mp.cos(a), s), [0, mp.pi / 2, mp.pi]) / mp.pi
    return mp.pi * Om(r) * Om(rp) / R * avg


def G_extremes(x, y, ZR):  # (1/R) Om Om int_{c-}^{c+} g(t) dt / sqrt((c+ - t)(t - c-))
    s = s_of(ZR)
    r, rp = sph(x), sph(y); A, B = AB(x, y); cp, cm = A + B, A - B
    # t = (cp+cm)/2 - (cp-cm)/2 cos(phi) removes the endpoint singularities; written out directly:
    I = mp.quad(lambda t: g(t, s) / mp.sqrt((cp - t) * (t - cm)), [cm, A, cp])
    return Om(r) * Om(rp) / R * I


def check_green():
    worst = 0
    for ZR in ZRS:
        e1 = e2 = 0
        for i in range(4):
            x, y = rand_pair(sep=True)
            ref = G_radial(x, y, ZR)
            cl = G_closed(x, y, ZR)
            e1 = max(e1, abs(cl - ref) / abs(ref))
            if i < 2:
                e2 = max(e2, abs(G_closed(x, y, ZR, direct=True) - ref) / abs(ref))
            if i == 0:
                e2 = max(e2, abs(G_extremes(x, y, ZR) - ref) / abs(ref))
        print(f"    ZR={ZR:6}: closed(A,B) {float(e1):.1e}, direct S^4 / extremes form {float(e2):.1e}")
        worst = max(worst, e1, e2)
    report("N6 Green's function vs partial-wave sum, 7 couplings incl. ZR = -1/8, -0.3, -2 (complex s)", worst, 1e-10)
    # Z = 0: elementary form Om Om / (8 pi R sqrt((1-c+)(1-c-))) = 1/(4 pi |x-y|)
    e = 0
    for _ in range(10):
        x, y = rand_pair()
        r, rp = sph(x), sph(y); A, B = AB(x, y)
        elem = Om(r) * Om(rp) / (8 * np.pi * R * np.sqrt((1 - A - B) * (1 - A + B)))
        e = max(e, abs(elem * 4 * np.pi * np.linalg.norm(x - y) - 1),
                abs(G_closed(x, y, 0) * 4 * mp.pi * np.linalg.norm(x - y) - 1))
    report("N7 Z = 0: closed form = elementary form = 1/(4 pi |x - y|)", e, 1e-12)


# --- addition formula
def a_coef(k, j, lam=1.5):  # Gegenbauer addition theorem coefficient (DLMF 18.18.8)
    return 2 ** (2 * j) * mp.factorial(k - j) * mp.rf(lam, j) ** 2 / mp.rf(2 * lam, k + j) * (2 * lam + 2 * j - 1) / (2 * lam - 1)


def psi_l(N, l, r):
    k = 2 * l + 1; sg = np.log(r / R) / 2
    return r ** -0.5 * lpmv(k, N, np.tanh(sg))


def norm2(N, l):
    return quad(lambda t: R * lpmv(2 * l + 1, N, np.tanh(t)) ** 2 / np.cosh(t) ** 2, -40, 40, epsabs=0, epsrel=1e-12, limit=200)[0]


def Pi_sum(N, x, y):
    r, rp = sph(x), sph(y); cosTh = x @ y / (r * rp)
    return sum((2 * l + 1) / (4 * np.pi) * eval_legendre(l, cosTh) * psi_l(N, l, r) * psi_l(N, l, rp) / norm2(N, l)
               for l in range((N - 1) // 2 + 1))


def Pi_closed(N, x, y, direct=False):
    r, rp = sph(x), sph(y); A, B = AB(x, y)
    a = np.linspace(0, 2 * np.pi, 600, endpoint=False)  # exact for polynomials of degree < 600
    c = A + B * np.cos(a) if not direct else np.array([cosg_direct(x, y, t) for t in a])
    return (2 * N + 1) * Om(r) * Om(rp) / (16 * np.pi * R ** 2) * eval_gegenbauer(N - 1, 1.5, c).mean()


def Pi_partial(N, x, y):  # explicit partial waves: Gegenbauer addition theorem + Hopf identity
    k = N - 1; r, rp = sph(x), sph(y); cosTh = x @ y / (r * rp)
    ct, ctp = (r - R) / (r + R), (rp - R) / (rp + R); st, stp = np.sqrt(1 - ct ** 2), np.sqrt(1 - ctp ** 2)
    tot = 0.0
    for l in range(k // 2 + 1):
        tot += float(a_coef(k, 2 * l)) * (st * stp) ** (2 * l) * eval_gegenbauer(k - 2 * l, 2 * l + 1.5, ct) \
            * eval_gegenbauer(k - 2 * l, 2 * l + 1.5, ctp) * eval_legendre(l, cosTh)
    return (2 * N + 1) * Om(r) * Om(rp) / (16 * np.pi * R ** 2) * tot


def check_addition():
    e = ep = 0
    for N in range(1, 9):
        for i in range(6):
            x, y = rand_pair()
            ref = Pi_sum(N, x, y)
            e = max(e, abs(Pi_closed(N, x, y) - ref) / abs(ref))
            if i < 2:
                e = max(e, abs(Pi_closed(N, x, y, True) - ref) / abs(ref))
            ep = max(ep, abs(Pi_partial(N, x, y) - Pi_closed(N, x, y)) / abs(Pi_closed(N, x, y)))
    report("N8 addition formula Pi_N, N = 1..8, vs explicit eigenfunction sum", e, 1e-9)
    report("N9 fibre-averaged Gegenbauer = explicit partial waves (Hopf identity), N = 1..8", ep, 1e-11)
    # normalization: int_0^inf Om^2 S_l^2 w r^2 dr = 4 R^2 (2l+1) / ((2N+1) a_{k,2l})
    e = 0
    for N in range(1, 9):
        k = N - 1
        for l in range(k // 2 + 1):
            def f(r):
                ct = (r - R) / (r + R); st = np.sqrt(1 - ct ** 2)
                S = st ** (2 * l) * eval_gegenbauer(k - 2 * l, 2 * l + 1.5, ct)
                return Om(r) ** 2 * S ** 2 * 2 / (r * (1 + r / R) ** 2) * r ** 2
            # substitute r = R e^{2 sigma}
            I = quad(lambda sg: f(R * np.exp(2 * sg)) * 2 * R * np.exp(2 * sg), -40, 40, epsabs=0, epsrel=1e-12, limit=400)[0]
            ref = 4 * R ** 2 * (2 * l + 1) / ((2 * N + 1) * float(a_coef(k, 2 * l)))
            e = max(e, abs(I / ref - 1))
    report("N10 normalization of the radial functions read off the addition formula", e, 1e-9)
    # Hopf identity numerically for larger l
    e = 0
    from scipy.special import eval_chebyu
    a = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    for l in range(0, 21):
        for c in rng.uniform(-1, 1, 5):
            e = max(e, abs(eval_chebyu(2 * l, c * np.cos(a)).mean() - eval_legendre(l, 2 * c * c - 1)),
                    abs(eval_chebyu(2 * l + 1, c * np.cos(a)).mean()))
    report("N11 Hopf identity, l = 0..20, random c (trapezoid exact for polynomials)", e, 1e-11)


if __name__ == "__main__":
    print("part 1: symbolic"); symbolic()
    print("part 2: numerical")
    check_unit_kernel(); check_geometry(); check_green(); check_addition()
    print("ALL HELD" if all(ALL) else "SOME FAILED", f"({sum(ALL)}/{len(ALL)})")
