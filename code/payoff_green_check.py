"""Numerical check of Theorem D of the paper: closed form of the Green's function of
-Lap psi = Z w_R psi on R^3, w_R = 2/(r (1 + r/R)^2) (Demkov-Ostrovsky mu = 1/2).

Independent side (no S^4, no KS):  radial problem per l, solved with Legendre
functions (sigma = ln(r/R)/2, kappa = 2l+1, s(s+1) = 2 Z R, psi = r^(-1/2) y(sigma),
y'' - kappa^2 y + s(s+1) sech^2 y = 0), summed over l with P_l(cos Theta).
Closed form side:  G = pi Om(u) Om(v) R^-1 (1/pi) int_0^pi Gunit(A + B cos a) da,
Gunit(x) = P_s'(-x)/(8 pi sin(pi s)) = d/dx[-P_s(-x)/(8 pi sin(pi s))].
Also: addition formula for the Madelung group N,
Pi_N(x,y) = (2N+1) Om(u) Om(v) / (16 pi R^2) * avg_a C^(3/2)_{N-1}(cos gamma_a),
against sum over the explicit eigenfunctions (n_r = N-1-2l >= 0).
Run: OMP_NUM_THREADS=1 nice -n 19 python3 payoff_green_check.py
"""
import numpy as np, mpmath as mp
from scipy.special import lpmv, eval_legendre, eval_gegenbauer
mp.mp.dps = 25
rng = np.random.default_rng(7)
R = 1.7

def sph(x):
    r = np.linalg.norm(x); return r, x[2] / r

def lift(x):  # a point u=(z1,z2) of the KS fibre over x (Hopf form of the paper)
    r = np.linalg.norm(x)
    z1 = np.sqrt((r + x[2]) / 2)
    z2 = np.conj((x[0] + 1j * x[1]) / 2) / z1
    return z1, z2

def S4(z1, z2):  # unit-sphere point of the stereographic map, |u|^2 = r
    n = abs(z1) ** 2 + abs(z2) ** 2
    s = np.sqrt(R) * 2 / (R + n)
    return np.array([s * z1.real, s * z1.imag, s * z2.real, s * z2.imag, (n - R) / (n + R)])

def rand_pair(sep=False):  # sep: radii differ by a factor >= 1.5 (radial l-sum converges like ratio^l)
    while True:
        x = rng.normal(size=3); y = rng.normal(size=3)
        x *= np.exp(rng.uniform(-1, 1)) * R / np.linalg.norm(x)
        y *= np.exp(rng.uniform(-1, 1)) * R / np.linalg.norm(y)
        if not sep or abs(np.log(np.linalg.norm(x) / np.linalg.norm(y))) > 0.4:
            return x, y

def Om(r): return 2 / (1 + r / R)

def cosg_direct(x, y, a):
    z1, z2 = lift(x); w1, w2 = lift(y); e = np.exp(1j * float(a))
    return float(S4(z1, z2) @ S4(e * w1, e * w2))

def AB(x, y):
    r, c = sph(x); rp, cp = sph(y)
    ct, ctp = (r - R) / (r + R), (rp - R) / (rp + R)
    st, stp = np.sqrt(1 - ct ** 2), np.sqrt(1 - ctp ** 2)
    cosTh = x @ y / (r * rp)
    return ct * ctp, st * stp * np.sqrt((1 + cosTh) / 2)

# ---------------- independent: radial sum
def G_radial(x, y, ZR, lmax=90):
    s = (-1 + mp.sqrt(1 + 8 * ZR)) / 2
    r, _ = sph(x); rp, _ = sph(y); cosTh = x @ y / (r * rp)
    sg = lambda rr: mp.log(rr / R) / 2
    tot = mp.mpf(0)
    for l in range(lmax):
        k = 2 * l + 1
        y1 = lambda t: mp.legenp(s, -k, -mp.tanh(t), type=2)  # regular at r=0
        y2 = lambda t: mp.legenp(s, -k, mp.tanh(t), type=2)   # decays at infinity
        t0 = mp.mpf(0.3)
        W = (mp.diff(y1, t0) * y2(t0) - y1(t0) * mp.diff(y2, t0)) / 2  # = r^2 (f1' f2 - f1 f2')
        lo, hi = (r, rp) if r < rp else (rp, r)
        g = y1(sg(lo)) * y2(sg(hi)) / mp.sqrt(lo * hi) / W
        term = (2 * l + 1) / (4 * mp.pi) * eval_legendre(l, cosTh) * g
        tot += term
        if l > 8 and abs(term) < mp.mpf(10) ** -18 * abs(tot): break
    return tot

# ---------------- closed form (Theorem D of the paper)
def Gu(xc, s):
    # d/dx[A P_s(-x)] = -A P_s'(-x), A = -1/(8 pi sin(pi s)); P_s' = derivative in z
    return (1 / (8 * mp.pi * mp.sin(mp.pi * s))) * mp.diff(lambda z: mp.legenp(s, 0, z, type=3), -xc)

def G_closed(x, y, ZR, direct=False):
    s = (-1 + mp.sqrt(1 + 8 * ZR)) / 2
    r, _ = sph(x); rp, _ = sph(y)
    A, B = AB(x, y)
    f = (lambda a: Gu(cosg_direct(x, y, a), s)) if direct else (lambda a: Gu(A + B * mp.cos(a), s))
    avg = (mp.quad(f, [0, mp.pi / 2, mp.pi]) / mp.pi) if not direct else (mp.quad(f, np.linspace(0, 2 * np.pi, 9)) / (2 * mp.pi))
    return mp.pi * Om(r) * Om(rp) / R * avg   # G_3 = pi Om Om G_g, G_g = G_unit / R

# ---------------- addition formula
def psi_l(N, l, r):
    k = 2 * l + 1; sg = np.log(r / R) / 2
    return r ** -0.5 * lpmv(k, N, np.tanh(sg))

def norm2(N, l):
    from scipy.integrate import quad
    return quad(lambda t: R * lpmv(2 * l + 1, N, np.tanh(t)) ** 2 / np.cosh(t) ** 2, -40, 40, epsabs=0, epsrel=1e-12)[0]

def Pi_sum(N, x, y):
    r, _ = sph(x); rp, _ = sph(y); cosTh = x @ y / (r * rp)
    return sum((2 * l + 1) / (4 * np.pi) * eval_legendre(l, cosTh) * psi_l(N, l, r) * psi_l(N, l, rp) / norm2(N, l)
               for l in range((N - 1) // 2 + 1))

def Pi_closed(N, x, y, direct=False):
    r, _ = sph(x); rp, _ = sph(y); A, B = AB(x, y)
    a = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    c = A + B * np.cos(a) if not direct else np.array([cosg_direct(x, y, t) for t in a])
    avg = eval_gegenbauer(N - 1, 1.5, c).mean()
    return (2 * N + 1) * Om(r) * Om(rp) / (16 * np.pi * R ** 2) * avg

if __name__ == "__main__":
    # geometry: over the fibre circle, direct cos(gamma) ranges exactly over [A-B, A+B]
    g = 0
    for _ in range(50):
        x, y = rand_pair(); A, B = AB(x, y)
        c = [cosg_direct(x, y, t) for t in np.linspace(0, 2 * np.pi, 721)]
        g = max(g, abs(max(c) - (A + B)), abs(min(c) + 0 - (A - B)))
    print("max |direct extremes - (A -/+ B)|:", g)
    errs = {}
    for ZR in (0.3, 1.7, 4.2, -0.05):
        e1, e2 = 0, 0
        for i in range(6):
            x, y = rand_pair(sep=True)
            ref = G_radial(x, y, ZR)
            cl = G_closed(x, y, ZR); cd = G_closed(x, y, ZR, direct=True)
            e1 = max(e1, abs(cl - ref) / abs(ref)); e2 = max(e2, abs(cd - ref) / abs(ref))
            if i == 0: print(f"  ZR={ZR} sample ref={mp.nstr(ref, 12)} closed={mp.nstr(cl, 12)}")
        print(f"ZR={ZR}: max rel err closed(A,B) {float(e1):.2e}, closed(direct S^4 points) {float(e2):.2e}")
        errs[ZR] = e1
    print("Green max rel error:", max(float(v) for v in errs.values()))
    ea = 0
    for N in (1, 2, 3, 4):
        for i in range(8):
            x, y = rand_pair()
            ref = Pi_sum(N, x, y); ea = max(ea, abs(Pi_closed(N, x, y) - ref) / abs(ref), abs(Pi_closed(N, x, y, True) - ref) / abs(ref))
        print(f"N={N}: running max rel err {ea:.2e}")
    print("Addition formula max rel error:", ea)
