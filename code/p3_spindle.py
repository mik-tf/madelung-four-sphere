#!/usr/bin/env python3
"""P3: the spindle spectrum against Demkov-Ostrovsky (DO) Eq. 25.

DO Eq. 25: v_N = R^2 mu^2 (N + 1/(2 mu)) (N + 1/(2 mu) - 1), with
N = n_r + l + 1 + (1/mu - 1) l, for the family V_mu = -v / (r^2 [(r/R)^mu + (R/r)^mu]^2).

Independent side: numerical shooting, no closed form used. With r = R e^x,
psi = r^(-1/2) y and sigma = mu x, the sector l gives
    -y'' + kappa^2 y = E sech^2(sigma) y,   kappa = (2 l + 1)/(2 mu),   E = v_N / (R^2 mu^2).
We integrate the decaying solutions from both ends, find the values of E where they
match at sigma = 0 (sign changes of y y' / (y^2 + y'^2)), and compare the lowest
three with E = v_N / (R^2 mu^2) from Eq. 25 (mu in 1/3, 1/2, 2/3, 1, 2; l = 0, 1).
Prints PASS or FAIL.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

L = 40.0


def match(E, kappa):
    def rhs(s, u):
        return [u[1], (kappa**2 - E / np.cosh(s)**2) * u[0]]
    sol = solve_ivp(rhs, [-L, 0.0], [1.0, kappa], rtol=1e-12, atol=1e-300, method="DOP853")
    y, yp = sol.y[0, -1], sol.y[1, -1]
    return y * yp / (y * y + yp * yp)  # zero for even (y' = 0) and odd (y = 0) states


def lowest_roots(kappa, count, emax):
    grid = np.linspace(1e-6, emax, 600)
    vals = [match(e, kappa) for e in grid]
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(match, grid[i], grid[i + 1], args=(kappa,), xtol=1e-12))
            if len(roots) == count:
                break
    return roots


def do_eq25(mu, l, nr, R=1.0):
    N = nr + l + 1 + (1.0 / mu - 1.0) * l
    return R**2 * mu**2 * (N + 1 / (2 * mu)) * (N + 1 / (2 * mu) - 1)


def main():
    worst = 0.0
    for mu in (1 / 3, 1 / 2, 2 / 3, 1.0, 2.0):
        for l in (0, 1):
            kappa = (2 * l + 1) / (2 * mu)
            want = [do_eq25(mu, l, nr) / mu**2 for nr in range(3)]
            got = lowest_roots(kappa, 3, want[-1] * 1.05)
            if len(got) != 3:
                print(f"mu={mu:.4f} l={l}: found {len(got)} roots, FAIL")
                return 1
            err = max(abs(g - w) / w for g, w in zip(got, want))
            worst = max(worst, err)
            print(f"mu={mu:.4f} l={l}: E = {[round(g, 8) for g in got]}  max rel diff {err:.1e}")
    ok = worst < 1e-7
    print(f"worst relative difference {worst:.1e}")
    print("P3:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
