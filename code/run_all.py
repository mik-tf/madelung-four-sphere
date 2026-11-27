#!/usr/bin/env python3
"""Run every check and print PASS/FAIL per prediction. Run from the repository root:
    OMP_NUM_THREADS=1 nice -n 19 python3 code/run_all.py
Needs Python 3, numpy, scipy, sympy, mpmath. Exit status is nonzero if anything fails.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, *args):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *args], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def lemma_status(out, n):
    """True if every line in the 'Lemma n:' block of proofs_check.py output starts with 'held'."""
    block, on = [], False
    for line in out.splitlines():
        m = re.match(r"Lemma (\d+):", line)
        if m:
            on = int(m.group(1)) == n
            continue
        if on and re.match(r"(held|FAILED)", line):
            block.append(line)
    return bool(block) and all(b.startswith("held") for b in block)


def maxerr(out, label):
    m = re.search(re.escape(label) + r"\s*([0-9.eE+-]+)", out)
    return float(m.group(1)) if m else float("inf")


rows = []
rc, out = run("proofs_check.py")
proofs_ok = {n: lemma_status(out, n) for n in range(1, 7)}
rows.append(("P1  couplings Z R = N(N+1)/2 (Lemma 3, sympy)", proofs_ok[3]))
rows.append(("P2a count by exact rank, k = 0..8 (Lemma 5, sympy)", proofs_ok[5]))
rows.append(("Lemma 1  KS lift (sympy)", proofs_ok[1]))
rows.append(("Lemma 2  conformal weight and covariance (sympy)", proofs_ok[2]))
rows.append(("P3b spindle Lemma 4, DO Eq. 25 symbolic", proofs_ok[4]))
rows.append(("P4  Kelvin = pole swap, Lemma 6 (sympy)", proofs_ok[6]))
rows.append(("proofs_check.py exit status", rc == 0))

rc, out = run("p2_enumeration.py")
rows.append(("P2  independent enumeration 1,1,4,4,9,9,16,16,25", rc == 0 and "P2: PASS" in out))
rc, out = run("p3_spindle.py")
rows.append(("P3  spindle spectrum vs DO Eq. 25 (numerical shooting)", rc == 0 and "P3: PASS" in out))
rc, out = run("p6_scalar_curvature.py")
rows.append(("P6  symbolic scalar curvature", rc == 0 and "P6: PASS" in out))
rc, out = run("payoff_green_check.py")
g = maxerr(out, "Green max rel error:")
a = maxerr(out, "Addition formula max rel error:")
e = maxerr(out, "max |direct extremes - (A -/+ B)|:")
rows.append((f"P5  Green's function (max rel err {g:.1e}), addition formula ({a:.1e})",
             rc == 0 and g < 1e-9 and a < 1e-9 and e < 1e-4))  # extremes are sampled on a 721-point circle grid

rc, out = run("theorem_d_check.py")
rows.append(("P5b Theorem D full checks (all real Z off the spectrum, Hopf identity, every N)",
             rc == 0 and "ALL HELD (24/24)" in out))

rc, out = run("p7_extra_checks.py")
rows.append(("P7  norm ratios 4/21, 4/105 (k = 2) and charge-1 sector sizes (sympy)", rc == 0 and "P7: PASS" in out))

rc, out = run("make_figures.py", "--check")
rows.append(("P8  figure layout check (no overlaps of text with text, lines, arrows, patches; inside canvas)",
             rc == 0 and out.count(": PASS") == 5))

rc, out = run("ds4_m1.py")
rows.append(("P9  de Sitter remark, N = 1..8 (exact sympy)", rc == 0 and "OVERALL PASS" in out))

width = max(len(r[0]) for r in rows)
print()
for name, ok in rows:
    print(f"{name.ljust(width)}  {'PASS' if ok else 'FAIL'}")
allok = all(ok for _, ok in rows)
print("\nALL PASS" if allok else "\nSOME FAILED")
sys.exit(0 if allok else 1)
