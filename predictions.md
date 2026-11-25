# Predictions

Each prediction is a line `P<n>:` and is tested by the script named in brackets. Run all with `python3 code/run_all.py`.

P1: the couplings are Z R = N(N + 1)/2, N = 1, 2, 3, ... [code/proofs_check.py, Lemma 3]
P2: the multiplicities are floor((N + 1)/2)^2 = 1, 1, 4, 4, 9, 9, 16, 16, 25, by independent enumeration of circle-invariant harmonics on R^5 [code/p2_enumeration.py]
P3: the spindle reproduces Demkov-Ostrovsky Eq. 25, v_N = R^2 mu^2 (N + 1/(2 mu)) (N + 1/(2 mu) - 1), for mu = 1/3, 1/2, 2/3, 1, 2 [code/p3_spindle.py, code/proofs_check.py Lemma 4]
P4: the Kelvin transform equals the pole swap of the 4-sphere [code/proofs_check.py, Lemma 6]
P5: the zero-energy Green's function and the addition formula are as in Theorem D [code/payoff_green_check.py]
P5b: Theorem D holds for every real Z off the spectrum (including Z R < -1/8 with complex s, and Z = 0), the Hopf identity holds for all l, and the addition formula holds for every N [code/theorem_d_check.py, 24 checks]
P6: there is no point symmetry beyond O(3) x Z_2 for mu other than 1; the scalar curvature is Scal = 4 (2 (1 - mu^2) cosh^2(mu x) + 6 mu^2) [code/p6_scalar_curvature.py]
P7: for k = 2 the ratio of the L^2(S^4) norm squared to 2 pi^2 times the Fischer norm squared of the Cauchy datum is 4/21 (l = 0) and 4/105 (l = 1); the circle-charge-1 sector sizes for k = 0..6 are 0, 2, 2, 6, 6, 12, 12 [code/p7_extra_checks.py]
P8: in each of the five figures no text overlaps another text, a line, an arrow, a marker or a patch outline, and every text lies inside the canvas [code/make_figures.py --check]

Falsifiers: a source shows the S^4 lift before this work; an eigenfunction without smooth lift or with extra spectrum; an independent count that differs from floor((N + 1)/2)^2; a spindle spectrum that differs from Eq. 25; a Green's function formula that fails at a new point.
