# The Demkov-Ostrovsky Madelung model as the circle-invariant conformal Laplacian of the round 4-sphere

Author: Michaël Perreault

Version 1.0.1.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24522163.svg)](https://doi.org/10.5281/zenodo.24522163)

The zero-energy Demkov-Ostrovsky problem for the Tietz field (the n + l Madelung grouping, mu = 1/2) is, through the Kustaanheimo-Stiefel map and stereographic projection, the circle-invariant conformal Laplacian of the round 4-sphere. This explains the group sizes and the Kelvin inversion, shows that the Demkov-Ostrovsky family is a family of spindles that is smooth after the lift exactly at mu = 1/2, and gives the zero-energy Green's function in closed form.

- **A.** The lift to the 4-sphere is a bijection onto circle-invariant solutions; the eigenvalues are Z_N = N(N + 1)/(2R); group N has floor((N + 1)/2)^2 states; the Kelvin transform is the pole swap.
- **B.** The Demkov-Ostrovsky family is a spindle with cone parameter 2 mu; it is smooth iff mu = 1/2.
- **C.** For mu other than 1 there are no point realisations of the finite transformations beyond O(3) x Z_2.
- **D.** The zero-energy Green's function is G_Z(x, y) = (pi/R) Omega(u) Omega(v) times the Kustaanheimo-Stiefel fibre average of the S^4 resolvent kernel (checked numerically).

Paper: source in [paper.md](paper.md) (full text). Predictions: [predictions.md](predictions.md).

## Licences

- Text (`paper.md`) and archived outputs: CC BY 4.0, see [LICENSE-paper.md](LICENSE-paper.md).
- Code in `code/`: MIT, see [LICENSE](LICENSE).

## Re-run the checks

Python 3 with numpy, scipy, sympy and mpmath. From the repository root:

```
OMP_NUM_THREADS=1 nice -n 19 python3 code/run_all.py
```

It runs every script below and prints PASS or FAIL per prediction (about 3 minutes on one core). Each script also runs alone from the repository root.

| Check | Command | What it checks |
|---|---|---|
| Lemmas 1 to 6 and Theorem B | `python3 code/proofs_check.py` | symbolic checks of the lemmas |
| Circle-invariant harmonic count | `python3 code/p2_enumeration.py` | independent count of circle-invariant harmonics on R^5 |
| Spindle spectrum | `python3 code/p3_spindle.py` | numerical spindle spectrum against Demkov-Ostrovsky Eq. 25 |
| Green's function | `python3 code/payoff_green_check.py` | Green's function and addition formula |
| Scalar curvature | `python3 code/p6_scalar_curvature.py` | symbolic scalar curvature of g_mu |
| Figure layout | `python3 code/make_figures.py --check` | figure layout check (no overlaps, inside the canvas) |

## Cite as

Perreault, M. The Demkov-Ostrovsky Madelung model as the circle-invariant conformal Laplacian of the round 4-sphere, version 1.0.1. https://github.com/mik-tf/madelung-four-sphere (DOI: [10.5281/zenodo.24522163](https://doi.org/10.5281/zenodo.24522163) (all versions; v1.0 is [10.5281/zenodo.24522164](https://doi.org/10.5281/zenodo.24522164))). See [CITATION.cff](CITATION.cff).
