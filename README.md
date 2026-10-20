# Rate-induced tipping with lagged forcing and inertia: code and data

Code and numerical data only (the manuscript is not included).

Requirements: Python 3, numpy, scipy, mpmath, matplotlib. Run every script from the `code/` folder.

## Quick checks (exact formulas)
- `t1_check.py`: recomputes the exact thresholds of Theorems 1 and 2 (Tables I and II) and compares them with bisection on the ODE.
- `sturm_check.py`, `watson_check.py`, `ratio_check.py`: checks of the Sturm criterion, Watson's formula and the Bessel-zero ratios (Appendix B).

## Long runs (full-system integrations; run times were not measured, expect minutes to hours)
- `t2_check.py`, `t2_remainder_scaling.py`: inertial correction, Table III and the remainder coefficients.
- `p3_map.py`, `p3_small_c_island.py`, `p3_window_continuation.py`, `p3_cstar_method1.py`, `p3_cstar_method2.py`, `p3_edge_manifold.py`: weak-damping map, window width, critical damping c*.
- `swing_check.py`: swing-equation example.
- `level_check.py`: sensitivity to the escape level (3, 4, 6).

## Tables and figures
- `make_tables.py` writes the LaTeX tables from the saved results in `data/`.
- `make_figs_aip.py` (and `make_figs.py`, `make_graphical_abstract.py`) write the figures.

Outputs go to `../data` and `../figs`. `data/raw` holds the raw outputs behind the c* tables and the remainder coefficients. These scripts overwrite existing result files.
