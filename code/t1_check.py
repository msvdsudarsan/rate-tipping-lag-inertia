import numpy as np, csv
from common import *
rows = []
for tau in [0.1, 0.25, 0.5, 1.0, 2.0, 4.0]:
    b = Lam_c_bessel(tau); o = Lam_c_ode(tau)
    rows.append((tau, b, o, abs(b-o)))
    print("tau %.2f  Bessel %.10f  ODE %.10f  diff %.1e" % (tau, b, o, abs(b-o)), flush=True)
with open('../data/table_T1_filtered.csv', 'w') as f:
    w = csv.writer(f); w.writerow(['tau', 'Lambda_c_Bessel', 'Lambda_c_ODE', 'abs_diff']); w.writerows(rows)
rows = []
for Q in [1.5, 2, 3, 5, 10, 100]:
    a = Lam_c_ramp(Q)
    lo, hi = 2.0, a*3
    for _ in range(50):
        m = .5*(lo+hi)
        if tips_ramp(Q, m): hi = m
        else: lo = m
    o = .5*(lo+hi)
    rows.append((Q, a, o, abs(a-o)))
    print("Q %6.1f  closed %.10f  ODE %.10f  diff %.1e" % (Q, a, o, abs(a-o)), flush=True)
with open('../data/table_T1prime_ramp.csv', 'w') as f:
    w = csv.writer(f); w.writerow(['Q', 'Lambda_c_closed', 'Lambda_c_ODE', 'abs_diff']); w.writerows(rows)
