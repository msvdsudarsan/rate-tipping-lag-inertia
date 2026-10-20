"""Sensitivity of the tipping classification to the escape level u=3,4,6 (and proximity-to-sink rule)."""
import numpy as np
from p3_map import tip
mism = tot = 0
for c in [0.03, 0.05, 0.08]:
    for tau in np.linspace(0.3, 1.5, 25):
        a, b, d = tip(tau, c, 3.0), tip(tau, c, 4.0), tip(tau, c, 6.0)
        tot += 1; mism += (not (a == b == d))
        if not (a == b == d): print("mismatch", c, tau, a, b, d)
print("points", tot, "mismatches", mism)
open('../data/level_check.txt', 'w').write("points %d mismatches %d\n" % (tot, mism))
