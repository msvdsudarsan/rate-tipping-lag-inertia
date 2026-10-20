"""Swing equation with lagged tanh ramp: critical ramp rate r_c(tau, c), refined by bisection."""
import numpy as np, csv
from scipy.integrate import solve_ivp
P0, P1 = 0.2, 0.9
ds = np.pi-np.arcsin(P1); Vs = -P1*ds-np.cos(ds)
def run(r, tau, c):
    T0 = 10/r+10*tau+12; d0 = np.arcsin(P0)
    Pin = lambda t: P0+(P1-P0)*0.5*(1+np.tanh(r*t))
    f = lambda t, s: [s[1], -c*s[1]+s[2]-np.sin(s[0]), (Pin(t)-s[2])/tau]
    ev = lambda t, s: s[0]-np.pi; ev.terminal = True
    Tend = 10/r+25*tau+12
    s = solve_ivp(f, [-T0, Tend], [d0, 0.0, P0], events=ev, rtol=1e-10, atol=1e-12, method='LSODA')
    if len(s.t_events[0]) > 0: return 1
    d, w, P = s.y[:, -1]
    E = w*w/2-P1*d-np.cos(d)
    if d < ds and E < Vs-1e-6: return 0
    s2 = solve_ivp(lambda t, y: [y[1], -c*y[1]+P1-np.sin(y[0])], [Tend, Tend+60/c], [d, w], events=ev, rtol=1e-10, atol=1e-12, method='LSODA')
    return 1 if len(s2.t_events[0]) > 0 else 0
rg = np.geomspace(0.03, 40, 70)
rows = []
for c in [0.5, 0.15, 0.05]:
    for tau in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]:
        tt = 1e-4 if tau == 0 else tau
        st = [run(r, tt, c) for r in rg]
        roots = []
        for i in range(len(rg)-1):
            if st[i] != st[i+1]:
                lo, hi = rg[i], rg[i+1]; slo = st[i]
                for _ in range(24):
                    m = np.sqrt(lo*hi)
                    if run(m, tt, c) == slo: lo = m
                    else: hi = m
                roots.append((round(float(np.sqrt(lo*hi)), 4), 'N>T' if st[i] == 0 else 'T>N'))
        print("c", c, "tau", tau, "fraction tipping", round(float(np.mean(st)), 3), roots, flush=True)
        rows.append((c, tau, round(float(np.mean(st)), 3), ';'.join('%s@%s' % (b, a) for a, b in roots)))
with open('../data/table_swing.csv', 'w') as f:
    w = csv.writer(f); w.writerow(['c', 'tau', 'frac_tip_over_rates', 'transitions']); w.writerows(rows)
