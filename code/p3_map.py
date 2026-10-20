"""Tipping map of the weakly damped 3-D system in the (c, tau) plane, Lambda = 3.
Tip: u=x+lambda reaches 4. No tip: state within 0.05 of the final stable state."""
import numpy as np, csv, sys
from multiprocessing import Pool
from scipy.integrate import solve_ivp
Lam = 3.0
def tip(tau, c, level=4.0, Lam=Lam):
    f = lambda t, s: [s[1], -c*s[1]+(s[0]+s[2])**2-1, (Lam-s[2])/tau]
    ev = lambda t, s: (s[0]+s[2])-level; ev.terminal = True
    ev2 = lambda t, s: 0.05-((s[0]+Lam+1)**2+s[1]**2+(s[2]-Lam)**2)**0.5; ev2.terminal = True
    s = solve_ivp(f, [0, 60+20*tau+14/max(c, .01)], [-1, 0, 0], events=[ev, ev2], rtol=1e-10, atol=1e-12, method='LSODA')
    return len(s.t_events[0]) > 0
def row(c):
    taus = np.linspace(0.30, 1.60, 131)
    st = np.array([tip(t, c) for t in taus])
    edges = []
    for i in range(len(taus)-1):
        if st[i] != st[i+1]:
            lo, hi = taus[i], taus[i+1]; s0 = st[i]
            for _ in range(30):
                m = .5*(lo+hi)
                if tip(m, c) == s0: lo = m
                else: hi = m
            edges.append(.5*(lo+hi))
    return c, st.astype(int).tolist(), edges
if __name__ == '__main__':
    cs = np.round(np.linspace(0.004, 0.076, 37), 5)
    with Pool(2) as p: res = p.map(row, cs)
    with open('../data/p3_map.csv', 'w') as f:
        w = csv.writer(f); w.writerow(['c'] + ['tau_%.4f' % t for t in np.linspace(0.30, 1.60, 131)])
        for c, st, e in res: w.writerow([c] + st)
    with open('../data/p3_edges.csv', 'w') as f:
        w = csv.writer(f); w.writerow(['c', 'edges...'])
        for c, st, e in res: w.writerow([c] + ['%.6f' % x for x in e])
    for c, st, e in res: print(c, ['%.5f' % x for x in e], flush=True)
