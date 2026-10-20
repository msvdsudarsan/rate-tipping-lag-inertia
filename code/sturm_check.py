"""Check (i) tipping <=> w has a zero (300 random nonnegative inputs) and (ii) the Sturm sufficient criterion."""
import numpy as np
from scipy.integrate import solve_ivp
rng = np.random.default_rng(12345)
def tip_ode(qf, tmax, thr=6.0):
    f = lambda t, y: [y[0]**2-1+qf(t)]
    ev = lambda t, y: y[0]-thr; ev.terminal = True
    s = solve_ivp(f, [0, tmax], [-1.0], events=ev, rtol=1e-11, atol=1e-13, method='DOP853', max_step=0.05)
    return len(s.t_events[0]) > 0
def tip_w(qf, tmax):
    f = lambda t, y: [y[1], (1-qf(t))*y[0]]
    ev = lambda t, y: y[0]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, [0, tmax], [1.0, 1.0], events=ev, rtol=1e-11, atol=1e-13, method='DOP853', max_step=0.05)
    return len(s.t_events[0]) > 0
def make():
    k = rng.integers(1, 4); A = rng.uniform(0.3, 6, k); c = rng.uniform(0.2, 8, k); s = rng.uniform(0.1, 2, k)
    return (lambda t: float(np.sum(A*np.exp(-((t-c)/s)**2)))), (A, c, s)
agree = n = used = viol = 0
for i in range(300):
    q, (A, c, s) = make(); tmax = c.max()+6*s.max()+30
    a = tip_ode(q, tmax); b = tip_w(q, tmax); n += 1; agree += (a == b)
    ts = np.linspace(0, tmax, 1500); qs = np.array([q(t) for t in ts]); found = False
    for L in np.linspace(0.3, 6, 25):
        m = int(L/(ts[1]-ts[0]))
        if m < 2: continue
        thr = 1+np.pi**2/L**2
        for j in range(0, len(ts)-m, 5):
            if qs[j:j+m+1].min() > thr: found = True; break
        if found: break
    if found:
        used += 1
        if not a: viol += 1
print("tipping<=>w-zero agreement: %d/%d; Sturm criterion triggered %d times, violations %d" % (agree, n, used, viol))
open('../data/sturm_check.txt', 'w').write("agreement %d/%d; sturm triggered %d; violations %d\n" % (agree, n, used, viol))
