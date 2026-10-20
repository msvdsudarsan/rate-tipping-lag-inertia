"""Inertial correction: theta*(c)=tau*(c)/c = theta0 + kappa/c^2 + O(c^-4).
theta0: root of Lambda_c(theta)=Lambda (Bessel); kappa: adjoint (Melnikov-type) formula;
reduced = 1-D slow-manifold model; full = 3-D system x'=v, v'=-cv+(x+lam)^2-1, tau*lam'=Lambda-lam."""
import numpy as np, csv, sys
from multiprocessing import Pool
from scipy.special import jv, jvp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from common import Lam_c_bessel

def theta0(Lam): return brentq(lambda t: Lam_c_bessel(t)-Lam, 0.05, Lam, xtol=1e-14)

def kappa(Lam, th):
    a = 2*np.sqrt(Lam*th)
    J  = lambda s: jv(2*th, a*np.exp(-s/(2*th)))
    dJ = lambda s: jvp(2*th, a*np.exp(-s/(2*th)))*a*np.exp(-s/(2*th))*(-1/(2*th))
    J0 = J(0.0); w = lambda s: J(s)/J0; wp = lambda s: dJ(s)/J0
    u0 = lambda s: -wp(s)/w(s)
    qd = lambda s: (Lam/th)*np.exp(-s/th)
    num = quad(lambda s: w(s)**2*2*u0(s)*(u0(s)**2-1+qd(s)), 0, 60*th, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    den = quad(lambda s: w(s)**2*Lam*th**-2*np.exp(-s/th)*(s/th-1), 0, 60*th, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return num/den, den

def tip_red(th, eps, Lam):
    f = lambda s, y: [(1-2*eps*y[0])*(y[0]**2-1+(Lam/th)*np.exp(-s/th))]
    ev = lambda s, y: y[0]-3; ev.terminal = True
    sol = solve_ivp(f, [0, 60+30*th], [-1.0], events=ev, rtol=1e-13, atol=1e-15, method='DOP853')
    return len(sol.t_events[0]) > 0

def tip_full(tau, c, Lam):
    f = lambda t, s: [s[1], -c*s[1]+(s[0]+s[2])**2-1, (Lam-s[2])/tau]
    ev = lambda t, s: (s[0]+s[2])-4; ev.terminal = True
    s = solve_ivp(f, [0, 30*c+20*tau], [-1, 0, 0], events=ev, rtol=1e-11, atol=1e-13, method='LSODA')
    return len(s.t_events[0]) > 0

def bis(fn, lo, hi, n=40):
    for _ in range(n):
        m = .5*(lo+hi)
        if fn(m): lo = m     # tipping below threshold, safe above
        else: hi = m
    return .5*(lo+hi)

def job(args):
    Lam, c, th0 = args
    red = bis(lambda t: tip_red(t, 1/c**2, Lam), 0.8*th0, 1.4*th0)
    full = bis(lambda t: tip_full(t*c, c, Lam), 0.8*th0, 1.4*th0, n=36)
    return Lam, c, red, full

if __name__ == '__main__':
    out = []; jobs = []; info = {}
    for Lam in [2.5, 3.0, 4.0]:
        th0 = theta0(Lam); k, den = kappa(Lam, th0); info[Lam] = (th0, k, den)
        print("Lambda", Lam, "theta0 %.11f kappa %.6f adjoint-denominator %.4f" % (th0, k, den), flush=True)
        for c in [10, 30, 100]: jobs.append((Lam, c, th0))
    with Pool(9) as p:
        for Lam, c, red, full in p.imap_unordered(job, jobs):
            th0, k, den = info[Lam]
            exp1 = th0 + k/c**2
            print("Lam %.1f c %3d  full %.8f  reduced %.8f  th0+k/c^2 %.8f  (full-th0)c^2 %.5f  (red-th0)c^2 %.5f" %
                  (Lam, c, full, red, exp1, (full-th0)*c*c, (red-th0)*c*c), flush=True)
            out.append((Lam, c, th0, k, den, full, red, exp1))
    out.sort()
    with open('../data/table_T2_inertia.csv', 'w') as f:
        w = csv.writer(f); w.writerow(['Lambda', 'c', 'theta0', 'kappa', 'adj_den', 'theta_full3D', 'theta_reduced', 'theta0+kappa/c2']); w.writerows(out)
