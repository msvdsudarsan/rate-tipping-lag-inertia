"""Shared helpers: first positive zero of J_nu (nu>-1), Riccati ODE threshold."""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
mp.mp.dps = 30

def j1(nu):
    """first positive zero of J_nu, nu>-1 (mpmath has no nu<0 support, so scan + findroot)."""
    f = lambda x: mp.besselj(nu, x)
    x = mp.mpf('0.02'); fx = f(x)
    while True:
        y = x + mp.mpf('0.05'); fy = f(y)
        if fx*fy < 0:
            return mp.findroot(f, (x, y), solver='anderson')
        x, fx = y, fy

def Lam_c_bessel(tau):
    nu = 2*mp.mpf(tau) - 1
    return float(j1(nu)**2/(4*tau))

def tips_filtered(tau, Lam, rtol=1e-12):
    """u'=u^2-1+(Lam/tau)exp(-t/tau), u(0)=-1. Tipping <=> u reaches 1 (then u'>0 forever)."""
    f = lambda t, u: [u[0]**2 - 1 + (Lam/tau)*np.exp(-t/tau)]
    ev = lambda t, u: u[0] - 1.0
    ev.terminal = True
    s = solve_ivp(f, [0, 60*tau + 60], [-1.0], events=ev, rtol=rtol, atol=1e-14, method='DOP853')
    return len(s.t_events[0]) > 0

def Lam_c_ode(tau, lo=2.0, hi=None, n=50):
    hi = hi or 4*max(tau, 1)+10
    for _ in range(n):
        m = 0.5*(lo+hi)
        if tips_filtered(tau, m): hi = m
        else: lo = m
    return 0.5*(lo+hi)

def Lam_c_ramp(Q):
    """ramp-and-stop: lambda=Q t up to Lambda; critical Lambda."""
    if Q <= 1: return np.inf
    w = np.sqrt(Q-1)
    return 2*Q*np.arctan(1/w)/w

def tips_ramp(Q, Lam):
    T = Lam/Q
    f = lambda t, u: [u[0]**2 - 1 + (Q if t < T else 0.0)]
    ev = lambda t, u: u[0] - 1.0
    ev.terminal = True
    s = solve_ivp(f, [0, T], [-1.0], events=ev, rtol=1e-13, atol=1e-15, method='DOP853', max_step=T/200)
    if len(s.t_events[0]) > 0: return True
    # after the ramp q=0: u'=u^2-1 -> tips iff u(T)>=1 ; edge u=1 is not tipping (strict)
    return s.y[0, -1] > 1.0
