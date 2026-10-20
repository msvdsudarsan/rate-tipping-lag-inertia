import numpy as np, sys
from scipy.integrate import solve_ivp
from scipy.optimize import root
Lam=3.0
def pt(logrho,tau,c,sg):
    s=(-c-np.sqrt(c*c+8))/2; T=26*tau
    xs=1-Lam*(1-np.exp(-T/tau))
    def f(t,y):
        u=y[0]+Lam*(1-np.exp(-t/tau)); return [y[1],-c*y[1]+u*u-1]
    d=sg*np.exp(logrho-abs(s)*T)
    sol=solve_ivp(f,[T,0],[xs+d,d*s],rtol=1e-12,atol=1e-15,method='DOP853')
    return sol.y[:,-1]
def guess(tau,c):
    best=(9,None,None)
    for sg in (1,-1):
        for lr in np.linspace(-3,9,1200):
            p=pt(lr,tau,c,sg); d=np.hypot(p[0]+1,p[1])
            if np.isfinite(d) and d<best[0]: best=(d,lr,sg)
    return best
c=float(sys.argv[1])
for tau0 in [float(a) for a in sys.argv[2:]]:
    d,lr,sg=guess(tau0,c)
    r=root(lambda z:pt(z[0],z[1],c,sg)+np.array([1.0,0.0]),[lr,tau0],tol=1e-13)
    print("c",c,"guess tau",tau0,"closest d",round(d,4),"-> solved tau",r.x[1],"logrho",r.x[0],"resid",np.abs(r.fun).max(),flush=True)
