import numpy as np
from scipy.integrate import solve_ivp
def tipr(th,eps,Lam):
    f=lambda s,y:[(1-2*eps*y[0])*(y[0]**2-1+(Lam/th)*np.exp(-s/th))]
    ev=lambda s,y:y[0]-3; ev.terminal=True
    sol=solve_ivp(f,[0,60+30*th],[-1.0],events=ev,rtol=1e-13,atol=1e-15,method='DOP853')
    return len(sol.t_events[0])>0
def bis(fn,lo,hi,n=45):
    for _ in range(n):
        m=.5*(lo+hi)
        if fn(m): lo=m
        else: hi=m
    return .5*(lo+hi)
for Lam,th0,k in [(3.0,0.56672992808,0.25387),(4.0,1.22371235414,0.41951)]:
    prev=None
    for c in [5,10,20,40,80]:
        tr=bis(lambda t:tipr(t,1/c**2,Lam),0.8*th0,1.4*th0)
        r=(tr-th0)*c*c-k
        print(Lam,c,tr,"resid*c",r*c,"resid*c^2",r*c*c,flush=True)
