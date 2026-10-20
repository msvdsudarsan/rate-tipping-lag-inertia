import numpy as np, sys
from scipy.integrate import solve_ivp
def tip(tau,c,Lam):
    f=lambda t,s:[s[1],-c*s[1]+(s[0]+s[2])**2-1,(Lam-s[2])/tau]
    ev=lambda t,s:(s[0]+s[2])-4; ev.terminal=True
    ev2=lambda t,s:0.05-((s[0]+Lam+1)**2+s[1]**2+(s[2]-Lam)**2)**0.5; ev2.terminal=True
    s=solve_ivp(f,[0,60+20*tau+14/max(c,.01)],[-1,0,0],events=[ev,ev2],rtol=1e-10,atol=1e-12,method='LSODA')
    return len(s.t_events[0])>0
def bis(a,b,c,Lam,n=36):
    sa=tip(a,c,Lam)
    for _ in range(n):
        m=.5*(a+b)
        if tip(m,c,Lam)==sa: a=m
        else: b=m
    return .5*(a+b)
def window(c,Lam,lo,hi):
    # find N->T (a2) then T->N (b3) in [lo,hi]
    taus=np.linspace(lo,hi,81); st=[tip(t,c,Lam) for t in taus]
    a2=b3=None
    for i in range(80):
        if (not st[i]) and st[i+1] and a2 is None: a2=bis(taus[i],taus[i+1],c,Lam)
        if st[i] and (not st[i+1]) and a2 is not None and b3 is None: b3=bis(taus[i],taus[i+1],c,Lam)
    return a2,b3
Lam=float(sys.argv[1]); lo,hi=float(sys.argv[2]),float(sys.argv[3]); c0=float(sys.argv[4]); c1=float(sys.argv[5])
cs=np.linspace(c0,c1,10)
for c in cs:
    a2,b3=window(c,Lam,lo,hi)
    print(Lam,round(c,5),a2,b3,None if a2 is None or b3 is None else b3-a2,flush=True)
    if a2 is None or b3 is None: break
    w=0.75*(b3-a2)+0.03; mid=.5*(a2+b3); lo,hi=max(mid-w,0.1),mid+w
