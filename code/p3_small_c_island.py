import numpy as np
from scipy.integrate import solve_ivp
Lam=3.0
def tip(tau,c):
    f=lambda t,s:[s[1],-c*s[1]+(s[0]+s[2])**2-1,(Lam-s[2])/tau]
    ev=lambda t,s:(s[0]+s[2])-4; ev.terminal=True
    ev2=lambda t,s:0.05-((s[0]+Lam+1)**2+s[1]**2+(s[2]-Lam)**2)**0.5; ev2.terminal=True
    s=solve_ivp(f,[0,60+20*tau+14/c],[-1,0,0],events=[ev,ev2],rtol=1e-11,atol=1e-13,method='LSODA')
    return len(s.t_events[0])>0 and len(s.t_events[0])>0 and s.t_events[0].size>0
def bis(a,b,c):
    sa=tip(a,c)
    for _ in range(34):
        m=.5*(a+b)
        if tip(m,c)==sa: a=m
        else: b=m
    return .5*(a+b)
for c in [0.02,0.015,0.01,0.0075,0.005]:
    taus=np.linspace(0.62,0.66,41); st=[tip(t,c) for t in taus]
    ch=[bis(taus[i],taus[i+1],c) for i in range(40) if st[i]!=st[i+1]]
    print(c,ch,flush=True)
