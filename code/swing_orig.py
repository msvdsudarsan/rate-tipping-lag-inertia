import numpy as np, json
from scipy.integrate import solve_ivp
P0,P1=0.2,0.9
ds=np.pi-np.arcsin(P1); Vs=-P1*ds-np.cos(ds)
def run(r,tau,c):
    T0=10/r+10*tau+12
    d0=np.arcsin(P0)
    Pin=lambda t:P0+(P1-P0)*0.5*(1+np.tanh(r*t))
    def f(t,s):
        d,w,P=s; return [w,-c*w+P-np.sin(d),(Pin(t)-P)/tau]
    ev=lambda t,s:s[0]-np.pi; ev.terminal=True
    # phase 1: until ramp essentially done, plus filter settled
    Tend=10/r+25*tau+12
    s=solve_ivp(f,[-T0,Tend],[d0,0.0,P0],events=ev,rtol=1e-10,atol=1e-12,method='LSODA')
    if len(s.t_events[0])>0: return 1
    d,w,P=s.y[:,-1]
    E=w*w/2-P1*d-np.cos(d)
    # frozen final system: no escape if E < Vs (saddle potential energy) with delta<ds
    if d<ds and E<Vs-1e-6: return 0
    # else continue
    s2=solve_ivp(lambda t,y:[y[1],-c*y[1]+P1-np.sin(y[0])],[Tend,Tend+60/c],[d,w],events=ev,rtol=1e-10,atol=1e-12,method='LSODA')
    return 1 if len(s2.t_events[0])>0 else 0
rg=np.geomspace(0.03,40,70)
out={}
for c in [0.5,0.15,0.05]:
    for tau in [0,0.25,0.5,1,2,4]:
        if tau==0:
            tt=1e-4
        else: tt=tau
        st=[run(r,tt,c) for r in rg]
        ch=[(round(float(np.sqrt(rg[i]*rg[i+1])),3),'N>T' if st[i]==0 else 'T>N') for i in range(len(rg)-1) if st[i]!=st[i+1]]
        print("c",c,"tau",tau,"tip frac",np.mean(st),ch,flush=True)
