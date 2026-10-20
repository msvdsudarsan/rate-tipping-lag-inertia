import numpy as np
from scipy.special import jv
from scipy.optimize import brentq
def zeros(nu,K=4):
    xmax=nu+ (K+1)*np.pi+ 6*max(nu,1)**(1/3)+8
    x=np.linspace(1e-6,xmax,12000); y=jv(nu,x); z=[]
    for i in range(len(x)-1):
        if y[i]*y[i+1]<0: z.append(brentq(lambda t:jv(nu,t),x[i],x[i+1],xtol=1e-13))
        if len(z)==K: break
    return z
nus=np.concatenate([np.linspace(-0.97,3,120),np.linspace(3.1,40,140)])
R=[];ok=True
prev=None
for n in nus:
    z=zeros(n,4); R.append([z[0]/z[k] for k in (1,2,3)])
R=np.array(R); d=np.diff(R,axis=0)
print("ratios j1/jk increasing in nu for k=2,3,4:",(d>0).all(axis=0), "min diffs",d.min(axis=0))
# Rayleigh sum check
n=1.7; z=zeros(n,60); print("Rayleigh sum partial",sum(1/t**2 for t in z),"vs",1/(4*(n+1)))
