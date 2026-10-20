import numpy as np, sys
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
def tesc(tau,c,Lam):
    f=lambda t,y:[y[1],-c*y[1]+(y[0]+y[2])**2-1,(Lam-y[2])/tau]
    ev=lambda t,y:(y[0]+y[2])-3.0; ev.terminal=True
    sol=solve_ivp(f,[0,80+30*tau],[-1,0,0],events=ev,rtol=1e-13,atol=1e-15,method='DOP853')
    return sol.t_events[0][0] if len(sol.t_events[0]) else 1e3
mu=lambda c:(-c+np.sqrt(c*c+8))/2
cfg={3.0:[(0.0516,0.918,0.315),(0.0552,0.911,0.2585),(0.0588,0.904,0.187),(0.0624,0.897,0.06)],
     2.5:[(0.048,0.6606,0.2376),(0.050,0.658,0.19),(0.052,0.6557,0.1354),(0.053,0.6545,0.093)],
     4.0:[(0.052,1.463,0.253),(0.056,1.4507,0.1753),(0.058,1.446,0.116),(0.059,1.4435,0.0711)],
     5.0:[(0.044,2.0409,0.2269),(0.048,2.0243,0.138),(0.0485,2.022,0.122),(0.0495,2.019,0.081)]}
Lam=float(sys.argv[1]); cs=[];Y=[]
for c,cen,W in cfg[Lam]:
    h=max(0.3*W,0.012)
    r=minimize_scalar(lambda t:tesc(t,c,Lam),bounds=(cen-h,cen+h),method='bounded',options={'xatol':1e-9})
    cs.append(c);Y.append(np.exp(-mu(c)*r.fun)); print(Lam,c,"minT",r.fun,"tau",r.x,flush=True)
k,b=np.polyfit(cs[-3:],Y[-3:],1); print("Lam",Lam,"c* (extrap last3) =",-b/k,flush=True)
