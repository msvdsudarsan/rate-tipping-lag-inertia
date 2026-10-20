import numpy as np, sys
exec(open('p3_window_continuation.py').read().split('Lam=float')[0])
Lam=float(sys.argv[1]); c=float(sys.argv[2]); dc=float(sys.argv[3])
a2,b3=window(c,Lam,0.05,1.3*Lam)
print("start",Lam,c,a2,b3,flush=True)
pts=[]
while a2 is not None and b3 is not None:
    pts.append((c,b3-a2)); 
    w=0.75*(b3-a2)+0.03; mid=.5*(a2+b3); lo,hi=max(mid-w,0.05),mid+w
    c+=dc
    a2,b3=window(c,Lam,lo,hi)
    print(Lam,round(c,5),a2,b3,flush=True)
P=np.array(pts[-3:]); cc=P[:,0]; w2=P[:,1]**2
k,b=np.polyfit(cc,w2,1)
print("RESULT Lam",Lam,"cstar_fit",-b/k,"slope",k,"lastpts",pts[-3:],flush=True)
