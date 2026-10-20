import mpmath as mp
mp.mp.dps=25
def jz(nu,k):
    # bracket via sign changes of J_nu on a fine grid
    f=lambda x: mp.besselj(nu,x)
    xs=[mp.mpf(i)/50 for i in range(1,50*(8+4*k))]
    c=0
    for a,b in zip(xs,xs[1:]):
        if f(a)*f(b)<0:
            c+=1
            if c==k: return mp.findroot(f,(a,b),solver='anderson')
for nu in [-0.9,-0.6,-0.3,0,0.5,2,5]:
  for k in [1,2,3]:
    j=jz(nu,k); h=mp.mpf('1e-6')
    d=(jz(nu+h,k)-jz(nu-h,k))/(2*h)
    I=2*j*mp.quad(lambda t: mp.besselk(0,2*j*mp.sinh(t))*mp.e**(-2*nu*t),[0,0.5,2,10,40])
    print(nu,k,mp.nstr(d,9),mp.nstr(I,9))
