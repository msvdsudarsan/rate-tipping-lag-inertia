import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, csv, sys
sys.path.insert(0,'.')
from common import Lam_c_bessel
fig=plt.figure(figsize=(13.28,5.31),dpi=120)
ax=fig.add_axes([0.05,0.16,0.26,0.66]); t=np.linspace(0,6,300)
for tau,col in [(0.05,'#d55e00'),(0.6,'#0072b2')]:
    ax.plot(t,3*(1-np.exp(-t/tau)),color=col,lw=2,label=r'lag $\tau=%g$'%tau)
ax.set_title('Lagged input $\\lambda(t)$',fontsize=13); ax.set_xlabel('time'); ax.legend(frameon=False)
ax2=fig.add_axes([0.40,0.16,0.26,0.66])
ts=np.linspace(0.05,3,60); L=[Lam_c_bessel(x) for x in ts]
ax2.plot(ts,L,color='k',lw=2); ax2.fill_between(ts,L,max(L)*1.05,color='#9ecae1',alpha=.6)
ax2.text(1.0,max(L)*0.95,'tips',ha='center'); ax2.text(1.8,2.3,'tracks',ha='center')
ax2.set_title(r'Exact threshold $\Lambda_c=j^2_{2\tau-1,1}/(4\tau)$',fontsize=13); ax2.set_xlabel(r'lag $\tau$'); ax2.set_ylabel(r'critical shift $\Lambda_c$'); ax2.set_ylim(1.8,max(L)*1.05)
ax3=fig.add_axes([0.74,0.16,0.24,0.66])
rd=[r for r in list(csv.reader(open('../data/p3_edges.csv')))[1:] if len(r)==4]
c=[float(r[0]) for r in rd]; a2=[float(r[2]) for r in rd]; b3=[float(r[3]) for r in rd]
ax3.fill_between(c,a2,b3,color='#9ecae1'); ax3.plot(c,a2,'k',lw=1); ax3.plot(c,b3,'k',lw=1)
ax3.set_title('Weak damping: closing\nre-entrant tipping window',fontsize=13); ax3.set_xlabel('damping $c$'); ax3.set_ylabel(r'lag $\tau$')
for a in (ax,ax2,ax3): a.spines[['top','right']].set_visible(False)
fig.savefig('../figs/graphical_abstract.png',dpi=120)
