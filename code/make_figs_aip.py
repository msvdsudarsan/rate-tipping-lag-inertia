import numpy as np, csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 150,
                     'pdf.fonttype': 42})
OI = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']
D = '../data/'; F = '../figs/'
# ---- Fig 1: exact thresholds
fig, ax = plt.subplots(1, 2, figsize=(6.69, 2.8))
taus = np.linspace(0.02, 4.5, 200); lc = [Lam_c_bessel(t) for t in taus]
ax[0].plot(taus, lc, color=OI[0], lw=1.6, label=r'$j_{2\tau-1,1}^2/(4\tau)$')
ax[0].axhline(2, color='0.5', ls=':', lw=1); ax[0].text(3.0, 2.12, r'$\Lambda=2$ (instantaneous step)', fontsize=8, color='0.35')
ax[0].plot(taus, taus, color='0.6', ls='--', lw=1); ax[0].text(3.3, 2.9, r'$\Lambda=\tau$', fontsize=8, color='0.4')
rows = list(csv.DictReader(open(D+'table_T1_filtered.csv')))
ax[0].plot([float(r['tau']) for r in rows], [float(r['Lambda_c_ODE']) for r in rows], 'o', mfc='none', color=OI[1], ms=5, label='ODE bisection')
ax[0].set_xlabel(r'lag $\tau$'); ax[0].set_ylabel(r'critical total shift $\Lambda_c$'); 
ax[0].fill_between(taus, lc, 8, color=OI[0], alpha=0.08); ax[0].text(0.15, 6.6, 'tipping', fontsize=8, color=OI[0]); ax[0].text(3.2, 1.35, 'tracking', fontsize=8, color='0.3')
ax[0].set_ylim(1.0, 8)
Qs = np.linspace(1.02, 20, 300)
ax[1].plot(Qs, [Lam_c_ramp(q) for q in Qs], color=OI[0], lw=1.6, label=r'$\frac{2Q}{\omega}\arctan\frac{1}{\omega}$')
rows = list(csv.DictReader(open(D+'table_T1prime_ramp.csv')))
ax[1].plot([float(r['Q']) for r in rows], [float(r['Lambda_c_ODE']) for r in rows], 'o', mfc='none', color=OI[1], ms=5, label='ODE bisection')
ax[1].axhline(2, color='0.5', ls=':', lw=1); ax[1].set_ylim(1.5, 8); ax[1].set_xlim(1, 20)
ax[1].set_xlabel(r'ramp rate $Q$'); ax[1].set_ylabel(r'critical total shift $\Lambda_c$'); ax[1].legend(frameon=False, fontsize=8)
ax[0].set_title('(a) lagged step', fontsize=9, loc='left'); ax[1].set_title('(b) ramp that stops', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(F+'aip_fig1_exact.pdf'); plt.close(fig)
# ---- Fig 2: inertial correction
rows = list(csv.DictReader(open(D+'table_T2_inertia.csv')))
fig, ax = plt.subplots(1, 2, figsize=(6.69, 2.3))
for L, col in zip([2.5, 3.0, 4.0], OI):
    rr = [r for r in rows if float(r['Lambda']) == L]; c = np.array([float(r['c']) for r in rr])
    th0 = float(rr[0]['theta0']); k = float(rr[0]['kappa']); full = np.array([float(r['theta_full3D']) for r in rr])
    ax[0].plot(c, (full-th0)*c**2, 'o-', color=col, ms=4, label=r'$\Lambda=%.1f$' % L); ax[0].axhline(k, color=col, ls='--', lw=0.9)
    ax[1].loglog(c, np.abs(full-th0-k/c**2), 'o-', color=col, ms=4)
ax[0].set_xscale('log'); ax[0].set_xlabel(r'damping $c$'); ax[0].set_ylabel(r'$c^2(\theta^*_{\rm full}-\theta_0)$'); ax[0].legend(frameon=False, fontsize=8, loc='center right', bbox_to_anchor=(1.0, 0.45))
cc = np.array([10, 100.]); ax[1].loglog(cc, 1.9/cc**4, color='0.5', ls=':'); ax[1].text(20, 2e-5, r'$\propto c^{-4}$', fontsize=8, color='0.4')
ax[1].set_xlabel(r'damping $c$'); ax[1].set_ylabel(r'$|\theta^*_{\rm full}-\theta_0-\kappa/c^2|$')
ax[0].set_title(r'(a) dashed lines: $\kappa(\Lambda)$', fontsize=9, loc='left'); ax[1].set_title('(b) remainder', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(F+'aip_fig2_inertia.pdf'); plt.close(fig)
print('fig1, fig2 done')
try:
    # ---- Fig 3: window map
    import csv as _csv
    rd = list(_csv.reader(open(D+'p3_map.csv'))); hdr = rd[0]; taus = np.array([float(h.split('_')[1]) for h in hdr[1:]])
    cs = np.array([float(r[0]) for r in rd[1:]]); M = np.array([[int(x) for x in r[1:]] for r in rd[1:]])
    ed = {}
    for r in list(_csv.reader(open(D+'p3_edges.csv')))[1:]:
        ed[float(r[0])] = [float(x) for x in r[1:]]
    fig, ax = plt.subplots(figsize=(3.37, 2.9))
    dc = cs[1]-cs[0]
    ax.imshow(M.T, origin='lower', aspect='auto', extent=[cs[0]-dc/2, cs[-1]+dc/2, taus[0], taus[-1]], cmap=matplotlib.colors.ListedColormap(['#f3f3f3', '#9ecae1']), interpolation='nearest')
    for k, (lab, col) in enumerate([('$a_1$', OI[1]), ('$a_2$', OI[2]), ('$b_3$', OI[3])]):
        xs = [c for c in cs if len(ed[c]) == 3]; ys = [ed[c][k] for c in xs]
        ax.plot(xs, ys, '-', color=col, lw=1.4, label=lab)
    ax.set_xlabel(r'damping $c$'); ax.set_ylabel(r'lag $\tau$'); ax.legend(frameon=True, fontsize=8, loc='upper right')
    ax.text(0.012, 0.42, 'tip', fontsize=8); ax.text(0.014, 1.5, 'no tip', fontsize=8); ax.text(0.03, 1.1, 'tip window', fontsize=8); ax.text(0.0105, 0.69, 'no-tip\nisland', fontsize=8)
    ax.axvline(float(open('../numbers.tex').read().split('cstar}{')[1].split('}')[0]), color='0.4', ls=':', lw=0.9)
    ax.set_xlim(0, cs[-1]); ax.set_ylim(0.3, 1.6)
    fig.tight_layout(); fig.savefig(F+'aip_fig3_window.pdf'); plt.close(fig)
except FileNotFoundError as e:
    print('skipped', e)
try:
    # ---- Fig 4: swing
    rows = list(_csv.DictReader(open(D+'table_swing.csv')))
    fig, ax = plt.subplots(figsize=(3.37, 2.6))
    T=[0.0,0.25,0.5,1.0,2.0,4.0]; pos=list(range(len(T)))
    for c, col in zip([0.15, 0.05], OI):
        xs=[]; ys=[]
        for p,t in zip(pos,T):
            r=[x for x in rows if float(x['c'])==c and float(x['tau'])==t][0]['transitions']
            if r=='' or '@' not in r:
                ax.plot([p],[9],'x',color=col,ms=6); continue
            xs.append(p); ys.append(float(r.split(';')[0].split('@')[1]))
        ax.plot(xs, ys, 'o-', color=col, ms=4, label='$c=%.2f$' % c)
    ax.set_yscale('log'); ax.set_xticks(pos); ax.set_xticklabels(['0','0.25','0.5','1','2','4'])
    ax.set_xlabel(r'lag $\tau$'); ax.set_ylabel(r'critical ramp rate $r_c$'); ax.legend(frameon=False, fontsize=8, loc='center left')
    ax.set_ylim(0.5,12)
    ax.text(4.5,9,'$\\times$: no tipping at any tested $r$',fontsize=8,ha='right',va='bottom',color='0.3')
    fig.tight_layout(); fig.savefig(F+'aip_fig4_swing.pdf'); plt.close(fig)
    print('fig3, fig4 done')
except FileNotFoundError as e:
    print('skipped', e)
