"""Generate LaTeX tables and numbers.tex from the CSV/raw outputs in ../data."""
import csv, re, numpy as np, os, glob
D = '../data/'; O = '../'
def tex_num(x, nd=10): return ('%.' + str(nd) + 'f') % x
# ---------- ramp table
rows = list(csv.DictReader(open(D+'table_T1prime_ramp.csv')))
with open(O+'tab_ramp.tex', 'w') as f:
    f.write('\\begin{tabular}{rccc}\\toprule\n$Q$ & $\\Lc(Q)$ closed form & $\\Lc(Q)$ ODE bisection & $|\\text{diff}|$\\\\\\midrule\n')
    for r in rows:
        f.write('%g & %s & %s & $%.0e$\\\\\n' % (float(r['Q']), tex_num(float(r['Lambda_c_closed'])), tex_num(float(r['Lambda_c_ODE'])), max(float(r['abs_diff']), 1e-16)))
    f.write('\\bottomrule\\end{tabular}\n')
rows = list(csv.DictReader(open(D+'table_T1_filtered.csv')))
with open(O+'tab_filtered.tex', 'w') as f:
    f.write('\\begin{tabular}{rccc}\\toprule\n$\\tau$ & $j_{2\\tau-1,1}^2/(4\\tau)$ & ODE bisection & $|\\text{diff}|$\\\\\\midrule\n')
    for r in rows:
        f.write('%g & %s & %s & $%.0e$\\\\\n' % (float(r['tau']), tex_num(float(r['Lambda_c_Bessel'])), tex_num(float(r['Lambda_c_ODE'])), float(r['abs_diff'])))
    f.write('\\bottomrule\\end{tabular}\n')
# ---------- inertia table
rows = list(csv.DictReader(open(D+'table_T2_inertia.csv')))
with open(O+'tab_inertia.tex', 'w') as f:
    f.write('\\begin{tabular}{rrcccccc}\\toprule\n$\\Lambda$ & $c$ & $\\theta_0$ & $\\kappa$ & $\\theta^*_{\\rm full}$ & $\\theta^*_{\\rm red}$ & $\\theta_0+\\kappa/c^2$ & scaled remainder\\\\\\midrule\n')
    for r in rows:
        th0, k, c = float(r['theta0']), float(r['kappa']), float(r['c']); full = float(r['theta_full3D']); red = float(r['theta_reduced'])
        rem = '$%.4f$' % ((full-th0)*c*c-k) if c < 50 else '--'
        f.write('%.1f & %d & %.9f & %.7f & %.9f & %.9f & %.9f & %s\\\\\n' % (float(r['Lambda']), c, th0, k, full, red, th0+k/c**2, rem))
    f.write('\\bottomrule\\end{tabular}\n')
# ---------- c* table
def method1():
    res = {}
    for L, fn in [(2.5, 'g16_2.5.out'), (4.0, 'g16_4.out'), (5.0, 'g16_5.out')]:
        for line in open(D+'raw/'+fn):
            if line.startswith('RESULT'): res[L] = float(line.split()[4])
    # Lambda=3 from the continuation table g15_3.out (c, a2, b3, width): fit W^2 linear over last three points
    P = []
    for line in open(D+'raw/g15_3.out'):
        p = line.split()
        if len(p) >= 5 and p[2] != 'None': P.append((float(p[1]), float(p[4])))
    P = np.array(P[-3:]); k, b = np.polyfit(P[:, 0], P[:, 1]**2, 1); res[3.0] = -b/k
    return res
def method2():
    res = {}; tauf = {}
    for fn in sorted(glob.glob(D+'raw/z4_*.out')):
        L = float(re.findall(r'z4_([0-9.]+)\.out', fn)[0]); lines = open(fn).read().strip().split('\n')
        for line in lines:
            if 'c* (extrap' in line: res[L] = float(line.split('=')[1])
        pts = [l for l in lines if 'minT' in l]; tauf[L] = float(pts[-1].split()[-1])
    return res, tauf
m1 = method1(); m2, tf = method2()
with open(O+'tab_cstar.tex', 'w') as f:
    f.write('\\begin{tabular}{cccccc}\\toprule\n$\\Lambda$ & $c^*$ Method~1 & $c^*$ Method~2 & rel.\\ difference & $\\tau_f$\\\\\\midrule\n')
    for L in sorted(m1):
        if L in m2:
            f.write('%.1f & %.5f & %.5f & %.2f\\%% & %.3f\\\\\n' % (L, m1[L], m2[L], 100*abs(round(m1[L],5)-round(m2[L],5))/round(m2[L],5), tf[L]))
    f.write('\\bottomrule\\end{tabular}\n')
# ---------- swing table
if os.path.exists(D+'table_swing.csv'):
    rows = list(csv.DictReader(open(D+'table_swing.csv')))
    taus = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]; cs = [0.5, 0.15, 0.05]
    with open(O+'tab_swing.tex', 'w') as f:
        f.write('\\begin{tabular}{r' + 'c'*len(taus) + '}\\toprule\n$c\\backslash\\tau$ & ' + ' & '.join('%g' % t for t in taus) + '\\\\\\midrule\n')
        for c in cs:
            cells = []
            for t in taus:
                r = [x for x in rows if float(x['c']) == c and float(x['tau']) == t][0]
                tr = r['transitions']
                cells.append('none' if tr == '' else tr.split(';')[0].split('@')[1] + ('$^\\dagger$' if (';' in tr or tr.split('@')[0] != 'N>T') else ''))
            f.write('%g & %s\\\\\n' % (c, ' & '.join(cells)))
        f.write('\\bottomrule\\end{tabular}\n')
# ---------- numbers.tex
nums = {}
nums['cstar'] = '%.4f' % m2.get(3.0, m1[3.0]) if 3.0 in m2 else '%.4f' % m1[3.0]
open(O+'numbers.tex', 'w').write('\\newcommand{\\cstar}{%s}\n\\newcommand{\\edgeA}{%s}\n\\newcommand{\\edgeB}{%s}\n' % (nums['cstar'], os.environ.get('EDGEA', '0.648236'), os.environ.get('EDGEB', '0.752738')))
print('tables written'); print(m1, m2, tf)

# ---------- window width table
import csv as _csv
_rows=[r for r in list(_csv.reader(open(D+'p3_edges.csv')))[1:] if len(r)==4]
with open(O+'tab_W.tex','w') as f:
    f.write('\\begin{tabular}{ccc}\\toprule\n$c$ & $W$ & $W^2/(c^*-c)$\\\\\\midrule\n')
    for r in _rows:
        c=float(r[0])
        if c in (0.046,0.05,0.054,0.058,0.06,0.062):
            W=float(r[3])-float(r[2]); f.write('%.3f & %.4f & %.2f\\\\\n'%(c,W,W*W/(0.0628-c)))
    f.write('\\bottomrule\\end{tabular}\n')
