import sys, os, math, json
import numpy as np
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, R)
import run_tests as rt, model, jacobian as jb

def forward(x):
    q, per, e, inc, angle, w, dt = x
    mu = (2*np.pi/per)**2*q**3; alpha = mu/q*(1-e)
    s = dt/q
    for _ in range(200):
        xx = alpha*s*s; c1, c3 = jb.stumpff(xx); c0 = np.cos(np.sqrt(xx)); c2 = (1-c0)/xx
        s = s - (mu*c3*s**3 + q*c1*s - dt)/(mu*c2*s*s + q*c0)
    xx = alpha*s*s; c1, c3 = jb.stumpff(xx); c2 = (1-np.cos(np.sqrt(xx)))/xx
    X = q - mu*c2*s*s; Y = np.sqrt(mu*q*(1+e))*c1*s
    lam = w + np.arctan2(Y, X); si = np.sin(inc/2)
    return np.array([np.log(q), e*np.cos(w), e*np.sin(w), si*np.cos(angle), si*np.sin(angle), np.cos(lam)/per, np.sin(lam)/per])

def inverse(p):
    q = np.exp(p[0]); per = 1/np.hypot(p[5], p[6]); e = np.hypot(p[1], p[2]); w = np.arctan2(p[2], p[1])
    inc = 2*np.arcsin(np.hypot(p[3], p[4])); angle = np.arctan2(p[4], p[3])
    nu = np.arctan2(p[6], p[5]) - w
    mu = (2*np.pi/per)**2*q**3; alpha = mu/q*(1-e)
    E = 2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2), np.sqrt(1+e)*np.cos(nu/2)); s = E/np.sqrt(alpha)
    xx = alpha*s*s; c1, c3 = jb.stumpff(xx)
    return np.array([q, per, e, inc, angle, w, mu*c3*s**3 + q*c1*s])

def lnj(p, npla, mul, jit):
    tot = 0
    for k in range(npla):
        x = inverse(p[7*k:7*k+7]); m = np.zeros((7, 7))
        for j in range(7):
            h = 1e-6*max(abs(x[j]), 1e-3); xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
            m[:, j] = (forward(xp) - forward(xm))/(2*h)
        tot += np.linalg.slogdet(m)[1]
    if mul: tot -= p[-1]
    if jit: tot -= p[7*npla+1]
    return tot

exe = sys.argv[1]
for case in sorted(os.listdir(R + '/cases')):
    res = rt.run_driver(exe, 'univ', case)
    if 'CHI2' not in res:
        print(case, 'ECHEC execution', res['stderr'][:100]); continue
    spec = json.load(open(f'{R}/cases/{case}/point.json'))
    py = model.chi2(spec['point'], spec['data']) if case != 'rel_seppa' else float('nan')
    path = res['MERIT'] - res['FMAP']
    npar = res['npar']; base = np.array([res['PM'][i] for i in range(1, npar+1)])
    j0 = lnj(base, res['npla'], res['multipla'], 'jit' in case); gaps = []
    for k in range(1, 5):
        if res['RAPQ'][k] <= 0: continue
        pt = np.array([res['PT'][1000*k+i] for i in range(1, npar+1)])
        exp = res['LNPRIOR'][k] - res['LNPRIOR0'] - (lnj(pt, res['npla'], res['multipla'], 'jit' in case) - j0)
        gap = jb.known_prior_gap(pt, base, res['npla'], 'ellip', res['multipla'], case == 'rel_rv_prior')  # B11 seul
        gaps.append(math.log(res['RAPQ'][k]) - exp - gap)
    print(f"{case:16s} chi2 MCMC {res['CHI2']:.6f} | LM {path:.6f} | Python {py:.6f} | prior : écart résiduel max {max(abs(g) for g in gaps):.1e}")
