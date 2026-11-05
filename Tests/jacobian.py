"""Jacobien orbite -> parametres du MCMC, calcule numeriquement.

Independant du Fortran : on recode seulement le passage des elements aux
parametres, puis on derive par differences finies. Sert a verifier RATIO_PRIORS :

    ln(RAPQ) = Delta ln(prior) - Delta ln(J)
"""

import numpy as np


def stumpff(x):
    """Fonctions c1 et c3 pour x > 0 (orbite liee)."""
    r = np.sqrt(x)
    return np.sin(r) / r, (r - np.sin(r)) / (x * r)


def ecc_func(e):
    return 10.0 * e * e / (1.0 + 10.0 * e * e)


def classical_forward(x):
    """(a, P, e, i, angle, w, M0) -> parametres classiques du MCMC."""
    a, per, e, inc, angle, w, m0 = x
    exq = np.sqrt(1.0 - e * e)
    u = m0
    for _ in range(100):
        u = u - (u - e * np.sin(u) - m0) / (1.0 - e * np.cos(u))
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(u / 2.0), np.sqrt(1.0 - e) * np.cos(u / 2.0))
    lam = w + nu
    si = np.sqrt(exq) * np.sin(inc / 2.0)
    return np.array([np.log(a), e / exq * np.cos(w), e / exq * np.sin(w), si * np.cos(angle), si * np.sin(angle),
                     np.cos(lam) / per, np.sin(lam) / per])


def classical_inverse(p):
    k, h = p[1], p[2]
    s = np.hypot(k, h)
    e = s / np.sqrt(1.0 + s * s)
    exq = e / s
    w = np.arctan2(h, k)
    si = np.sqrt((p[3] ** 2 + p[4] ** 2) / exq)
    inc = 2.0 * np.arcsin(si)
    angle = np.arctan2(p[4], p[3])
    per = 1.0 / np.hypot(p[5], p[6])
    nu = np.arctan2(p[6], p[5]) - w
    u = 2.0 * np.arctan2(np.sqrt(1.0 - e) * np.sin(nu / 2.0), np.sqrt(1.0 + e) * np.cos(nu / 2.0))
    m0 = u - e * np.sin(u)
    return np.array([np.exp(p[0]), per, e, inc, angle, w, m0])


def universal_forward(x, t0_minus_tp=None):
    """(q, Pq, e, i, angle, w, dt) -> parametres universels, avec dt = T0 - tp."""
    q, per, e, inc, angle, w, dt = x
    mu = (2.0 * np.pi / per) ** 2 * q**3
    alpha = mu / q * (1.0 - e)
    s = dt / q
    for _ in range(200):
        xx = alpha * s * s
        c1, c3 = stumpff(xx)
        c0 = np.cos(np.sqrt(xx))
        c2 = (1.0 - c0) / xx
        t = mu * c3 * s**3 + q * c1 * s
        s = s - (t - dt) / (mu * c2 * s * s + q * c0)
    si = np.sin(inc / 2.0)
    return np.array([np.log(q), e * np.cos(w), e * np.sin(w), si * np.cos(angle), si * np.sin(angle), np.log(per),
                     s * ecc_func(e)])


def universal_inverse(p):
    q = np.exp(p[0])
    per = np.exp(p[5])
    e = np.hypot(p[1], p[2])
    w = np.arctan2(p[2], p[1])
    si = np.hypot(p[3], p[4])
    inc = 2.0 * np.arcsin(si)
    angle = np.arctan2(p[4], p[3])
    s = p[6] / ecc_func(e)
    mu = (2.0 * np.pi / per) ** 2 * q**3
    xx = mu / q * (1.0 - e) * s * s
    c1, c3 = stumpff(xx)
    dt = mu * c3 * s**3 + q * c1 * s
    return np.array([q, per, e, inc, angle, w, dt])


def ln_jacobian_block(p, version):
    """ln |det(d parametres / d elements)| pour une orbite."""
    inverse, forward = (classical_inverse, classical_forward) if version == "ellip" else (universal_inverse, universal_forward)
    x = inverse(p)
    matrix = np.zeros((7, 7))
    for j in range(7):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        matrix[:, j] = (forward(xp) - forward(xm)) / (2.0 * h)
    return np.linalg.slogdet(matrix)[1]


def ln_jacobian(p, npla, version, multipla, jitter):
    """ln J de tous les parametres, avec ln(masse) et ln(jitter) comme variables."""
    total = sum(ln_jacobian_block(p[7 * k:7 * k + 7], version) for k in range(npla))
    if multipla:
        total += -p[-1]  # d(ln m0)/d(m0) = 1/m0
    if jitter:
        total += -p[7 * npla + 1]  # d(ln j)/d(j) = 1/j
    return total


def scales(p, npla, version):
    """Pour chaque orbite : ln(a ou q), ln(P ou Pq), ln(masse dynamique)."""
    out = []
    for k in range(npla):
        b = p[7 * k : 7 * k + 7]
        per = 1.0 / np.hypot(b[5], b[6]) if version == "ellip" else np.exp(b[5])
        out.append((b[0], np.log(per), np.log((2.0 * np.pi / per) ** 2) + 3.0 * b[0]))
    return out


def known_prior_gap(point, base, npla, version, multipla, user_prior_star):
    """Ecart deja identifie entre RATIO_PRIORS et COMPUTE_MAP + jacobien.

    Version universelle (B10) : facteur (q/q')(P/P') par orbite.
    Masses ajustees (B11) : un facteur M'/M ou m0'/m0 par prior concerne, moins
    le m0'/m0 deja compte dans le jacobien de ln(m0).
    """
    a, b = scales(point, npla, version), scales(base, npla, version)
    gap = 0.0
    if version == "univ":
        gap -= sum((a[k][0] - b[k][0]) + (a[k][1] - b[k][1]) for k in range(npla))
    if multipla:
        dm0 = point[-1] - base[-1]
        dmass = [a[k][2] - b[k][2] for k in range(npla)]
        # Priors de masse positive : m0, puis M1 - m0, M2 - M1, etc.
        count_m0 = 2
        gap += (count_m0 - 1) * dm0 + dmass[0]
        for k in range(1, npla):
            gap += dmass[k - 1] + dmass[k]
        if user_prior_star:
            gap += dm0
    return gap
