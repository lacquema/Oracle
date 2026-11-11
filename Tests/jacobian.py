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


def universal_forward(x):
    """(q, Pq, e, i, angle, w, dt) -> parametres universels du MCMC (M011), avec dt = T0 - tp.

    Parametres 6 et 7 = (cos, sin)(lambda)/Pq, lambda = w + anomalie vraie a T0.
    """
    q, per, e, inc, angle, w, dt = x
    # Orbite liee : anomalie vraie a T0 par l'equation de Kepler elliptique (depart de Danby),
    # equivalente a la variable universelle s et sure pres de e = 1, ou un Newton en s parti
    # de dt/q pouvait ne pas converger (cas quasi paraboliques, 2026-11-11)
    mu = (2.0 * np.pi / per) ** 2 * q**3
    a = q / (1.0 - e)
    m0 = np.sqrt(mu / a**3) * dt
    m0 = np.arctan2(np.sin(m0), np.cos(m0))
    u = m0 + 0.85 * e * np.sign(np.sin(m0))
    for _ in range(100):
        u = u - (u - e * np.sin(u) - m0) / (1.0 - e * np.cos(u))
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(u / 2.0), np.sqrt(1.0 - e) * np.cos(u / 2.0))
    lam = w + nu
    si = np.sin(inc / 2.0)
    return np.array([np.log(q), e * np.cos(w), e * np.sin(w), si * np.cos(angle), si * np.sin(angle),
                     np.cos(lam) / per, np.sin(lam) / per])


def universal_inverse(p):
    """Inverse de universal_forward, pour une orbite liee."""
    q = np.exp(p[0])
    per = 1.0 / np.hypot(p[5], p[6])
    e = np.hypot(p[1], p[2])
    w = np.arctan2(p[2], p[1])
    si = np.hypot(p[3], p[4])
    inc = 2.0 * np.arcsin(si)
    angle = np.arctan2(p[4], p[3])
    nu = np.arctan2(p[6], p[5]) - w
    mu = (2.0 * np.pi / per) ** 2 * q**3
    alpha = mu / q * (1.0 - e)
    u = 2.0 * np.arctan2(np.sqrt(1.0 - e) * np.sin(nu / 2.0), np.sqrt(1.0 + e) * np.cos(nu / 2.0))
    s = u / np.sqrt(alpha)
    c1, c3 = stumpff(alpha * s * s)
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
        per = 1.0 / np.hypot(b[5], b[6])  # Les deux versions depuis M011
        out.append((b[0], np.log(per), np.log((2.0 * np.pi / per) ** 2) + 3.0 * b[0]))
    return out


def known_prior_gap(point, base, npla, version, multipla, user_prior_star):
    """Ecart deja identifie entre RATIO_PRIORS et COMPUTE_MAP + jacobien.

    B10 corrigee par M011, B11 par M013 : plus aucun ecart attendu.
    """
    return 0.0
