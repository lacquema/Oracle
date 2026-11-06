#!/usr/bin/env python3
"""Refait l'ajustement de Gaia BH3 du notebook de l'ESA, independamment, et compare.

    python3 Tests/gaia/bh3_esa.py [CHEMIN/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.csv.gz]

Reference : github.com/esa/gaia-bhthree, branche gaia-dr4-prerelease (licence ESA-PL permissive),
notebook Gaia-DR4-prerelease_BH3_fit_astrometric_orbit.ipynb. Leur choix :
CCD gardes si used_by_agis_al, un point par CCD, modele de Holl+2023 sans terme de couleur,
ni perspective, ni jitter. Le fichier de donnees n'est pas dans le depot (ESA, pre-publication DR4).

Sert a verifier la lecture des donnees (colonnes, unites, temps) avant le code Gaia d'ORACLE.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

BH3 = 4318465066420528000
DEFAULT = os.path.expanduser("~/Oracle.env/GaiaDR4/GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.csv.gz")

# Valeurs affichees par le notebook de l'ESA (valeur, chiffres significatifs affiches)
ESA_SINGLE = {"ra": 0.71001, "dec": 0.26549, "plx": 0.98329, "mura": -30.94602, "mudec": -148.54974}
ESA_ORBIT = {"P": 4183.7, "e": 0.72782, "M0": 339.187, "ra": 3.4304, "dec": 2.5298, "plx": 1.65980,
             "mura": -28.3435, "mudec": -155.1704, "TIA": 2.4222, "TIB": 10.560, "TIF": 20.547, "TIG": -16.909}


def read(path):
    """Une ligne par CCD utilise par la solution de Gaia (used_by_agis_al)."""
    raw = pd.read_csv(path)
    raw = raw[raw.source_id == BH3]

    def split(value):
        return [x.strip() for x in str(value).strip("()").split(",")]

    rows = []
    for _, tr in raw.iterrows():
        cols = [split(tr[c]) for c in
                ("obs_time_tcb", "scan_pos_angle", "centroid_pos_al", "centroid_pos_error_al", "used_by_agis_al")]
        for time, angle, w, sigma, used in zip(*cols):
            if w in ("", "nan", "None") or used != "true":
                continue
            # temps TCB en ns depuis J2010, corrige du temps de trajet barycentrique
            rows.append((float(time) + tr.obs_time_bary_corr, float(angle), tr.parallax_factor_al,
                         float(w), float(sigma)))
    d = np.array(rows)
    t = d[:, 0] / 1e9 / 86400.0 / 365.25 - 7.5   # annees depuis J2017.5
    psi = np.radians(d[:, 1])
    return t, np.sin(psi), np.cos(psi), d[:, 2], d[:, 3], d[:, 4]


def single_star(t, s, c, fplx):
    return np.c_[s, c, fplx, t * s, t * c]


def with_orbit(t, s, c, fplx, period, ecc, m0):
    """Holl+2023 : w = 5 parametres + (A X + F Y) cos psi + (B X + G Y) sin psi."""
    mean = np.radians(m0) + 2.0 * np.pi * t * 365.25 / period
    ea = mean + ecc * np.sin(mean)
    for _ in range(80):
        ea -= (ea - ecc * np.sin(ea) - mean) / (1.0 - ecc * np.cos(ea))
    x, y = np.cos(ea) - ecc, np.sqrt(1.0 - ecc * ecc) * np.sin(ea)
    return np.c_[single_star(t, s, c, fplx), x * c, x * s, y * c, y * s]


def linear(design, w, sigma):
    a = design / sigma[:, None]
    p, *_ = np.linalg.lstsq(a, w / sigma, rcond=None)
    return p, a @ p - w / sigma


def compare(label, names, values, reference):
    print(f"== {label}")
    good = True
    for name, value in zip(names, values):
        ref = reference[name]
        digits = len(f"{ref}".split(".")[1]) if "." in f"{ref}" else 0
        # 1,5 unite du dernier chiffre affiche : l'ESA lit certaines colonnes en simple precision
        ok = abs(value - ref) <= 1.5 * 10.0 ** (-digits)
        good &= ok
        print(f"  {name:6s} {value:12.5f}   ESA {ref:<11}  {'OK' if ok else 'ÉCART'}")
    return good


def main():
    t, s, c, fplx, w, sigma = read(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
    print(f"{len(t)} mesures CCD (ESA : 558)")
    p, _ = linear(single_star(t, s, c, fplx), w, sigma)
    good = compare("Étoile seule", list(ESA_SINGLE), p, ESA_SINGLE)

    def residuals(q):
        return linear(with_orbit(t, s, c, fplx, *q), w, sigma)[1]

    fits = [least_squares(residuals, [p0, e0, m0], bounds=([500, 0, -720], [20000, 0.97, 720]))
            for p0 in (1667, 3000, 4200, 6000) for e0 in (0.3, 0.7) for m0 in (0, 90, 180, 270)]
    best = min(fits, key=lambda f: f.cost)
    period, ecc, m0 = best.x
    lin, res = linear(with_orbit(t, s, c, fplx, *best.x), w, sigma)
    good &= compare("Avec orbite", list(ESA_ORBIT), [period, ecc, m0 % 360.0, *lin], ESA_ORBIT)
    print(f"  chi2 = {np.sum(res ** 2):.1f} pour {len(t)} mesures")
    print("Bilan :", "accord avec l'ESA" if good else "ÉCART avec l'ESA")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
