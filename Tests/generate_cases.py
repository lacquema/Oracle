#!/usr/bin/env python3
"""Genere Tests/cases/ : donnees simulees et fichiers d'entree du Fortran.

A ne relancer que pour ajouter un cas : les fichiers produits sont suivis par
git et la reference des tests en depend.
"""

import json
import math
import os

import numpy as np

import cases_definition as cd
from model import System

HERE = os.path.dirname(os.path.abspath(__file__))


def make_data(name, spec, rng):
    """Tire les donnees d'un cas et renvoie les lignes du fichier et les valeurs."""
    system = System(cd.MSTAR_TRUE, spec["orbits"])
    plx = 1000.0 / cd.DIST_PC
    text = ["2400000.   JD offset", "-" * 40]
    values = {}
    separator = False

    if spec.get("rel"):
        rho = spec.get("corr", 0.0)
        seppa = spec.get("position_format", 1) == 2
        rows, lines = [], []
        for ipla in range(len(spec["orbits"])):
            for t in np.linspace(54000.0, 57500.0, 12):
                x, y = system.relative_position(ipla, t)
                xo = float(x * plx + rng.normal(0, cd.SIG_POS_MAS))
                yo = float(y * plx + rng.normal(0, cd.SIG_POS_MAS))
                rows.append([ipla + 1, float(t), xo, yo, cd.SIG_POS_MAS, cd.SIG_POS_MAS, rho])
                if seppa:
                    sep = math.hypot(xo, yo)
                    pa = math.degrees(math.atan2(yo, xo))
                    dpa = math.degrees(cd.SIG_POS_MAS / sep)
                    lines.append(f"{ipla + 1} {t:.3f} {sep:.6f} {pa:.6f} {cd.SIG_POS_MAS:.6f} {dpa:.6f} {rho:.3f}")
                else:
                    corr = f" {rho:.3f}" if "corr" in spec else ""
                    lines.append(f"{ipla + 1} {t:.3f} {xo:.6f} {yo:.6f} {cd.SIG_POS_MAS:.6f} {cd.SIG_POS_MAS:.6f}{corr}")
        text += ["1", "I JD X Y dX dY Corr", "-" * 40] + lines
        separator = True
        # En (Sep, PA), la conversion passe par io.f : pas de recalcul independant.
        if not seppa:
            values["rel"] = [[r[0], r[1], round(r[2], 6), round(r[3], 6), r[4], r[5], r[6]] for r in rows]

    if spec.get("rv"):
        if separator:
            text.append("-" * 40)
        rows = []
        for t in np.linspace(54000.0, 57500.0, 40):
            v = float(system.star_rv(cd.V0_TRUE, t) + rng.normal(0, cd.SIG_RV_KMS))
            rows.append([float(t), round(v, 6), cd.SIG_RV_KMS])
        text += ["2", "JD RV dRV", "-" * 40] + [f"{r[0]:.3f} {r[1]:.6f} {r[2]:.6f}" for r in rows]
        values["rv"] = [[round(r[0], 3), r[1], r[2]] for r in rows]
        separator = True

    if spec.get("rvp"):
        if separator:
            text.append("-" * 40)
        rows = []
        for t in np.linspace(54500.0, 57000.0, 12):
            v = float(system.planet_rv(0, t) + rng.normal(0, cd.SIG_RVP_KMS))
            rows.append([1, float(t), round(v, 6), cd.SIG_RVP_KMS])
        text += ["3", "I JD RV dRV", "-" * 40] + [f"{r[0]} {r[1]:.3f} {r[2]:.6f} {r[3]:.6f}" for r in rows]
        values["rvp"] = [[r[0], round(r[1], 3), r[2], r[3]] for r in rows]

    if "rel" in values:
        values["rel"] = [[r[0], round(r[1], 3), r[2], r[3], r[4], r[5], r[6]] for r in values["rel"]]
    return "\n".join(text) + "\n", values


def make_input(spec, points, univ):
    """Ecrit le fichier d'entree, dans l'ordre lu par INIT_SIMU et INIT_DATA."""
    isdata = [int(bool(spec.get(k))) for k in ("rel", "rv", "rvp")] + [0]
    fmt = spec.get("position_format", 1)
    lines = [
        "2",
        " ".join(map(str, isdata)),
        str(len(points)),
        "data.dat",
        "1d-7",
        f"2 {fmt} {int('corr' in spec)}",
    ]
    if spec.get("rv"):
        lines.append(str(int("jitter" in spec)))
    lines += [
        f"{cd.DIST_PC} pc",
        str(cd.MSTAR_POINT),
        "res.dat",
        "dump.dat",
        "100000000",
    ]
    if "prior" in spec:
        lines += ["1", spec["prior"]]
    else:
        lines.append("0")
    lines.append(str(cd.T0))
    if spec.get("rv"):
        lines.append(f"{cd.V0_POINT} {spec.get('jitter', 0.0)}")
    for o in points:
        scale = o["a"] * (1.0 - o["e"]) if univ else o["a"]
        lines.append(
            f"{o['m_mjup']!r} mj {scale!r} {o['e']!r} {o['i']!r} {o['node']!r} {o['w']!r} {o['tp']!r}"
        )
    if not univ:
        lines.append("0.0 1000000.0")
    lines += ["0.0 1000.0", f"0.0 {spec.get('emax', 0.95)}"]
    return "\n".join(lines) + "\n"


def main():
    rng = np.random.default_rng(20261105)
    for name, spec in cd.CASES.items():
        folder = os.path.join(HERE, "cases", name)
        os.makedirs(folder, exist_ok=True)
        text, values = make_data(name, spec, rng)
        points = [cd.offset(o) for o in spec["orbits"]]
        with open(os.path.join(folder, "data.dat"), "w") as f:
            f.write(text)
        with open(os.path.join(folder, "input_ellip.txt"), "w") as f:
            f.write(make_input(spec, points, False))
        with open(os.path.join(folder, "input_univ.txt"), "w") as f:
            f.write(make_input(spec, points, True))
        point = dict(mstar=cd.MSTAR_POINT, dist_pc=cd.DIST_PC, orbits=points)
        if spec.get("rv"):
            point["v0"] = cd.V0_POINT
            point["jitter"] = spec.get("jitter", 0.0)
        with open(os.path.join(folder, "point.json"), "w") as f:
            json.dump(dict(point=point, data=values), f, indent=1)
    print(f"{len(cd.CASES)} cas generes")


if __name__ == "__main__":
    main()
