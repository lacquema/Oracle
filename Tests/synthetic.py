#!/usr/bin/env python3
"""Systemes synthetiques : le vrai programme doit retrouver une verite connue.

Chaque systeme est invente ici, ses donnees sont tirees avec Tests/model.py,
puis le vrai programme (compile avec OpenMP, comme en production) tourne jusqu'au
bout. Ce test est long, il n'est donc pas dans `make test`.

    python3 Tests/synthetic.py liste
    python3 Tests/synthetic.py lance NOM [--univ] [--rapide] [--threads 6] [--limite 3600]
    python3 Tests/synthetic.py analyse NOM [--univ]

Le cas `sans_donnees` n'a aucune mesure : le MCMC doit alors redonner le prior.
"""

import argparse
import gzip
import math
import os
import shutil
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(ROOT, "Tests")
sys.path.insert(0, TESTS)
import model  # noqa: E402

SUB = os.environ.get("ORACLE_SUB", os.path.join(ROOT, "Code", "mcmc", "sub"))
MAIN = os.path.join(ROOT, "Code", "mcmc", "main")
WORK = os.environ.get("ORACLE_SYNTH", os.path.join(TESTS, "build_synthetic"))

ORBIT_A = dict(m_mjup=30.0, a=5.0, e=0.35, i=40.0, node=120.0, w=50.0, tp=55500.0)
ORBIT_B1 = dict(m_mjup=300.0, a=3.0, e=0.25, i=50.0, node=70.0, w=30.0, tp=55000.0)
ORBIT_B2 = dict(m_mjup=500.0, a=12.0, e=0.15, i=55.0, node=80.0, w=210.0, tp=56500.0)

SYSTEMS = {
    # Une orbite bien couverte, avec vitesses de l'etoile : masses separees
    "rel_vr_etoile": dict(mstar=1.2, orbits=[ORBIT_A], rel=30, rv=30, span=6000.0),
    # Deux orbites, relative seule, prior sur la masse totale (comme GG Tau)
    "deux_orbites": dict(mstar=0.6, orbits=[ORBIT_B1, ORBIT_B2], rel=30, span=16000.0,
                         prior=("1 1 1", 1.363, 0.05)),
    # Une orbite, relative et vitesses de la planete, prior sur la masse totale
    "rel_vr_planete": dict(mstar=1.2, orbits=[ORBIT_A], rel=30, rvp=12, span=6000.0,
                           prior=("1 1", 1.2286, 0.05)),
    # Orbite peu couverte, proche de 51 Eri : arc de 10 ans sur environ 26 ans
    "arc_court": dict(mstar=1.5, orbits=[dict(m_mjup=3.0, a=10.0, e=0.55, i=135.0, node=60.0, w=200.0,
                                             tp=62000.0)], rel=25, rvp=4, span=3650.0),
    # Aucune mesure : le MCMC doit redonner le prior
    "sans_donnees": dict(mstar=1.0, orbits=[ORBIT_A], span=6000.0),
}

DIST_PC = 40.0
SIG_POS = 1.0
SIG_RV = 0.01
SIG_RVP = 0.5
V0 = 3.0
T_START = 54000.0
T0 = 56000.0
RANGE_P = (100.0, 100000.0)
RANGE_A = (0.5, 100.0)
RANGE_E = (0.0, 0.9)


def write_data(name, spec, folder):
    """Tire les donnees (graine fixe) et ecrit le fichier lu par le Fortran."""
    rng = np.random.default_rng(abs(hash_name(name)))
    system = model.System(spec["mstar"], spec["orbits"])
    plx = 1000.0 / DIST_PC
    lines = ["2400000. JD offset", "-" * 40]
    if spec.get("rel"):
        lines += ["1", "I JD X Y dX dY", "-" * 40]
        for ipla in range(len(spec["orbits"])):
            for t in np.linspace(T_START, T_START + spec["span"], spec["rel"]):
                x, y = system.relative_position(ipla, t)
                xo = float(x) * plx + rng.normal(0, SIG_POS)
                yo = float(y) * plx + rng.normal(0, SIG_POS)
                lines.append(f"{ipla + 1} {t:.3f} {xo:.6f} {yo:.6f} {SIG_POS:.6f} {SIG_POS:.6f}")
    elif not spec.get("rv"):
        lines += ["1", "I JD X Y dX dY", "-" * 40]
    if spec.get("rv"):
        lines += ["-" * 40, "2", "JD RV dRV", "-" * 40]
        for t in np.linspace(T_START, T_START + spec["span"], spec["rv"]):
            v = float(system.star_rv(V0, t)) + rng.normal(0, SIG_RV)
            lines.append(f"{t:.3f} {v:.6f} {SIG_RV:.6f}")
    if spec.get("rvp"):
        lines += ["-" * 40, "3", "I JD RV dRV", "-" * 40]
        for t in np.linspace(T_START + spec["span"] / 2, T_START + spec["span"], spec["rvp"]):
            v = float(system.planet_rv(0, t)) + rng.normal(0, SIG_RVP)
            lines.append(f"1 {t:.3f} {v:.6f} {SIG_RVP:.6f}")
    with open(os.path.join(folder, "data.dat"), "w") as f:
        f.write("\n".join(lines) + "\n")


def hash_name(name):
    """Graine stable d'une execution a l'autre."""
    return sum((i + 1) * ord(c) for i, c in enumerate(name)) + 20261105


def write_input(spec, folder, univ):
    """Ecrit l'entree : point de depart volontairement decale de la verite."""
    isdata = [int(bool(spec.get(k))) for k in ("rel", "rv", "rvp")] + [0]
    if not any(isdata):
        isdata[0] = 1
    choice = "2" if spec.get("rel") else "3"
    lines = [choice, " ".join(map(str, isdata)), str(len(spec["orbits"])), "data.dat", "1d-7", "2 1 0"]
    if spec.get("rv"):
        lines.append("0")
    lines += [f"{DIST_PC} pc", str(spec["mstar"] * 1.05), "results.dat", "dump.dat", "100000000"]
    if "prior" in spec:
        coefficients, mean, sdev = spec["prior"]
        lines += ["1", coefficients, "1", f"{mean} {sdev} ms"]
    else:
        lines.append("0")
    lines.append(str(T0))
    if spec.get("rv"):
        lines.append(f"{V0 + 0.05} 0.0")
    for o in spec["orbits"]:
        a = o["a"] * 1.03
        e = o["e"] + 0.03
        scale = a * (1.0 - e) if univ else a
        lines.append(
            f"{o['m_mjup'] * 1.2} mj {scale} {e} {o['i'] + 4} {o['node'] + 4} {o['w'] + 4} {o['tp'] + 40}"
        )
    if not univ:
        lines.append(f"{RANGE_P[0]} {RANGE_P[1]}")
    lines += [f"{RANGE_A[0]} {RANGE_A[1]}", f"{RANGE_E[0]} {RANGE_E[1]}", "exit", ""]
    with open(os.path.join(folder, "input.txt"), "w") as f:
        f.write("\n".join(lines))


def fast_sources(ninfer):
    """Copie le code et reduit le nombre d'orbites gardees. Le depot n'est pas touche."""
    root = os.path.join(WORK, f"src_{ninfer}")
    for name, source in (("main", MAIN), ("sub", SUB)):
        target = os.path.join(root, name)
        if os.path.isdir(target):
            shutil.rmtree(target)
        shutil.copytree(source, target, ignore=shutil.ignore_patterns("~*", "*.mod", "*.o"))
    path = os.path.join(root, "sub", "commons.f")
    with open(path, newline="") as f:
        text = f.read()
    old = "NINFER = 50000"
    if text.count(old) != 1:
        raise SystemExit("NINFER introuvable dans commons.f")
    with open(path, "w", newline="") as f:
        f.write(text.replace(old, f"NINFER = {ninfer}"))
    return os.path.join(root, "main"), os.path.join(root, "sub")


def build(univ, ninfer=None):
    """Compile le programme avec les memes options que le Makefile (-O3 -fopenmp)."""
    kind = ("univ" if univ else "ellip") + (f"_{ninfer}" if ninfer else "")
    folder = os.path.join(WORK, "bin", kind)
    os.makedirs(folder, exist_ok=True)
    main_dir, sub_dir = fast_sources(ninfer) if ninfer else (MAIN, SUB)
    base = ["gfortran", "-w", "-O3", "-fopenmp", "-J", folder, "-I", folder, "-I", sub_dir]
    main = "astrom_univ_mcmco.f" if univ else "astrom_mcmco.f"
    kep = "kepu.f" if univ else "kepellip.f"
    objects = [os.path.join(folder, "main.o")]
    # gfortran cherche INCLUDE et les modules (.mod) d'abord dans le dossier courant :
    # compiler depuis main_dir pour INCLUDE, puis depuis folder, qui ne contient que le bon data.mod
    subprocess.run([*base, "-c", os.path.join(main_dir, main), "-o", objects[0]], check=True, capture_output=True,
                   cwd=main_dir)
    for source in ["utils.f", "mcmc.f", "mrqfit.f", "io.f", kep]:
        obj = os.path.join(folder, source.replace(".f", ".o"))
        subprocess.run([*base, "-c", os.path.join(sub_dir, source), "-o", obj], check=True, capture_output=True,
                       cwd=folder)
        objects.append(obj)
    exe = os.path.join(folder, "oracle")
    subprocess.run(["gfortran", "-w", "-O3", "-fopenmp", *objects, "-o", exe], check=True, capture_output=True,
                   cwd=folder)
    return exe


def case_folder(name, univ):
    return os.path.join(WORK, name + ("_univ" if univ else "_ellip"))


def launch(name, univ, threads, limit, ninfer=None):
    spec = SYSTEMS[name]
    folder = case_folder(name, univ)
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    write_data(name, spec, folder)
    write_input(spec, folder, univ)
    exe = build(univ, ninfer)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), STACKSIZE="1000000", GFORTRAN_UNBUFFERED_ALL="1")
    with open(os.path.join(folder, "input.txt")) as stdin, open(os.path.join(folder, "out.txt"), "w") as out:
        try:
            code = subprocess.run([exe], stdin=stdin, stdout=out, stderr=subprocess.STDOUT, cwd=folder, env=env,
                                  timeout=limit).returncode
        except subprocess.TimeoutExpired:
            code = "limite de temps atteinte"
    print(f"{name} : fin, code {code}")


def read_results(folder):
    path = os.path.join(folder, "results.dat.gz")
    if not os.path.exists(path):
        return None, None
    with gzip.open(path, "rt") as f:
        lines = f.read().split("\n")
    last = max(i for i, line in enumerate(lines) if line.strip())
    count, nparam, nbody = [int(float(x)) for x in lines[last - 1].split()]
    data = np.fromstring(lines[last], sep=" ").reshape(nbody, nparam, count)
    names = {13: "a P e w i W tp m V0 Jit m0 Chi2 Map".split(), 11: "a P e w i W tp m m0 Chi2 Map".split()}[nparam]
    return names, data


def angle_offset(samples, truth):
    """Ecart angulaire signe a la verite, en degres, dans [-180, 180]."""
    return (np.degrees(samples) - truth + 180.0) % 360.0 - 180.0


def report_line(label, samples, truth):
    """Retourne l'ecart a la verite en unites d'ecart-type du posterior."""
    q16, q50, q84 = np.percentile(samples, [16, 50, 84])
    sigma = 0.5 * (q84 - q16)
    z = (q50 - truth) / sigma if sigma > 0 else float("inf")
    inside = np.percentile(samples, 0.15) <= truth <= np.percentile(samples, 99.85)
    print(f"  {label:<12} vrai {truth:<10.5g} trouvé {q50:<10.5g} +{q84 - q50:<9.3g} -{q50 - q16:<9.3g}"
          f" écart {z:+6.2f} sigma  {'OK' if inside else 'HORS 3 sigma'}")
    return inside


def analyse_prior(names, data):
    """Sans donnees : chaque variable doit suivre le prior du code."""
    get = lambda n: data[0, names.index(n)][1:]  # noqa: E731
    a, per, e, i, tp = get("a"), get("P"), get("e"), get("i"), get("tp")
    print(f"  {len(a)} orbites. Pour une loi uniforme sur [0,1] : moyenne 0,5 et écart-type 0,2887")
    checks = {
        "ln a": (np.log(a) - math.log(RANGE_A[0])) / (math.log(RANGE_A[1]) - math.log(RANGE_A[0])),
        "ln P": (np.log(per) - math.log(RANGE_P[0])) / (math.log(RANGE_P[1]) - math.log(RANGE_P[0])),
        "e": (e - RANGE_E[0]) / (RANGE_E[1] - RANGE_E[0]),
        "cos i": 0.5 * (1.0 + np.cos(i)),
    }
    good = True
    for label, values in checks.items():
        mean, std = values.mean(), values.std()
        ok = abs(mean - 0.5) < 0.02 and abs(std - 0.2887) < 0.02
        good = good and ok
        print(f"  {label:<6} moyenne {mean:.4f}  écart-type {std:.4f}  min {values.min():.3f} max {values.max():.3f}"
              f"  {'OK' if ok else 'ÉCART'}")
    n = 2.0 * np.pi / per
    phase = (n * (T0 - tp)) % (2.0 * np.pi) / (2.0 * np.pi)
    mean, std = phase.mean(), phase.std()
    ok = abs(mean - 0.5) < 0.02 and abs(std - 0.2887) < 0.02
    good = good and ok
    print(f"  phase  moyenne {mean:.4f}  écart-type {std:.4f}  {'OK' if ok else 'ÉCART'}")
    return good


def analyse(name, univ):
    spec = SYSTEMS[name]
    folder = case_folder(name, univ)
    names, data = read_results(folder)
    print(f"== {name} ({'universel' if univ else 'classique'})")
    with open(os.path.join(folder, "out.txt"), errors="replace") as f:
        text = f.read()
    print(f"  Terminé : {'oui' if 'Computation finished' in text else 'non'}")
    if names is None:
        print("  Pas encore de résultats.")
        return False
    if name == "sans_donnees":
        return analyse_prior(names, data)
    good = True
    ok = data[0, names.index("Chi2")] < 1e4
    system = model.System(spec["mstar"], spec["orbits"])
    for k, orbit in enumerate(spec["orbits"]):
        print(f"  Orbite {k + 1}")
        raw = lambda n: data[k, names.index(n)][ok]  # noqa: E731,B023

        def get(n, raw=raw):
            # En variables universelles, les colonnes a et P contiennent q et Pq
            if univ and n in ("a", "P"):
                factor = 1.0 - raw("e")
                return raw("a") / factor if n == "a" else raw("P") / factor**1.5
            return raw(n)
        true_period = 2.0 * math.pi / math.sqrt(system.sigma[k] / orbit["a"] ** 3)
        good &= report_line("a (au)", get("a"), orbit["a"])
        good &= report_line("P (jours)", get("P"), true_period)
        good &= report_line("e", get("e"), orbit["e"])
        good &= report_line("i (deg)", np.degrees(get("i")), orbit["i"])
        if spec.get("rv") or spec.get("rvp"):
            good &= report_line("Omega (écart)", angle_offset(get("W"), orbit["node"]), 0.0)
            good &= report_line("omega (écart)", angle_offset(get("w"), orbit["w"]), 0.0)
        if len(spec["orbits"]) > 1 or spec.get("rv"):
            good &= report_line("m (Mjup)", get("m") * model.MJUP / model.MJUP, orbit["m_mjup"])
    if len(spec["orbits"]) > 1 or spec.get("rv"):
        good &= report_line("m0 (Msun)", data[0, names.index("m0")][ok], spec["mstar"])
    else:
        mass = (2.0 * np.pi / data[0, names.index("P")][ok]) ** 2 * data[0, names.index("a")][ok] ** 3 / model.SMAS
        good &= report_line("M tot (Msun)", mass, system.sigma[0] / model.SMAS)
    if spec.get("rv"):
        good &= report_line("V0 (km/s)", data[0, names.index("V0")][ok], V0)
    chi = data[0, names.index("Chi2")][ok]
    print(f"  chi2 médian {np.median(chi):.1f}, minimal {chi.min():.1f}")
    return good


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["liste", "lance", "analyse"])
    parser.add_argument("nom", nargs="?")
    parser.add_argument("--univ", action="store_true")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--limite", type=int, default=3600)
    parser.add_argument("--rapide", action="store_true", help="20 000 orbites au lieu de 400 000 (NINFER = 2500)")
    args = parser.parse_args()
    if args.action == "liste":
        print("\n".join(SYSTEMS))
    elif args.action == "lance":
        launch(args.nom, args.univ, args.threads, args.limite, 2500 if args.rapide else None)
    else:
        return 0 if analyse(args.nom, args.univ) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
