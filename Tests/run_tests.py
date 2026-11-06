#!/usr/bin/env python3
"""Tests de verification d'ORACLE. Usage : make test

Chaque ligne affiche OK, ECHEC ou CONNU (anomalie deja listee dans
MODIFICATIONS.md, section "Anomalies connues"). Le code de sortie vaut 0
seulement s'il n'y a aucun ECHEC. Voir Tests/README.md.

    --update-reference  refait la reference (seulement apres un changement
                        volontaire des resultats, a justifier dans le journal)
"""

import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(ROOT, "Tests")
sys.path.insert(0, TESTS)
import jacobian  # noqa: E402
import model  # noqa: E402

# ORACLE_SUB permet d'essayer une correction sur une copie, hors du depot
SUB = os.environ.get("ORACLE_SUB", os.path.join(ROOT, "Code", "mcmc", "sub"))
MAIN = os.path.join(ROOT, "Code", "mcmc", "main")
BUILD = os.environ.get("ORACLE_BUILD", os.path.join(TESTS, "build"))
REFERENCE = os.path.join(TESTS, "reference", "reference.json")
JOURNAL = os.path.join(ROOT, "MODIFICATIONS.md")
JOURNAL_BASE = "bcced4b"  # Dernier commit avant le debut du journal
VERSIONS = {"ellip": (1, "kepellip.f"), "univ": (2, "kepu.f")}
COMMON_SOURCES = ["utils.f", "mcmc.f", "mrqfit.f", "io.f"]
GFORTRAN = os.environ.get("COMPILF", "gfortran")

# Anomalies connues (voir MODIFICATIONS.md). Tant qu'elles existent, le test
# concerne est marque CONNU. S'il passe, il est marque ECHEC pour obliger a
# retirer la ligne ici et a mettre le journal a jour.
# B5 et B6 corrigees par M005 : plus d'ecart connu sur le chi2 universel.
KNOWN_CHI2_UNIV = {}
KNOWN_GRADIENT = {}  # B7 corrigee par M009


def known_memory(stderr):
    """Identifie un debordement de tableau deja liste."""
    # B3 corrigee par M008, B4 par M006 : plus aucun debordement connu
    return None


TOL_REFERENCE = 1e-10
TOL_INDEPENDENT = 1e-7
TOL_GRADIENT = 1e-4
TOL_PRIOR = 1e-4


class Report:
    def __init__(self):
        self.counts = {"OK": 0, "ECHEC": 0, "CONNU": 0}

    def line(self, status, name, detail=""):
        self.counts[status] += 1
        if status != "OK" or VERBOSE:
            print(f"  {status:<6} {name}" + (f" : {detail}" if detail else ""))

    def check(self, condition, name, detail="", known=None):
        if known:
            if condition:
                self.line("ECHEC", name, f"passe alors que {known} est listee : mettre a jour KNOWN et le journal")
            else:
                self.line("CONNU", name, f"{known} ; {detail}")
        else:
            self.line("OK" if condition else "ECHEC", name, detail)

    def section(self, title):
        print(f"\n{title}")


VERBOSE = "-v" in sys.argv


# ---------------------------------------------------------------- compilation

def quiet(command):
    """Lance une compilation et n'affiche sa sortie qu'en cas d'erreur.

    Depuis BUILD : gfortran cherche les modules (.mod) d'abord dans le dossier courant,
    et un ancien data.mod peut trainer a la racine du depot.
    """
    os.makedirs(BUILD, exist_ok=True)
    result = subprocess.run(command, capture_output=True, text=True, cwd=BUILD)
    if result.returncode:
        print(result.stdout + result.stderr)
        raise SystemExit("Compilation impossible : " + " ".join(command))


def compile_version(version, flags, folder):
    """Compile le programme de test et les routines courantes dans folder."""
    os.makedirs(folder, exist_ok=True)
    method, kep = VERSIONS[version]
    objects = []
    base = [GFORTRAN, "-w", *flags, "-J", folder, "-I", folder, "-I", SUB]
    driver_obj = os.path.join(folder, "driver.o")
    quiet([*base, "-c", os.path.join(TESTS, "driver.f"), "-o", driver_obj])
    objects.append(driver_obj)
    for source in COMMON_SOURCES + [kep]:
        obj = os.path.join(folder, source.replace(".f", ".o"))
        quiet([*base, "-c", os.path.join(SUB, source), "-o", obj])
        objects.append(obj)
    longitude = os.path.join(folder, "longitude.o")
    quiet([*base, "-c", os.path.join(TESTS, f"longitude_{version}.f"), "-o", longitude])
    objects.append(longitude)
    exe = os.path.join(folder, "driver")
    quiet([GFORTRAN, "-w", *flags, *objects, "-o", exe])
    return exe, objects


def build_all():
    """Compile deux fois : comme le Makefile (-O3) et avec controles (-O0)."""
    if os.path.isdir(BUILD):
        shutil.rmtree(BUILD)
    exes = {}
    for version in VERSIONS:
        exes[version, "O3"] = compile_version(version, ["-O3"], os.path.join(BUILD, version, "O3"))[0]
        exes[version, "check"] = compile_version(
            version, ["-O0", "-g", "-fcheck=all", "-fbacktrace"], os.path.join(BUILD, version, "check")
        )[0]
    return exes


# ---------------------------------------------------------------- execution

def run_driver(exe, version, case):
    """Lance le programme de test sur un cas et lit les lignes TEST."""
    folder = os.path.join(TESTS, "cases", case)
    method = VERSIONS[version][0]
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(os.path.join(folder, "data.dat"), tmp)
        with open(os.path.join(folder, f"input_{version}.txt")) as stdin:
            result = subprocess.run(
                [exe, str(method)], stdin=stdin, cwd=tmp, capture_output=True, text=True, timeout=120
            )
    values = {"returncode": result.returncode, "stderr": result.stderr.strip()[-300:]}
    for line in result.stdout.splitlines():
        if not line.startswith("TEST "):
            continue
        parts = line.split()
        if parts[1] == "META":
            values.update(npar=int(parts[2]), npla=int(parts[3]), multipla=parts[5] == "T")
        elif len(parts) == 4:
            values.setdefault(parts[1], {})[int(parts[2])] = float(parts[3])
        else:
            values[parts[1]] = float(parts[2])
    return values


def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)


# ---------------------------------------------------------------- tests

def test_journal(report):
    """Toute ligne changee depuis JOURNAL_BASE doit etre reperee et journalisee."""
    report.section("Journal et reperes dans le code")
    journal = open(JOURNAL, encoding="utf-8").read()
    entries = re.findall(r"^## (M\d{3})\b.*$", journal, flags=re.M)
    without_marker = set(re.findall(r"^## (M\d{3})\b.*\(sans repère\)\s*$", journal, flags=re.M))
    report.check(len(entries) == len(set(entries)), "numeros du journal uniques")

    marker = re.compile(r"(>>>|<<<) \[MOD (M\d{3})\]")
    git = ["git", "-C", ROOT]
    tracked = subprocess.run(
        [*git, "ls-files", "--cached", "--others", "--exclude-standard", "Code/mcmc", "Interface", "Makefile"],
        capture_output=True, text=True, check=True,
    ).stdout.split("\n")
    files = [f for f in tracked if f.endswith((".f", ".py")) or f == "Makefile"]
    seen = set()
    for name in files:
        path = os.path.join(ROOT, name)
        if not os.path.isfile(path):
            continue
        with open(path, encoding="utf-8", errors="replace", newline="") as f:
            lines = f.read().split("\n")
        # Lignes couvertes par un bloc de reperes
        covered, current, ok_blocks = set(), None, True
        for number, text in enumerate(lines, 1):
            match = marker.search(text)
            if match:
                kind, ident = match.groups()
                seen.add(ident)
                covered.add(number)
                if kind == ">>>":
                    ok_blocks = ok_blocks and current is None
                    current = ident
                else:
                    ok_blocks = ok_blocks and current == ident
                    current = None
            elif current:
                covered.add(number)
        report.check(ok_blocks and current is None, f"blocs bien fermes : {name}")
        # Lignes changees par rapport a la base
        known = subprocess.run([*git, "cat-file", "-e", f"{JOURNAL_BASE}:{name}"], capture_output=True).returncode == 0
        if not known:
            changed = set(range(1, len(lines)))
        else:
            diff = subprocess.run(
                [*git, "diff", "-U0", JOURNAL_BASE, "--", name], capture_output=True, text=True, errors="replace"
            ).stdout
            changed = set()
            for match in re.finditer(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", diff, flags=re.M):
                start, count = int(match.group(1)), int(match.group(2) or 1)
                changed |= set(range(start, start + count)) if count else {start}
        missing = sorted(changed - covered)
        report.check(not missing, f"changements repérés : {name}", f"lignes hors bloc : {missing[:12]}")
    for ident in sorted(seen - set(entries)):
        report.check(False, f"repère {ident} absent du journal")
    for ident in sorted(set(entries) - without_marker - seen):
        report.check(False, f"entrée {ident} du journal sans repère dans le code")


def test_cases(report, exes, update):
    results = {}
    report.section("Exécution (compilation -O3, puis avec contrôles -O0 -fcheck=all)")
    for case in sorted(os.listdir(os.path.join(TESTS, "cases"))):
        if not os.path.isdir(os.path.join(TESTS, "cases", case)):
            continue
        for version in VERSIONS:
            label = f"{case}/{version}"
            res = run_driver(exes[version, "O3"], version, case)
            results[label] = res
            report.check(res["returncode"] == 0 and "CHI2" in res, f"exécution {label}", res["stderr"])
            chk = run_driver(exes[version, "check"], version, case)
            stderr = " ".join(chk["stderr"].split())[:110]
            report.check(chk["returncode"] == 0, f"contrôles mémoire {label}", stderr, known=known_memory(chk["stderr"]))
            if "CHI2" in res and "CHI2" in chk:
                report.check(close(res["CHI2"], chk["CHI2"], TOL_REFERENCE), f"-O3 = -O0 {label}",
                             f"{res['CHI2']} / {chk['CHI2']}")

    usable = {k: v for k, v in results.items() if "CHI2" in v}

    report.section("Référence : chi2 et MAP inchangés depuis la référence enregistrée")
    current = {k: {"CHI2": v["CHI2"], "MAP": v["MAP"], "MERIT": v["MERIT"]} for k, v in usable.items()}
    if update:
        os.makedirs(os.path.dirname(REFERENCE), exist_ok=True)
        with open(REFERENCE, "w") as f:
            json.dump(current, f, indent=1, sort_keys=True)
        print(f"  Référence réécrite : {len(current)} entrées")
    with open(REFERENCE) as f:
        reference = json.load(f)
    for label in sorted(set(reference) | set(current)):
        if label not in current or label not in reference:
            report.check(False, f"référence {label}", "cas absent d'un côté")
            continue
        for key in ("CHI2", "MAP", "MERIT"):
            report.check(close(current[label][key], reference[label][key], TOL_REFERENCE), f"référence {label} {key}",
                         f"{current[label][key]!r} au lieu de {reference[label][key]!r}")

    report.section("Modèle Python indépendant : même chi2 que le Fortran")
    for label, res in sorted(usable.items()):
        case = label.split("/")[0]
        spec = json.load(open(os.path.join(TESTS, "cases", case, "point.json")))
        if case == "rel_seppa":
            continue
        value = model.chi2(spec["point"], spec["data"])
        known = KNOWN_CHI2_UNIV.get(case) if label.endswith("/univ") else None
        report.check(close(value, res["CHI2"], TOL_INDEPENDENT), f"indépendant {label}",
                     f"Python {value:.9f} / Fortran {res['CHI2']:.9f}", known=known)

    report.section("Variables classiques et universelles : même chi2 pour la même orbite")
    for label, res in sorted(usable.items()):
        case, version = label.split("/")
        if version == "ellip" and f"{case}/univ" in usable:
            other = usable[f"{case}/univ"]["CHI2"]
            report.check(close(res["CHI2"], other, TOL_INDEPENDENT), f"classique = universel {case}",
                         f"{res['CHI2']:.9f} / {other:.9f}", known=KNOWN_CHI2_UNIV.get(case))

    report.section("Cohérence interne des routines")
    for label, res in sorted(usable.items()):
        path = res["MERIT"] - res["FMAP"]
        known = KNOWN_CHI2_UNIV.get(label.split("/")[0]) if label.endswith("/univ") else None
        report.check(close(path, res["CHI2"], 1e-9), f"CHI2SEUL = MRQCOF {label}", f"{res['CHI2']} / {path}", known=known)
        report.check(close(-2.0 * res["MAP"], res["MERIT"], 1e-6), f"MAP = chi2 + Fmap {label}",
                     f"{-2.0 * res['MAP']} / {res['MERIT']}", known=known)

    report.section("Rapport des priors du MCMC (RATIO_PRIORS) contre jacobien numérique et COMPUTE_MAP")
    for label, res in sorted(usable.items()):
        case, version = label.split("/")
        jitter = "jit" in case
        npar = res["npar"]
        base = np.array([res["PM"][i] for i in range(1, npar + 1)])
        lnj0 = jacobian.ln_jacobian(base, res["npla"], version, res["multipla"], jitter)
        bad_exact, bad_known = [], []
        for k in range(1, 5):
            if res["RAPQ"][k] <= 0:  # Masse negative au point deplace : point refuse, rien a comparer
                continue
            point = np.array([res["PT"][1000 * k + i] for i in range(1, npar + 1)])
            lnj = jacobian.ln_jacobian(point, res["npla"], version, res["multipla"], jitter)
            expected = res["LNPRIOR"][k] - res["LNPRIOR0"] - (lnj - lnj0)
            value = math.log(res["RAPQ"][k])
            gap = jacobian.known_prior_gap(point, base, res["npla"], version, res["multipla"], case == "rel_rv_prior")
            if abs(value - expected) > TOL_PRIOR * 0.1:
                bad_exact.append(f"point {k}: {value:.6f} / {expected:.6f}")
            if abs(value - expected - gap) > TOL_PRIOR * 0.1:
                bad_known.append(f"point {k}: écart {value - expected:.6f} / identifié {gap:.6f}")
        codes = "+".join(c for c, on in (("B10", version == "univ"), ("B11", res["multipla"])) if on)
        report.check(not bad_known, f"rapport des priors, écart expliqué {label}", "; ".join(bad_known[:2]))
        report.check(not bad_exact, f"rapport des priors {label}", "; ".join(bad_exact[:1]), known=codes or None)

    report.section("Dérivées du Levenberg-Marquardt contre différences finies")
    for label, res in sorted(usable.items()):
        bad = []
        for i in res["BETA"]:
            analytic, numeric = res["BETA"][i], res["BETAFD"][i]
            # Chaque parametre a sa tolerance : ecart relatif, plus le bruit d'arrondi
            # des differences finies (arrondi du chi2 divise par le pas).
            step = 1e-6 * max(abs(res["P"][i]), 1e-2)
            noise = 1e-11 * abs(res["MERIT"]) / step
            if abs(analytic - numeric) > TOL_GRADIENT * abs(numeric) + noise:
                bad.append(f"p{i}: {analytic:.6g} / {numeric:.6g}")
        report.check(not bad, f"gradient {label}", "; ".join(bad[:5]), known=KNOWN_GRADIENT.get(label))

    report.section("Reconstruction des éléments (MRQELEMENTS)")
    for label, res in sorted(usable.items()):
        if not res["multipla"]:
            continue
        mstar = math.exp(res["PLAST"])
        report.check(close(res["MSTAR"], mstar, 1e-9), f"masse de l'étoile {label}",
                     f"{res['MSTAR'] / model.SMAS:.6f} au lieu de {mstar / model.SMAS:.6f} Msun")  # B1 : M006
        expected = res["COVLAST"] * mstar
        report.check(close(res["DMSTAR"], expected, 1e-9), f"erreur sur la masse {label}",
                     f"{res['DMSTAR']:.4g} au lieu de {expected:.4g}")  # B2 : M006


def test_priors(report, exes):
    """Priors log-normal, lineaire et fixe, ecrits comme le fait l'interface."""
    report.section("Lecture des priors sur les masses")
    folder = os.path.join(TESTS, "cases", "rel_rv")
    with open(os.path.join(folder, "input_ellip.txt")) as f:
        lines = f.read().split("\n")
    position = lines.index("100000000") + 1
    for kind, values in (("2", "1.0 0.1 ms"), ("3", "0.5 2.0 ms"), ("4", "1.0 ms")):
        stdin = lines[:position] + ["1", "1 0", kind, values] + lines[position + 1:]
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(os.path.join(folder, "data.dat"), tmp)
            result = subprocess.run([exes["ellip", "O3"], "1"], input="\n".join(stdin), cwd=tmp,
                                    capture_output=True, text=True, timeout=120)
        message = " ".join(result.stderr.split())[:80]
        report.check(result.returncode == 0, f"prior de type {kind}", message)  # B8 : M007, lecture seule


def test_prior_mstar(report, exes):
    """Prior sur m0 sans masse d'etoile ajustee : il doit agir ou etre refuse."""
    report.section("Prior sur la masse de l'étoile, une orbite sans vitesses de l'étoile")
    folder = os.path.join(TESTS, "cases", "rel")
    with open(os.path.join(folder, "input_ellip.txt")) as f:
        lines = f.read().split("\n")
    position = lines.index("100000000") + 1
    stdin = lines[:position] + ["1", "1 0", "1", "1.0 0.1 ms"] + lines[position + 1:]
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(os.path.join(folder, "data.dat"), tmp)
        result = subprocess.run([exes["ellip", "O3"], "1"], input="\n".join(stdin), cwd=tmp,
                                capture_output=True, text=True, timeout=120)
    found = re.search(r"TEST MAP\s+(\S+)", result.stdout)
    without = run_driver(exes["ellip", "O3"], "ellip", "rel")["MAP"]
    # Masse totale du point : 1.03 Msun, donc un terme de prior faible s'il agit.
    # S'il est evalue sur une masse nulle, il vaut -0.5*(1.0/0.1)^2 = -50.
    term = float(found.group(1)) - without if found else float("nan")
    report.check(result.returncode != 0 or abs(term) < 10.0, "prior sur m0 pris en compte",
                 f"terme de prior {term:.3f}, soit la valeur pour une masse nulle", known="B9")

    # Cas ou le prior doit agir : la valeur attendue est le terme gaussien seul.
    mjup = 1.0 / 1047.355
    checks = (
        ("rel", "1 1", 1.02 + 11.0 * mjup, "masse totale, une orbite"),
        ("rel_rv", "1 0", 1.02, "masse de l'étoile, avec vitesses de l'étoile"),
        ("rel_2orb", "1 0 0", 1.02, "masse de l'étoile, deux orbites"),
        ("rel_2orb", "1 1 1", 1.02 + 16.5 * mjup, "masse totale, deux orbites"),
    )
    for case, coefficients, mass, title in checks:
        folder = os.path.join(TESTS, "cases", case)
        with open(os.path.join(folder, "input_ellip.txt")) as f:
            lines = f.read().split("\n")
        position = lines.index("100000000") + 1
        stdin = lines[:position] + ["1", coefficients, "1", "1.0 0.1 ms"] + lines[position + 1:]
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(os.path.join(folder, "data.dat"), tmp)
            result = subprocess.run([exes["ellip", "O3"], "1"], input="\n".join(stdin), cwd=tmp,
                                    capture_output=True, text=True, timeout=120)
        found = re.search(r"TEST MAP\s+(\S+)", result.stdout)
        without = run_driver(exes["ellip", "O3"], "ellip", case)["MAP"]
        term = float(found.group(1)) - without if found else float("nan")
        expected = -0.5 * ((mass - 1.0) / 0.1) ** 2
        report.check(abs(term - expected) < 1e-8, f"prior actif : {title}", f"{term:.6f} au lieu de {expected:.6f}")


def test_smoke(report):
    """Levenberg-Marquardt complet sur GG Tau, avec le vrai programme."""
    report.section("Essai réel : GG Tau en Levenberg-Marquardt")
    folder = os.path.join(BUILD, "smoke")
    os.makedirs(folder, exist_ok=True)
    base = [GFORTRAN, "-w", "-O3", "-J", folder, "-I", folder, "-I", SUB]
    main_obj = os.path.join(folder, "main.o")
    quiet([*base, "-c", os.path.join(MAIN, "astrom_mcmco.f"), "-o", main_obj])
    objects = [main_obj]
    for source in COMMON_SOURCES + ["kepellip.f"]:
        obj = os.path.join(folder, source.replace(".f", ".o"))
        quiet([*base, "-c", os.path.join(SUB, source), "-o", obj])
        objects.append(obj)
    exe = os.path.join(folder, "astrom_mcmco")
    quiet([GFORTRAN, "-w", "-O3", *objects, "-o", exe])
    stdin = "\n".join([
        "4", "1 0 0 0", "1", "data.dat", "1d-7", "1 1 1", "145.00 pc", "0.8", "sim.dat", "dump.dat",
        "10000000", "0", "56583", "400 mj 4 0.5 40 110 290 57754", "0.00 10000.00", "0.00 10.00", "0.00 0.90", "",
    ])
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(os.path.join(ROOT, "Example", "GGTauAb12.tar")) as tar:
            member = tar.extractfile("GGTauAb12/data.dat")
            with open(os.path.join(tmp, "data.dat"), "wb") as f:
                f.write(member.read())
        result = subprocess.run([exe], input=stdin, cwd=tmp, capture_output=True, text=True, timeout=120)
    found = re.findall(r"Chi2\+Fmap =\s*([0-9.]+)", result.stdout)
    value = float(found[-1]) if found else float("nan")
    report.check(result.returncode == 0 and abs(value - 25.0651) < 1e-2, "GG Tau : chi2 + Fmap = 25.065",
                 f"obtenu {value}")


def main():
    update = "--update-reference" in sys.argv
    report = Report()
    print("Compilation des routines courantes...")
    exes = build_all()
    test_journal(report)
    test_cases(report, exes, update)
    test_priors(report, exes)
    test_prior_mstar(report, exes)
    test_smoke(report)
    c = report.counts
    print(f"\nBilan : {c['OK']} OK, {c['CONNU']} anomalies connues, {c['ECHEC']} ÉCHEC")
    return 1 if c["ECHEC"] else 0


if __name__ == "__main__":
    sys.exit(main())
