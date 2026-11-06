#!/usr/bin/env python3
"""Refait les copies corrigees et la variante lambda, sans toucher au depot.

    python3 Tests/parametrisation/prepare.py           # copies seules
    python3 Tests/parametrisation/prepare.py --verifie # plus la verification de la variante lambda

Copies creees dans Tests/build_synthetic/variantes/ :
    B5B6/sub           kepu.f corrige de B5 et B6
    B5B6B10/sub        idem, plus B10 (meme prior que COMPUTE_MAP)
    lambda/main, sub   base B5B6B10, parametres 6 et 7 = (cos, sin)(lambda)/Pq

Voir Tests/parametrisation/README.md.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTS)
CODE = os.path.join(ROOT, "Code", "mcmc")
OUT = os.path.join(TESTS, "build_synthetic", "variantes")
IGNORE = shutil.ignore_patterns("~*", "*.mod", "*.o")


def copy(source, target):
    if os.path.isdir(target):
        shutil.rmtree(target)
    shutil.copytree(source, target, ignore=IGNORE)


def run(*args):
    subprocess.run([sys.executable, *args], check=True)


def prepare():
    for name, extra in (("B5B6", []), ("B5B6B10", ["B10"])):
        sub = os.path.join(OUT, name, "sub")
        copy(os.path.join(CODE, "sub"), sub)
        run(os.path.join(HERE, "patch_B5B6.py"), sub, *extra)
    lam = os.path.join(OUT, "lambda")
    copy(os.path.join(OUT, "B5B6B10", "sub"), os.path.join(lam, "sub"))
    copy(os.path.join(CODE, "main"), os.path.join(lam, "main"))
    run(os.path.join(HERE, "patch_lambda.py"), os.path.join(lam, "sub"), os.path.join(lam, "main"))
    return lam


def verify(lam):
    """Compile le programme de test sur la variante lambda et lance check_lambda.py."""
    sub = os.path.join(lam, "sub")
    build = os.path.join(lam, "test")
    os.makedirs(build, exist_ok=True)
    base = ["gfortran", "-w", "-O3", "-J", build, "-I", build, "-I", sub]
    objects = []
    sources = [os.path.join(HERE, "driver_lambda.f")]
    sources += [os.path.join(sub, s) for s in ("utils.f", "mcmc.f", "mrqfit.f", "io.f", "kepu.f")]
    sources += [os.path.join(HERE, "longitude_lambda.f")]
    for source in sources:
        obj = os.path.join(build, os.path.basename(source).replace(".f", ".o"))
        subprocess.run([*base, "-c", source, "-o", obj], check=True, cwd=build)
        objects.append(obj)
    exe = os.path.join(build, "driver")
    subprocess.run(["gfortran", "-w", "-O3", *objects, "-o", exe], check=True, cwd=build)
    run(os.path.join(HERE, "check_lambda.py"), exe)


if __name__ == "__main__":
    lam = prepare()
    print("Copies prêtes dans", OUT)
    if "--verifie" in sys.argv:
        verify(lam)
