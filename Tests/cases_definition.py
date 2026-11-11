"""Definition des cas de test : vraie solution et point d'evaluation.

Les donnees sont tirees une fois pour toutes (graine fixe) puis enregistrees
dans Tests/cases/. Les routines sont evaluees a un point volontairement decale
de la vraie solution, pour que le gradient du chi2 ne soit pas nul.
"""

INNER = dict(m_mjup=10.0, a=2.0, e=0.3, i=60.0, node=30.0, w=70.0, tp=56000.0)
OUTER = dict(m_mjup=5.0, a=8.0, e=0.1, i=65.0, node=40.0, w=200.0, tp=57000.0)
# Quasi parabolique : passage au periastre pendant les mesures (point d'evaluation e = 0.99)
PARAB = dict(m_mjup=10.0, a=3.0, e=0.97, i=50.0, node=30.0, w=70.0, tp=56200.0)


def offset(orbit):
    """Point d'evaluation : la vraie orbite legerement decalee."""
    return dict(
        m_mjup=orbit["m_mjup"] * 1.10,
        a=orbit["a"] * 1.02,
        e=orbit["e"] + 0.02,
        i=orbit["i"] + 2.0,
        node=orbit["node"] + 3.0,
        w=orbit["w"] + 4.0,
        tp=orbit["tp"] + 20.0,
    )


# types : relative, vitesses de l'etoile, vitesses de la planete
# position_format : 1 = (Dec, RA), 2 = (Sep, PA)
CASES = {
    "rel": dict(orbits=[INNER], rel=True),
    "rel_corr": dict(orbits=[INNER], rel=True, corr=0.3),
    "rel_seppa": dict(orbits=[INNER], rel=True, corr=0.2, position_format=2),
    "rel_2orb": dict(orbits=[INNER, OUTER], rel=True),
    "rel_rv": dict(orbits=[INNER], rel=True, rv=True),
    "rel_rv_jit": dict(orbits=[INNER], rel=True, rv=True, jitter=0.02),
    "rel_rv_prior": dict(orbits=[INNER], rel=True, rv=True, prior="1 0\n1\n1.0 0.1 ms"),
    "rel_rv_2orb_jit": dict(orbits=[INNER, OUTER], rel=True, rv=True, jitter=0.02),
    "rel_rvp": dict(orbits=[INNER], rel=True, rvp=True),
    "rel_rv_rvp": dict(orbits=[INNER], rel=True, rv=True, rvp=True),
    "rv_seul": dict(orbits=[INNER], rv=True),
    # Orbite quasi parabolique (2026-11-11), classique et universel, borne sur e relevee.
    # Doit rester le dernier cas : les donnees sont tirees dans l'ordre de la liste.
    "rel_rv_jit_parab": dict(orbits=[PARAB], rel=True, rv=True, jitter=0.02, emax=0.995),
}

MSTAR_TRUE = 1.0
MSTAR_POINT = 1.02
DIST_PC = 50.0
V0_TRUE = 5.0
V0_POINT = 5.01
T0 = 56000.0
SIG_POS_MAS = 1.0
SIG_RV_KMS = 0.01
SIG_RVP_KMS = 0.5
