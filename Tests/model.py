"""Modele orbital ecrit independamment du Fortran.

Sert a generer les donnees simulees et a recalculer le chi2 sans utiliser
aucune routine de Code/mcmc. Seules les constantes viennent de commons.f.

Unites : jour, au, masse telle que G = 1 (comme le Fortran).
"""

import numpy as np

KG = 0.01720209895
SMAS = KG * KG
MJUP = SMAS / 1047.355
AU_KM = 1.495978707e8
DAY = 86400.0
MPS = DAY / AU_KM  # 1 km/s en au/jour


def kepler(mean_anomaly, e):
    """Resout M = u - e sin u par Newton."""
    mean_anomaly = np.array(mean_anomaly, dtype=float)
    # Depart de Danby (1987) : converge aussi pres de e = 1 (cas quasi paraboliques)
    u = mean_anomaly + 0.85 * e * np.sign(np.sin(mean_anomaly))
    for _ in range(100):
        u = u - (u - e * np.sin(u) - mean_anomaly) / (1.0 - e * np.cos(u))
    return u


def orbit_state(orbit, sigma, t):
    """Position (x vers le nord, y vers l'est) et vitesse le long de la visee.

    orbit : dictionnaire d'elements (a en au, angles en degres, tp en jours).
    sigma : masse dynamique de l'orbite, en unites du code.
    """
    a, e, tp = orbit["a"], orbit["e"], orbit["tp"]
    i, node, w = np.radians([orbit["i"], orbit["node"], orbit["w"]])
    n = np.sqrt(sigma / a**3)
    u = kepler(n * (np.asarray(t, dtype=float) - tp), e)
    cu, su = np.cos(u), np.sin(u)
    exq = np.sqrt(1.0 - e * e)
    rx = a * (cu - e)
    ry = a * exq * su
    e1x = np.cos(w) * np.cos(node) - np.sin(w) * np.cos(i) * np.sin(node)
    e1y = np.cos(w) * np.sin(node) + np.sin(w) * np.cos(i) * np.cos(node)
    e2x = -np.sin(w) * np.cos(node) - np.cos(w) * np.cos(i) * np.sin(node)
    e2y = -np.sin(w) * np.sin(node) + np.cos(w) * np.cos(i) * np.cos(node)
    x = rx * e1x + ry * e2x
    y = rx * e1y + ry * e2y
    # z = r sin(w + v) sin i, vers l'observateur ; vz = dz/dt
    udot = n / (1.0 - e * cu)
    vz = a * np.sin(i) * udot * (-np.sin(w) * su + exq * np.cos(w) * cu)
    return x, y, vz


class System:
    """Systeme hierarchique en coordonnees de Jacobi."""

    def __init__(self, mstar_msun, orbits):
        self.orbits = orbits
        self.sigma = []
        self.mfrac = []
        total = mstar_msun * SMAS
        for orbit in orbits:
            mass = orbit["m_mjup"] * MJUP
            total += mass
            self.sigma.append(total)
            self.mfrac.append(mass / total)

    def relative_position(self, ipla, t):
        """Position de l'orbite ipla (0..) par rapport a l'etoile, en au."""
        x, y, _ = orbit_state(self.orbits[ipla], self.sigma[ipla], t)
        x, y = np.array(x), np.array(y)
        for j in range(ipla):
            xj, yj, _ = orbit_state(self.orbits[j], self.sigma[j], t)
            x = x + self.mfrac[j] * xj
            y = y + self.mfrac[j] * yj
        return x, y

    def star_rv(self, v0_kms, t):
        """Vitesse radiale de l'etoile, en km/s."""
        v = 0.0
        for j, orbit in enumerate(self.orbits):
            _, _, vz = orbit_state(orbit, self.sigma[j], t)
            v = v + self.mfrac[j] * vz
        return v0_kms + v / MPS

    def planet_rv(self, ipla, t):
        """Vitesse radiale de l'orbite ipla par rapport a l'etoile, en km/s."""
        _, _, v = orbit_state(self.orbits[ipla], self.sigma[ipla], t)
        v = -np.array(v)
        for j in range(ipla):
            _, _, vz = orbit_state(self.orbits[j], self.sigma[j], t)
            v = v - self.mfrac[j] * vz
        return v / MPS


def chi2(point, data):
    """Chi2 d'un jeu de parametres sur des donnees en mas et km/s."""
    system = System(point["mstar"], point["orbits"])
    plx = 1000.0 / point["dist_pc"]
    total = 0.0
    for row in data.get("rel", []):
        ipla, t, xo, yo, sx, sy, rho = row
        x, y = system.relative_position(int(ipla) - 1, t)
        dx = (xo - x * plx) / sx
        dy = (yo - y * plx) / sy
        total += (dx * dx + dy * dy - 2 * rho * dx * dy) / (1 - rho * rho)
    jitter = point.get("jitter", 0.0)
    for t, v, sv in data.get("rv", []):
        total += (v - system.star_rv(point["v0"], t)) ** 2 / (sv * sv + jitter * jitter)
    for ipla, t, v, sv in data.get("rvp", []):
        total += (v - system.planet_rv(int(ipla) - 1, t)) ** 2 / (sv * sv)
    return float(total)
