# ORACLE et les autres codes d'ajustement orbital

Bilan établi le 2026-11-06 à partir des documentations en ligne et des dépôts des codes. Ces outils évoluent vite : revérifier avant de citer quoi que ce soit dans la thèse. « Non indiqué » signifie que la documentation lue n'en dit rien, pas que la fonction n'existe pas.

## Vue d'ensemble

| | ORACLE | orbitize! | orvara | Octofitter | BINARYS |
|---|---|---|---|---|---|
| Dernière version vue | branche `gaia-astrometrie-absolue` | 3.4.0 (juillet 2026) | 1.1.4 (août 2022) | 9.0.0 (août 2026) | 3.0 (README) |
| Langage | Fortran (OpenMP) + interface Python | Python, Kepler en C (GPU possible) | Python + Cython | Julia | C++ (TMB) piloté depuis R |
| Licence | | BSD-3 | BSD-3 | MIT | non indiquée |
| Échantillonnage | MCMC un paramètre à la fois, 8 chaînes, Gelman-Rubin ; LM au départ | OFTI, MCMC à tempering parallèle (ptemcee), échantillonnage emboîté (dynesty, MultiNest) | MCMC à tempering parallèle (ptemcee) | HMC/NUTS, tempering parallèle non réversible (Pigeons), rejet | descente de gradient avec différentiation automatique ; MCMC possible (tmbstan) |
| e ≥ 1 | **oui** (variables universelles) | non (e dans [0, 1]) | non indiqué, code elliptique | **annoncé** (« hyperbolic orbits »), sans exemple trouvé | non (formules elliptiques) |
| Plusieurs compagnons | oui, coordonnées de Jacobi | Kepler par planète ; N-corps (rebound) en option | Kepler emboîtés (masse intérieure) | Kepler superposés, hiérarchie (Jacobi possible) déclarée ; N-corps en option | non indiqué (un compagnon) |
| Astrométrie relative | oui | oui | oui | oui | oui |
| VR étoile / compagnon | oui / oui | oui / oui | oui / oui | oui / oui | oui / oui |
| Hipparcos (mesures) | non | oui | oui | oui | oui (IAD et transits) |
| Accélérations Hipparcos-Gaia (HGCA) | non | oui | oui | oui | non indiqué |
| Gaia DR4 par époque | en cours | non indiqué | non indiqué | **oui** (`GaiaDR4AstromObs`, pré-publication) | **oui** (utilisé pour BH3) |
| Paramètres linéaires | échantillonnés | échantillonnés | ajustés par moindres carrés à chaque pas | zéro des VR éliminé analytiquement | ajustés |

## Détails par code

### orbitize! (Blunt et al. 2020, 2024)
- **Échantillonneurs :** OFTI (tirage puis rejet, efficace pour un arc court et peu de données), MCMC à tempering parallèle (20 températures et 1000 marcheurs conseillés), échantillonnage emboîté. La documentation prévient que le MCMC « may not work well for multi-modal distributions ».
- **Départ des chaînes :** tiré uniformément dans le prior. C'est l'inverse d'ORACLE, où toutes les chaînes partent de la solution du LM (anomalie B12).
- **Convergence :** inspection visuelle des chaînes et temps d'autocorrélation, sans test formel.
- **Paramétrisations :** `a, e, i, ω, Ω, τ`, ou avec la période, ou position et vitesse cartésiennes (`XYZ`, avec des limites), ou des priors fondés sur les observables (O'Neil et al. 2019).
- **Gaia :** DR2 et EDR3 seulement, et seulement combinés à Hipparcos. Rien sur DR4.

### orvara (Brandt et al. 2021)
- Pas de nouvelle version depuis 2022.
- Ajuste ensemble les VR, l'astrométrie relative, les mesures d'Hipparcos et le catalogue d'accélérations Hipparcos-Gaia.
- **Paramètres :** `√e·sin ω`, `√e·cos ω`, et la **longitude moyenne à l'époque de référence** (proche de l'idée de nos paramètres 6 et 7).
- **Paramètres linéaires** (parallaxe, mouvement propre, zéro des VR) : ajustés par moindres carrés à chaque pas, ce qui rend le code rapide.
- **Plusieurs compagnons :** chaque orbite utilise la masse intérieure, comme des coordonnées de Jacobi.
- Pas de diagnostic de convergence automatique. Parallélisation déclarée médiocre.

### Octofitter (Thompson et al. 2023)
- **Le plus complet et le plus actif** : version 9 en août 2026, documentation mise à jour le 22 octobre 2026.
- **Échantillonneurs :** HMC/NUTS avec gradients automatiques, et tempering parallèle non réversible (Pigeons) pour les distributions à plusieurs solutions. Départ préparé par une approximation variationnelle (Pathfinder). Diagnostics `rhat` et `ess`.
- **Paramétrisations :** au choix, par groupes : taille `a` ou `P` ; forme `(e, ω)`, `(√e cos ω, √e sin ω)` ou `(e cos ω, e sin ω)` ; phase `tp`, anomalie moyenne ou **angle de position `θ` à une époque** ; ou **position et vitesse cartésiennes à une époque** (`x, y, z, vx, vy, vz`).
- **Orbites hyperboliques :** citées comme fonctionnalité, et la version 9 annonce « better hyperbolic orbits and Cartesian state vectors ». Mais la FAQ donne `e` dans [0, 1), et aucun exemple n'a été trouvé. À vérifier en pratique.
- **Gaia DR4 :** type `GaiaDR4AstromObs` pour les données de la pré-publication (12 sources). Modèle : décalages, mouvement propre, parallaxe et jitter. La réduction (regroupement par transit) reste à faire par l'utilisateur. La documentation dit ne pas modeler la part d'erreur commune aux CCD d'un même transit.
- Hipparcos, catalogue d'accélérations, Gaia DR2 et DR3, images, interférométrie, photométrie.
- Annonce un ajustement d'astrométrie « en quelques secondes sur un portable ».

### BINARYS (Leclerc et al. 2023, équipe de Carine Babusiaux)
- **Code :** C++ avec TMB, piloté depuis R. Ajustement par descente de gradient avec différentiation automatique, erreurs par la matrice de covariance (`sdreport`). MCMC possible avec `tmbstan`.
- **Données :** mesures d'Hipparcos (deux réductions), solutions Gaia DR2, EDR3 et DR3 (y compris les solutions d'étoiles non simples de DR3), **mesures Gaia DR4 par époque** (version 3.0, pour BH3), astrométrie relative, VR des deux composantes.
- **Modèle :** photocentre avec la fraction de luminosité `β` ajustée, accélération de perspective (version 2.1). Pas de terme de couleur dans sa documentation, alors que la réduction `GaiaTools.R` reçue d'Octave en contient un.
- Un seul compagnon, orbites elliptiques.
- Licence non indiquée. Le dépôt GitLab est dans un groupe « ipag-public », mais Octave présente `GaiaTools.R` comme la partie privée.

### kepmodel (Delisle et Ségransan 2022)
- Utilisé par le notebook ESA de BH3. Site de documentation inaccessible le 2026-11-06 : rien de vérifié.

## Ce que ça change pour ORACLE

1. **Les orbites non liées ne sont plus un argument exclusif.** Octofitter annonce les orbites hyperboliques et une paramétrisation cartésienne. Il faut tester s'il donne un résultat juste sur un cas comme `arc_parabole` avant d'affirmer quoi que ce soit dans un sens ou dans l'autre. orbitize!, orvara et BINARYS restent limités aux ellipses.
2. **ORACLE ressemble déjà aux autres sur deux points.** La longitude à l'époque de référence (nos paramètres 6 et 7) se retrouve dans orvara (longitude moyenne) et dans Octofitter (angle `θ`). L'échantillonnage en position et vitesse cartésiennes, envisagé ici pour plus tard, existe dans orbitize! et Octofitter.
3. **L'échantillonneur est le point le plus faible.** orbitize! fait partir ses chaînes du prior, Octofitter d'une approximation variationnelle, et les deux proposent le tempering parallèle. ORACLE fait partir toutes les chaînes du même point et peut se déclarer convergé trop tôt (B12).
4. **Pour Gaia DR4, deux références existent déjà** : BINARYS (déjà utilisé sur BH3) et Octofitter (`GaiaDR4AstromObs`). On peut comparer ORACLE aux deux sur BH3. Aucun des deux ne modélise les erreurs communes aux CCD d'un transit : c'est la question 4 posée à Carine.
5. **Ce qui reste propre à ORACLE**, d'après cette lecture :
   - les variables universelles avec lambda, testées sur des orbites quasi paraboliques ;
   - plusieurs orbites en coordonnées de Jacobi avec des priors sur des combinaisons de masses, et toutes les données ensemble (relative, VR des deux corps, Gaia) dans un même ajustement ;
   - l'interface graphique.

   Octofitter couvre une grande partie de cette liste. La niche la plus nette reste les orbites peu couvertes et très excentriques, à condition de le montrer par une comparaison sur les mêmes données.

## Comparaisons à faire

- BH3 : ORACLE contre le notebook ESA (kepmodel), BINARYS et Octofitter.
- **`arc_parabole` (Tests/synthetic.py) : ORACLE lambda contre Octofitter** (décidé avec Antoine le 2026-11-06). Lancer Octofitter dans deux paramétrisations : éléments (`a` ou `P`, `e`, ...) et position-vitesse cartésiennes à une époque, avec un prior qui laisse `e` dépasser 1. Questions : Octofitter accepte-t-il vraiment `e ≥ 1` pendant l'échantillonnage (sa FAQ dit [0, 1)) ? Retrouve-t-il la vérité (e = 0,98) sans biais, avec une distribution de `e` comparable à lambda (0,94 à 1,04 à 95 %) ? Installer Julia (une à deux heures). Les données sont dans `Tests/build_synthetic/peri/parab_lambda/arc_parabole_univ/data.dat`, ou se régénèrent avec `synthetic.py`.
- Un compagnon imagé connu (β Pic b par exemple) contre orbitize! et orvara.

## Sources lues le 2026-11-06

- orbitize! : https://orbitize.readthedocs.io/en/latest/ (pages sampler, basis, priors, kepler, nbody, hipparcos, gaia, driver, system, manual, tutoriels MCMC et multiplanète), https://github.com/sblunt/orbitize
- orvara : https://github.com/t-brandt/orvara (README, CHANGES.md, setup.py, main.py, orbit.pyx), arXiv:2109.06761. Le site readthedocs renvoyait une erreur 404.
- Octofitter : https://sefffal.github.io/Octofitter.jl/dev/ (pages samplers, thiele-innes, gaia-dr4-prerelease, hipparcos, g23h, rv-multi-planet, orbit-computation, v9-migration, api, faq, priors, pma), https://github.com/sefffal/Octofitter.jl et ses releases
- BINARYS : https://gricad-gitlab.univ-grenoble-alpes.fr/ipag-public/gaia/binarys (README, DetailedDescription.md), arXiv:2209.04210 (version 1), https://www.cosmos.esa.int/web/gaia/gaia-bhthree
- kepmodel : https://obswww.unige.ch/~delisle/kepmodel/doc/ inaccessible
