# Idées à voir plus tard : ajustement, analyse des solutions et dynamique

Notes de discussion, novembre 2026. Rien n'est décidé. À reprendre après les corrections, F1 (vitesses absolues) et Gaia DR4.

## Le point fort visé

Faire le lien entre l'ajustement et la dynamique (idée d'Antoine). Les fitters (orbitize!, orvara, Octofitter, kima, BINARYS) s'arrêtent à la distribution a posteriori. Les codes de dynamique (REBOUND, SWIFT/HJS, ASTROIDE) partent de conditions initiales qu'on leur donne. Entre les deux, chacun écrit ses propres scripts. ORACLE peut couvrir ce lien, avec une interface conçue pour ça.

Code proche à étudier avant de se lancer : l'Exo-Striker (Trifonov). C'est une interface graphique avec un ajustement képlérien ou N-corps et des tests de stabilité, mais surtout pour les vitesses radiales et les transits. ORACLE se distingue par l'imagerie, l'astrométrie, Gaia DR4, les orbites peu couvertes ou hyperboliques et les systèmes hiérarchiques.

## Analyse des solutions dans l'interface (idées d'Antoine)

- **Isoler une famille de solutions** et voir son effet sur les autres distributions. Antoine a commencé cet outil. C'est très utile pour les distributions multimodales, par exemple deux familles en `(Ω, ω)` ou en `e`.
- **Appliquer des formules aux solutions** pour obtenir des grandeurs dérivées, comme la position probable du bord d'un disque autour d'une orbite. Antoine a commencé un fichier séparé pour GG Tau.

Remarques de Claude :
- Rendre les formules générales : une expression écrite par l'utilisateur, évaluée sur tout l'échantillon, qui donne une nouvelle colonne traitée comme les autres (histogramme, corner plot, sélection). On fournit une petite bibliothèque prête : périastre, apoastre, rayon de Hill, bords de disque interne et externe (Holman et Wiegert 1999 pour les binaires), rapports de périodes, inclinaison mutuelle.
- Sélection et formules doivent se combiner : on sélectionne une famille, puis on regarde la grandeur dérivée sur cette famille seulement.
- Garder les poids. Si un jour on repondère les solutions (prior ajouté après coup, critère de stabilité), chaque solution porte un poids et tous les graphiques en tiennent compte.
- Exporter une sélection vers une simulation dynamique (voir plus bas) : c'est le pont naturel entre analyse et simulation.
- F2 (choisir le MAP dans l'analyse) entre dans la même logique.

## Contraindre l'ajustement par la dynamique : trois niveaux

1. **Diagnostics après l'ajustement.** On intègre un échantillon de solutions avec REBOUND ou ASTROIDE : fraction stable, MEGNO, cartes de stabilité, résonances, effet sur un disque. Peu coûteux et parallélisable, c'est le niveau le plus utile tout de suite.
2. **Prior dynamique dans le MCMC.** Des critères analytiques, calculés en quelques microsecondes, ajoutés dans `RATIO_PRIORS` et `COMPUTE_MAP` : stabilité de Hill, stabilité AMD (Laskar et Petit 2017), écart en rayons de Hill mutuels, croisement d'orbites. Pour les systèmes serrés, des classifieurs entraînés comme SPOCK (Tamayo et al. 2020) donnent une probabilité de stabilité. À utiliser comme option, jamais par défaut, et à comparer avec et sans.
3. **Modèle N-corps dans la vraisemblance** (HJS à la place du képlérien). Ce n'est nécessaire que si les interactions se voient sur la durée des données. C'est rare en imagerie, possible pour des systèmes compacts avec vitesses radiales ou Gaia. C'est cher, la convergence est difficile, et c'est le plus gros changement dans le code d'Hervé. Avant de l'envisager, le tester sur des données synthétiques N-corps pour voir quand le képlérien se trompe.

Les critères liés aux disques se généralisent mal (Antoine) : physique du disque, définition du bord, gaz ou poussière. Ils restent des mesures optionnelles, avec une tolérance qui inclut l'incertitude du modèle.

## Connecter l'interface à d'autres codes

- Fitters (orbitize!, Octofitter, orvara) en entrée, codes de dynamique (REBOUND, ASTROIDE) en sortie.
- La vraie difficulté, ce sont les conventions : le ω de la planète ou celui de l'étoile, le sens de Ω, les coordonnées de Jacobi ou astrocentriques, l'époque de référence, les unités. Il faut un format interne unique (solutions, poids, masses, système de coordonnées, époque) et une couche de conversion vérifiée par des tests aller-retour. C'est le même travail que celui fait sur les paramétrisations d'ORACLE.
- Commencer petit : ORACLE et REBOUND, puis un seul autre fitter (orbitize!, le plus utilisé). Chaque code ajouté demande de la maintenance.

## Risque et ordre proposé

Le risque, c'est l'éparpillement. Ordre proposé :
1. corrections ;
2. F1 ;
3. Gaia DR4 ;
4. analyse des solutions (sélection de familles, formules, poids) et diagnostics du niveau 1 avec REBOUND ou ASTROIDE ;
5. priors du niveau 2 ;
6. autres codes ;
7. niveau 3 seulement pour des cas précis.

Avec Gaia DR4, les masses individuelles seront mieux contraintes, ce qui rend les diagnostics dynamiques bien plus parlants.
