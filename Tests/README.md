# Tests de vérification

```
make test
```

La commande recompile les routines actuelles dans `Tests/build/`, sans toucher à `Code/bin` ni `Code/lib`. Elle demande `gfortran` et `numpy`. Ajouter `PYTHON3=...` si besoin.

Chaque test affiché vaut `ECHEC` (à corriger) ou `CONNU` (anomalie listée dans `MODIFICATIONS.md`). La dernière ligne donne le bilan. Le code de sortie vaut 0 s'il n'y a aucun `ECHEC`. `python3 Tests/run_tests.py -v` affiche aussi les tests réussis.

## Ce qui est vérifié

| Test | Question posée |
|---|---|
| Journal et repères | Toute ligne changée depuis `bcced4b` est-elle dans un bloc repéré et listée dans le journal ? |
| Exécution et contrôles | Les routines tournent-elles sans erreur, y compris avec les contrôles de bornes de tableaux (`-fcheck=all`) ? |
| Référence | χ², MAP et χ²+Fmap sont-ils identiques à ceux du code de départ (à 1e-10 près) ? |
| Modèle indépendant | Un modèle Python écrit séparément (`model.py`) donne-t-il le même χ² que le Fortran ? |
| Classique et universel | `kepellip.f` et `kepu.f` donnent-ils le même χ² pour la même orbite ? |
| Cohérence interne | `CHI2SEUL` et `MRQCOF` donnent-ils le même χ² ? `COMPUTE_MAP` est-il cohérent avec `MRQCOF` ? |
| Rapport des priors | `RATIO_PRIORS` correspond-il à `COMPUTE_MAP` et au jacobien calculé numériquement par `jacobian.py` ? |
| Dérivées | Le gradient du Levenberg-Marquardt correspond-il aux différences finies ? |
| Éléments | `MRQELEMENTS` redonne-t-il la masse de l'étoile et son erreur ? |
| Priors | Les priors sur les masses se lisent-ils, écrits comme le fait l'interface ? |
| Types de priors | Pour chaque type de prior de masse (0, 1, 2, 3, 4 et 6), le MAP (`USER_PRIOR_MAP`), le MCMC (`USER_PRIOR_RATIO`) et le L-M (`FMAP_PRIOR`) donnent-ils la densité voulue, et les masses hors du support (bornes, masse négative) sont-elles refusées ? |
| Essai réel | Le vrai programme retrouve-t-il la solution de GG Tau ? |

Les tests s'évaluent à paramètres fixés, sans hasard, sauf l'essai réel qui compare à une tolérance large.

## Les fichiers

- `cases/` : données simulées et entrées, à ne pas régénérer sans raison (`generate_cases.py`, graine fixe).
- `reference/reference.json` : valeurs du code de départ.
- `driver.f` : programme qui appelle les vraies routines et imprime leurs résultats.
- `priors.f` : programme qui appelle directement les trois routines de prior de masse, pour un prior et une liste de masses.
- `model.py` : modèle indépendant, `cases_definition.py` : définition des cas.
- `synthetic.py` : systèmes synthétiques lancés avec le vrai programme jusqu'au bout (long, hors de `make test`).
- `parametrisation/` : comparaison classique, universelle et variante lambda, avec les scripts pour refaire les copies corrigées et les résultats. Voir son `README.md`.
- `gaia/bh3_esa.py` : refait l'ajustement de Gaia BH3 du notebook de l'ESA à partir des données brutes de la pré-publication DR4 (fichier hors du dépôt), pour vérifier la lecture des données avant le code Gaia.

## Quand refaire la référence

Seulement si un changement doit modifier les résultats. Dans ce cas : `python3 Tests/run_tests.py --update-reference`, puis expliquer pourquoi dans `MODIFICATIONS.md`. Sinon, un écart à la référence signale une erreur.

## Essayer une correction hors du dépôt

Copier `Code/mcmc/sub` ailleurs, y faire la modification, puis :

```
ORACLE_SUB=/chemin/copie ORACLE_BUILD=/chemin/build python3 Tests/run_tests.py
```

## Ce qui n'est pas encore vérifié

- La couverture statistique du MCMC sur de nombreux tirages. Le retour à la solution est testé à part par `Tests/synthetic.py`, qui est long et n'est pas dans `make test`.
- Les colonnes de `results.dat` : vérifiées une fois à la main le 2026-11-05 (types `1 0 0 0`, `1 1 0 0`, `1 0 1 0`, `1 1 1 0`, `0 1 0 0`, deux orbites), pas encore dans `make test`.
- Le code Gaia, qui n'existe pas encore : ses tests viendront avec lui.
