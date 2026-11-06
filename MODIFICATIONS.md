# Journal des modifications d'ORACLE

Toute modification du code depuis le commit `bcced4b` est listée ici. Le code de calcul de référence est celui de ce commit.

## Comment lire ce journal

- Chaque modification a un numéro `M001`, `M002`, etc., du plus ancien au plus récent.
- Le même numéro entoure les lignes concernées dans le code. Pour tout retrouver : `grep -rn "\[MOD M" Code Interface Makefile`.

```fortran
c>>> [MOD M001] Résumé de la modification
c    AVANT : ancienne ligne, recopiée telle quelle
        nouvelle ligne
c<<< [MOD M001]
```

En Python et dans le `Makefile`, le repère s'écrit `# >>> [MOD M001]` et `# <<< [MOD M001]`.

`make test` vérifie que toute ligne changée depuis `bcced4b` est bien dans un bloc repéré et que tout repère figure ici. Une entrée sans repère porte la mention `(sans repère)` : c'est réservé à ce qui ne se commente pas, comme la suppression de fichiers.

## Anomalies connues, non corrigées

Ces anomalies existent dans le code de `bcced4b`. `make test` les affiche comme `CONNU`. Une fois corrigée, une anomalie passe dans une entrée `M...` ci-dessous. Aucune n'est corrigée à ce jour.

| N° | Où | Problème | Effet constaté |
|---|---|---|---|
| B1 | `kepellip.f:1617`, `kepu.f:2161` | `IF (NPLA.GT.1)` au lieu de `IF (MULTIPLA)` | Une orbite avec vitesses radiales : la masse de l'étoile est remplacée par la masse totale, et celle de la planète garde sa valeur d'entrée. Le MCMC démarre avec une planète de masse quasi nulle. Sur un cas simulé, χ² de départ 6227 au lieu de 38. |
| B2 | `kepellip.f:1553`, `kepu.f:2117` | `COV(NPLA,NPLA)` au lieu de `COV(NPAR,NPAR)` | Mauvaise erreur affichée sur la masse de l'étoile quand elle est ajustée. |
| B3 | `kepu.f:2179-2181` | `DMUDP(DKAL+1-NEL)` et `DMUDP(DKAL+6-NEL)` pour la première orbite, soit des indices négatifs | Écriture hors du tableau dans `MRQELEMENTS`, à chaque exécution de la version universelle. Effet réel non déterminé. |
| B4 | `kepellip.f:1556`, `kepu.f:2120`, `STORE_SOLUTION` (`io.f:1406-1407`) | Lecture de `P(NEL*NPLA+1)` et `P(NPAR)` sans tester `ISDATA(2)` et `MULTIPLA` | Types `1 0 1 0`, une orbite : lecture hors du tableau. Dans `results.dat`, `V0` contient en fait χ²/`MPS` et `m0` vaut environ 3379 : colonnes sans signification. Vérifié sur un run réel et sur un cas simulé. Les colonnes orbitales ne sont pas touchées. |
| B5 | `kepu.f:592` | `GM` n'est rempli que si `MULTIPLA` | Une orbite sans vitesses de l'étoile : `CHI2SEUL`, donc le MCMC universel, utilise une masse non initialisée. Sur un cas simulé, χ² de 18659 au lieu de 279. |
| B6 | `kepu.f:611-612` | `sin i` manque dans `AMX` et `AMY`, et `AMX` utilise `SIGMA` au lieu de `GM` | `CHI2SEUL` universel : mauvaise amplitude des vitesses radiales. Sur un cas simulé, χ² de 735 au lieu de 377. |
| B7 | `MRQVFITP` de `kepu.f`, cause non analysée | Dérivées de vitesse de la planète incorrectes pour les paramètres 1 et 6 | Gradient du LM faux de 4 % sur `p1` dans le cas `rel_rvp`. |
| B8 | `io.f`, `ENTER_PRIOR` | `READ(5,'(a)')` pour des nombres, priors de types 2, 3 et 4 | Le programme s'arrête à la lecture des priors log-normal, linéaire et fixe, écrits comme le fait l'interface. Seul le type 1 fonctionne. |

| B9 | `COMPUTE_MAP`, `RATIO_PRIORS` (`kepellip.f:188`, `kepellip.f:556`), même schéma dans `kepu.f` | Une orbite sans vitesses de l'étoile (`MULTIPLA` faux), prior sur `m0` (coefficients `1 0`) : le terme `BCOF(0)*MSTAR` est ignoré et le prior s'évalue sur une masse nulle, sans avertissement | Le prior n'agit pas, avec ou sans vitesses de la planète. Confirmé sur des runs réels lancés avec et sans ce prior : mêmes résultats, `Map` décalée d'une constante. Exemple : 51 Eri (une orbite, astrométrie relative et vitesses de la planète, prior 1,55 ± 0,03 Msun), où `a` et `P` ajustés donnent `a³/P²` ≈ 1,36 Msun, loin de la valeur du prior. Les coefficients `1 1` donnent bien un prior sur la masse totale (vérifié par `make test`). |

| B10 | `RATIO_PRIORS` de `kepu.f` (facteur `(Q/QT)*(UPQT/UPQ)`) | Le rapport contient `(q/q')(P/P')` en plus du prior `1/(q·P)` de `COMPUTE_MAP` et du commentaire. Le jacobien du commentaire oublie que les paramètres sont `ln q` et `ln P_q` | Le MCMC universel échantillonne en `1/(q²·P²)` au lieu de `1/(q·P)`. Mesuré exactement par le jacobien numérique. Intention à confirmer avec Hervé. |
| B11 | `RATIO_PRIORS` des deux versions, boucle sur les priors de masse | Chaque prior ajoute `M'/M` et `m0'/m0`, y compris les priors de masse positive et plusieurs fois pour la même masse. `COMPUTE_MAP` ne les contient pas | Quand les masses sont ajustées, le prior réellement échantillonné diffère de celui qu'indique `COMPUTE_MAP`, et favorise les grandes masses. Sur un run réel à deux orbites, `m0` 0,563 Msun contre 0,536 avec un prior plat sur les masses, estimé par repondération ; orbites externes presque inchangées. Prior voulu à confirmer avec Hervé. |

**Fragilité constatée :** un job tué pendant l'écriture de `dump.dat` laisse un fichier illisible à la reprise (`End of file`), car le fichier est écrit sur place.

**À savoir, sans être une anomalie du calcul :** avec une orbite et sans vitesses de l'étoile, les masses individuelles ne sont pas ajustées. La colonne `m` contient alors la masse totale, dans l'unité de la planète, et `m0` la valeur d'entrée. Une valeur comme 1079 pour une planète en Mjup vient de là (1,03 Msun).

**Déjà corrigé avant ce journal :** une colonne `m` infinie apparaît dans les runs continués avec un ancien binaire : `PLA%MUNIT` n'était pas gardé dans le fichier de sauvegarde. Le commit `83ee116` du 2025-12-12 l'a corrigé.

Correction de B5 et B6 essayée sur une copie : `kepu.f` retrouve alors le χ² du modèle indépendant dans tous les cas de test. La correction de B1 a été essayée de la même façon, sur le cas simulé relatif + vitesses radiales.

## Modifications

## M004 — Tests de paramétrisation et lecture Gaia (sans repère)

- **Date :** 2026-11-06
- **Fichiers :** `Tests/parametrisation/` et `Tests/gaia/` ajoutés ; `Tests/synthetic.py` complété (systèmes `arc_tres_court`, `arc_periastre`, `arc_parabole`, bornes propres à chaque système, analyse en `q` pour `e` proche de 1, variable `ORACLE_MAIN`) ; `Tests/README.md` mis à jour.
- **Pourquoi :** garder la comparaison des paramétrisations classique, universelle et lambda, qui a conduit à choisir lambda (voir `Tests/parametrisation/README.md`), et vérifier la lecture des données Gaia DR4 sur BH3 avant d'écrire le code Gaia.
- **Vérification :** `prepare.py` refait des copies identiques à celles des runs du 2026-11-05 et 06 ; la variante lambda passe sa vérification sur les 11 cas ; `bh3_esa.py` retrouve le notebook de l'ESA au dernier chiffre affiché près.
- **Effet sur les résultats :** aucun, les routines de calcul du dépôt n'ont pas changé. Les corrections B5, B6, B10 et la variante lambda restent dans des copies hors du dépôt.

## M003 — Tests de vérification (sans repère)

- **Date :** 2026-11-05
- **Fichiers :** dossier `Tests/` ajouté (voir `Tests/README.md`), ce journal ajouté.
- **Pourquoi :** pouvoir prouver qu'une modification ne change rien d'autre que ce qu'elle doit changer.
- **Vérification :** trois erreurs volontaires, introduites sur une copie ou retirées aussitôt, sont détectées : signe faux dans `POSFITS`, dérivée faussée de 0,1 %, fichier ajouté sans repère.
- **Systèmes synthétiques :** `Tests/synthetic.py` fait tourner le vrai programme sur des systèmes inventés et compare à la vérité. Ce test est long et hors de `make test`.
- **Référence :** `Tests/reference/reference.json` enregistre le comportement du code de `bcced4b`, anomalies comprises.
- **Effet sur les résultats :** aucun, les routines de calcul n'ont pas changé.

## M002 — Cible `make test`

- **Date :** 2026-11-05
- **Fichiers :** `Makefile`, et `Tests/build/` ajouté au `.gitignore`.
- **Pourquoi :** lancer les tests de vérification en une commande.
- **Effet sur les résultats :** aucun.

## M001 — Correction de la cible `mod` du `Makefile`

- **Date :** 2026-11-05
- **Fichier :** `Makefile`
- **Pourquoi :** la cible testait `main/data.mod` mais supprimait `./data.mod`, ce qui aurait fait échouer `make` sur un clone neuf une fois `data.mod` retiré du suivi.
- **Vérification :** `make compile` réussi depuis une copie propre du dépôt.
- **Effet sur les résultats :** aucun.

## M000 — Rangement du dépôt (sans repère)

- **Date :** 2026-11-05
- **Commit :** `1128db9`
- **Fichiers :** `.gitignore` ajouté ; `Code/bin/.gitkeep` et `Code/lib/.gitkeep` ajoutés. Retirés du suivi, conservés sur le disque : binaires, bibliothèques, `data.mod`, `.pyc`, `.DS_Store`, copies de sauvegarde (`~io.f`, `WidgetPlot~.py`) et `.last_path`.
- **Pourquoi :** éviter que chaque compilation modifie des fichiers suivis.
- **Conséquence :** après un changement de branche ou un clone, lancer `make compile`.
- **Effet sur les résultats :** aucun, le code de calcul n'a pas changé.
