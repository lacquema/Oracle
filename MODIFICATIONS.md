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

Ces anomalies existent dans le code de `bcced4b`. `make test` les affiche comme `CONNU`. Une fois corrigée, une anomalie passe dans une entrée `M...` ci-dessous. Corrigées à ce jour : B5 et B6 (M005), B1, B2 et B4 (M006), B8 (M007).

| N° | Où | Problème | Effet constaté |
|---|---|---|---|
| B1 | `kepellip.f:1617`, `kepu.f:2161` | `IF (NPLA.GT.1)` au lieu de `IF (MULTIPLA)` | **Corrigée par M006.** Une orbite avec vitesses radiales : la masse de l'étoile est remplacée par la masse totale, et celle de la planète garde sa valeur d'entrée. Le MCMC démarre avec une planète de masse quasi nulle. Sur un cas simulé, χ² de départ 6227 au lieu de 38. |
| B2 | `kepellip.f:1553`, `kepu.f:2117` | `COV(NPLA,NPLA)` au lieu de `COV(NPAR,NPAR)` | **Corrigée par M006.** Mauvaise erreur affichée sur la masse de l'étoile quand elle est ajustée. |
| B3 | `kepu.f:2179-2181` | `DMUDP(DKAL+1-NEL)` et `DMUDP(DKAL+6-NEL)` pour la première orbite, soit des indices négatifs | Écriture hors du tableau dans `MRQELEMENTS`, à chaque exécution de la version universelle. Effet réel non déterminé. |
| B4 | `kepellip.f:1556`, `kepu.f:2120`, `STORE_SOLUTION` (`io.f:1406-1407`) | Lecture de `P(NEL*NPLA+1)` et `P(NPAR)` sans tester `ISDATA(2)` et `MULTIPLA` | **Corrigée par M006.** Types `1 0 1 0`, une orbite : lecture hors du tableau. Dans `results.dat`, `V0` contient en fait χ²/`MPS` et `m0` vaut environ 3379 : colonnes sans signification. Vérifié sur un run réel et sur un cas simulé. Les colonnes orbitales ne sont pas touchées. |
| B5 | `kepu.f:592` | `GM` n'est rempli que si `MULTIPLA` | **Corrigée par M005.** Une orbite sans vitesses de l'étoile : `CHI2SEUL`, donc le MCMC universel, utilisait une masse non initialisée. Sur un cas simulé, χ² de 18659 au lieu de 279. |
| B6 | `kepu.f:611-612` | `sin i` manque dans `AMX` et `AMY`, et `AMX` utilise `SIGMA` au lieu de `GM` | **Corrigée par M005.** `CHI2SEUL` universel : mauvaise amplitude des vitesses radiales. Sur un cas simulé, χ² de 735 au lieu de 377. |
| B7 | `MRQVFITP` de `kepu.f`, cause non analysée | Dérivées de vitesse de la planète incorrectes pour les paramètres 1 et 6 | Gradient du LM faux de 4 % sur `p1` dans le cas `rel_rvp`. |
| B8 | `io.f`, `ENTER_PRIOR` | `READ(5,'(a)')` pour des nombres, priors de types 2, 3 et 4 | **Corrigée par M007 (lecture seulement, voir B13).** Le programme s'arrête à la lecture des priors log-normal, linéaire et fixe, écrits comme le fait l'interface. Seul le type 1 fonctionne. |
| B9 | `COMPUTE_MAP`, `RATIO_PRIORS` (`kepellip.f:188`, `kepellip.f:556`), même schéma dans `kepu.f` | Une orbite sans vitesses de l'étoile (`MULTIPLA` faux), prior sur `m0` (coefficients `1 0`) : le terme `BCOF(0)*MSTAR` est ignoré et le prior s'évalue sur une masse nulle, sans avertissement | Le prior n'agit pas, avec ou sans vitesses de la planète. Confirmé sur des runs réels lancés avec et sans ce prior : mêmes résultats, `Map` décalée d'une constante. Exemple : 51 Eri (une orbite, astrométrie relative et vitesses de la planète, prior 1,55 ± 0,03 Msun), où `a` et `P` ajustés donnent `a³/P²` ≈ 1,36 Msun, loin de la valeur du prior. Les coefficients `1 1` donnent bien un prior sur la masse totale (vérifié par `make test`). |
| B10 | `RATIO_PRIORS` de `kepu.f` (facteur `(Q/QT)*(UPQT/UPQ)`) | Le rapport contient `(q/q')(P/P')` en plus du prior `1/(q·P)` de `COMPUTE_MAP` et du commentaire. Le jacobien du commentaire oublie que les paramètres sont `ln q` et `ln P_q` | Le MCMC universel échantillonne en `1/(q²·P²)` au lieu de `1/(q·P)`. Mesuré exactement par le jacobien numérique. Intention à confirmer avec Hervé. |
| B11 | `RATIO_PRIORS` des deux versions, boucle sur les priors de masse | Chaque prior ajoute `M'/M` et `m0'/m0`, y compris les priors de masse positive et plusieurs fois pour la même masse. `COMPUTE_MAP` ne les contient pas | Quand les masses sont ajustées, le prior réellement échantillonné diffère de celui qu'indique `COMPUTE_MAP`, et favorise les grandes masses. Sur un run réel à deux orbites, `m0` 0,563 Msun contre 0,536 avec un prior plat sur les masses, estimé par repondération ; orbites externes presque inchangées. Prior voulu à confirmer avec Hervé. |
| B12 | `GELMAN_RUBIN` (`mcmc.f:347`), départ des chaînes | Les 8 chaînes partent du même point (solution du LM). Au premier calcul, l'écart entre chaînes est quasi nul, donc `r` ≈ 1 et `ts` est énorme : le critère d'inférence (`RINF`, `TSINF`) est atteint avant que les chaînes aient exploré | Constaté le 2026-11-06 sur `rel_vr_planete` en classique (relative + vitesses de la planète, `--rapide`) : inférence dès le premier contrôle dans 5 runs sur 5, avec le code de `a83b9d1` comme après M006. Distributions trop étroites ou asymétriques, `P` entre 1,1 et 3,1 sigma de la vérité selon le run. Même mécanisme probable pour la fausse convergence de `arc_periastre` en universel (`Tests/parametrisation/README.md`). Correction à discuter avec Hervé : départs dispersés, période de chauffe, ou critère sur la seconde moitié des chaînes. |
| B13 | `USER_PRIOR_MAP` et `USER_PRIOR_RATIO` (`mcmc.f:710-812`), **à vérifier** | Repéré à la lecture le 2026-11-06, non testé : type 3 (linéaire) : rapport toujours égal à 1 et aucun terme dans le MAP, les bornes lues ne semblent servir nulle part ; type 4 (fixe) : gaussienne de largeur relative 1e-3 dans le MAP mais 1e-6 dans le rapport ; type 6 (Heaviside, prior par défaut `MPRIOR(0)`, `io.f:385`) : dans `USER_PRIOR_MAP`, `EX = EXP(X)` avec `X` non initialisé, et `EX` est un argument, donc la masse de l'appelant serait écrasée | Effet non mesuré. Le cas 6 toucherait le MAP de tous les runs ; pourtant `make test` trouve `MAP = χ² + Fmap` cohérent : à comprendre avant de conclure. |

**Fragilité constatée :** un job tué pendant l'écriture de `dump.dat` laisse un fichier illisible à la reprise (`End of file`), car le fichier est écrit sur place.

**À savoir, sans être une anomalie du calcul :** avec une orbite et sans vitesses de l'étoile, les masses individuelles ne sont pas ajustées. La colonne `m` contient alors la masse totale, dans l'unité de la planète, et `m0` la valeur d'entrée. Une valeur comme 1079 pour une planète en Mjup vient de là (1,03 Msun).

**Déjà corrigé avant ce journal :** une colonne `m` infinie apparaît dans les runs continués avec un ancien binaire : `PLA%MUNIT` n'était pas gardé dans le fichier de sauvegarde. Le commit `83ee116` du 2025-12-12 l'a corrigé.

Correction de B5 et B6 essayée sur une copie : `kepu.f` retrouve alors le χ² du modèle indépendant dans tous les cas de test. La correction de B1 a été essayée de la même façon, sur le cas simulé relatif + vitesses radiales.

## Modifications

## M007 — Corrige B8 : lecture des priors de masse log-normal, linéaire et fixe

- **Date :** 2026-11-06
- **Fichier :** `Code/mcmc/sub/io.f`, routine `ENTER_PRIOR` (trois blocs repérés `M007`).
- **Changement :** `READ(5,'(a)')` devient `READ(5,*)` pour les types 2, 3 et 4, comme pour le type 1 juste au-dessus. Le format `'(a)'` lit du texte : le programme s'arrêtait en lisant des nombres.
- **Vérification :** `make test` : les trois priors, écrits comme le fait l'interface, se lisent sans erreur. **Seule la lecture est vérifiée** : l'effet de ces priors dans le MCMC n'est pas testé, et la lecture du code a révélé des points suspects (B13).
- **Effet sur les résultats :** les priors de types 2, 3 et 4 deviennent utilisables. Aucun run existant ne les utilisait, puisque le programme s'arrêtait.

## M006 — Corrige B1, B2 et B4 : masse de l'étoile, son erreur, et paramètres lus sans exister

- **Date :** 2026-11-06
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `Code/mcmc/sub/kepu.f` (routines `MRQELEMENTS` et `MRQRENEW`, mêmes changements dans les deux), `Code/mcmc/sub/io.f` (routine `STORE_SOLUTION`). Blocs repérés `M006`.
- **B1 (`MRQELEMENTS`) :** `IF (NPLA.GT.1)` devient `IF (MULTIPLA)`. La masse de l'étoile est un paramètre dès que `MULTIPLA` est vrai, donc aussi avec une seule orbite et des vitesses de l'étoile. Avant, ce cas remplaçait la masse de l'étoile par la masse totale après le LM, et la planète gardait sa masse d'entrée : le MCMC partait d'une planète de masse presque nulle.
- **B2 (`MRQELEMENTS`) :** l'erreur sur la masse de l'étoile lisait `COV(NPLA,NPLA)` au lieu de `COV(NPAR,NPAR)`, la case de `ln m0`.
- **B4 :** `V0` et le jitter n'existent comme paramètres qu'avec les vitesses de l'étoile (`ISDATA(2)`), et `m0` seulement si `MULTIPLA`. Avant, le code les lisait ou les écrivait dès qu'il y avait des vitesses radiales (`RADVEL`), y compris avec les seules vitesses de la planète, hors du tableau des paramètres :
  - `MRQELEMENTS` : `IF (RADVEL)` devient `IF (ISDATA(2))` pour `V0` et le jitter ;
  - `MRQRENEW` (nouveau tirage des paramètres quand le LM tombe sur une orbite impossible) : même changement, l'écriture de `V0` débordait ;
  - `STORE_SOLUTION` : `m0` lue seulement si `MULTIPLA`, `V0` et jitter seulement si `ISDATA(2)`. Sinon, ces colonnes de `results.dat` gardent une valeur fixe, celle de `STAR%MASS`, `STAR%V0` et `STAR%SIGJV` au moment de l'écriture : elles ne sont pas ajustées.
- **Vérification :**
  - `make test` : 349 OK, 0 échec, référence inchangée. Les tests de masse marqués B1 et B2 passent maintenant (une orbite avec vitesses de l'étoile, et deux orbites) ; le contrôle des tableaux (`-fcheck=all`) ne s'arrête plus sur `P` dans `rel_rvp` (B4) mais plus loin, sur B3, toujours présent ;
  - `rel_vr_etoile` (`Tests/synthetic.py`, `--rapide`) : le MCMC part de χ² = 86,4 au lieu de 7804, dans les deux versions, et retrouve tous les paramètres à moins de 1,7 sigma (m = 30,4 ± 0,4 Mjup pour 30) ;
  - `rel_vr_planete` en classique (types `1 0 1 0`) : dans `results.dat`, `V0` et `Jit` valent 0 et `m0` est constante, au lieu de χ²/`MPS` et d'environ 3379. Ce run a révélé B12, présent aussi avant M006 (voir la table).
- **Effet sur les résultats :** aucun χ² ni MAP ne change (la référence des tests est inchangée). Changent : le point de départ du MCMC après le LM pour une orbite avec vitesses de l'étoile (B1), l'erreur affichée sur `m0` après le LM (B2), et les colonnes `V0`, `Jit` et `m0` de `results.dat` quand ces paramètres n'existent pas (B4).

## M005 — Corrige B5 et B6 : masse et amplitudes des vitesses radiales dans `CHI2SEUL` universel

- **Date :** 2026-11-06
- **Fichier :** `Code/mcmc/sub/kepu.f`, routine `CHI2SEUL` (deux blocs repérés `M005`).
- **B5 :** avec une seule orbite et sans vitesses de l'étoile (`MULTIPLA` faux), le tableau `GM` des masses dynamiques n'était jamais rempli. Ajout de `ELSE` / `GM = MUF`, où `MUF` contient alors la masse totale de l'orbite (voir `ELEMENTS`).
- **B6 :** les amplitudes `AMX`, `AMY` des vitesses radiales oubliaient `sin i`, et `AMX` utilisait `SIGMA` (un scalaire, la masse cumulée de la dernière orbite) au lieu de `GM`. Formules remplacées par `AMX = -GM*SI*SOM` et `AMY = SQRT(QQ*GM*(1+EXC))*SI*COM`.
- **Pourquoi :** le MCMC universel calculait un χ² faux dès qu'il y avait des vitesses radiales, ou avec une orbite sans vitesses de l'étoile. `kepellip.f` n'est pas concerné.
- **Vérification :**
  - le code est identique, aux repères près, à la copie corrigée testée le 2026-11-05 (`Tests/parametrisation/patch_B5B6.py`) ;
  - `make test` : 324 OK, 0 échec. Dans `Tests/reference/reference.json`, refaite pour cette correction, seules les 10 entrées universelles touchées changent (`rel`, `rel_corr`, `rel_seppa`, `rel_rv`, `rel_rv_jit`, `rel_rv_prior`, `rel_rv_rvp`, `rel_rv_2orb_jit`, `rel_rvp`, `rv_seul`). Les 12 autres sont identiques, dont tous les cas classiques et `rel_2orb` universel. Le χ² universel est maintenant égal à celui du modèle Python et de la version classique dans tous les cas ; `KNOWN_CHI2_UNIV` est vidé dans `Tests/run_tests.py` ;
  - systèmes synthétiques en universel (`Tests/synthetic.py`, `--rapide`) : `rel_vr_etoile` donne m = 30,4 ± 0,4 Mjup pour 30 vrais (19,5 avant) et retrouve tous les paramètres à moins de 1,6 sigma ; `rel_vr_planete`, auparavant bloqué à χ² ≈ 2e11, va au bout et retrouve tout à moins de 1,9 sigma.
- **Effet sur les résultats :** les runs universels avec vitesses radiales, ou à une orbite sans vitesses de l'étoile, changent. Les runs classiques et les runs universels à astrométrie relative seule avec plusieurs orbites ne changent pas.

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
