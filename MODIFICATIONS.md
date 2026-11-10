# Journal des modifications d'ORACLE

Toute modification du code depuis le commit `bcced4b` est listée ici. Le code de calcul de référence est celui de ce commit.

## Comment lire ce journal

- Chaque modification a un numéro `M001`, `M002`, etc., du plus ancien au plus récent.
- Le même numéro entoure les lignes concernées dans le code d'Hervé (`Code/mcmc` et `Makefile`). Pour tout retrouver : `grep -rn "\[MOD M" Code Makefile`. L'interface (`Interface/`) est le code d'Antoine : ses changements sont décrits ici, sans repère dans le code.

```fortran
c>>> [MOD M001] Résumé de la modification
c    AVANT : ancienne ligne, recopiée telle quelle
        nouvelle ligne
c<<< [MOD M001]
```

En Python et dans le `Makefile`, le repère s'écrit `# >>> [MOD M001]` et `# <<< [MOD M001]`.

`make test` vérifie que toute ligne changée depuis `bcced4b` est bien dans un bloc repéré et que tout repère figure ici. Une entrée sans repère porte la mention `(sans repère)` : c'est réservé à ce qui ne se commente pas, comme la suppression de fichiers.

## Anomalies connues, non corrigées

Ces anomalies existent dans le code de `bcced4b`. `make test` les affiche comme `CONNU`. Une fois corrigée, une anomalie passe dans une entrée `M...` ci-dessous. Corrigées à ce jour : B5 et B6 (M005), B1, B2 et B4 (M006), B8 (M007), B3 (M008), B7 (M009), B13 (M010), B10 (M011), B9 (M012), B11 (M013), B12 (M015), B15 (M016), B18 (M020), B19 (M022).

| N° | Où | Problème | Effet constaté |
|---|---|---|---|
| B1 | `kepellip.f:1617`, `kepu.f:2161` | `IF (NPLA.GT.1)` au lieu de `IF (MULTIPLA)` | **Corrigée par M006.** Une orbite avec vitesses radiales : la masse de l'étoile est remplacée par la masse totale, et celle de la planète garde sa valeur d'entrée. Le MCMC démarre avec une planète de masse quasi nulle. Sur un cas simulé, χ² de départ 6227 au lieu de 38. |
| B2 | `kepellip.f:1553`, `kepu.f:2117` | `COV(NPLA,NPLA)` au lieu de `COV(NPAR,NPAR)` | **Corrigée par M006.** Mauvaise erreur affichée sur la masse de l'étoile quand elle est ajustée. |
| B3 | `kepu.f:2179-2181` | `DMUDP(DKAL+1-NEL)` et `DMUDP(DKAL+6-NEL)` pour la première orbite, soit des indices négatifs | **Corrigée par M008.** Précision : les indices négatifs viennent de la remise à zéro faite après le bloc `IF (I.EQ.1)`, pas du bloc `ELSE`. Écriture hors du tableau dans `MRQELEMENTS`, à chaque exécution de la version universelle. Effet réel non déterminé. |
| B4 | `kepellip.f:1556`, `kepu.f:2120`, `STORE_SOLUTION` (`io.f:1406-1407`) | Lecture de `P(NEL*NPLA+1)` et `P(NPAR)` sans tester `ISDATA(2)` et `MULTIPLA` | **Corrigée par M006.** Types `1 0 1 0`, une orbite : lecture hors du tableau. Dans `results.dat`, `V0` contient en fait χ²/`MPS` et `m0` vaut environ 3379 : colonnes sans signification. Vérifié sur un run réel et sur un cas simulé. Les colonnes orbitales ne sont pas touchées. |
| B5 | `kepu.f:592` | `GM` n'est rempli que si `MULTIPLA` | **Corrigée par M005.** Une orbite sans vitesses de l'étoile : `CHI2SEUL`, donc le MCMC universel, utilisait une masse non initialisée. Sur un cas simulé, χ² de 18659 au lieu de 279. |
| B6 | `kepu.f:611-612` | `sin i` manque dans `AMX` et `AMY`, et `AMX` utilise `SIGMA` au lieu de `GM` | **Corrigée par M005.** `CHI2SEUL` universel : mauvaise amplitude des vitesses radiales. Sur un cas simulé, χ² de 735 au lieu de 377. |
| B7 | `MRQVFITP` de `kepu.f` | Dérivées de vitesse de la planète incorrectes pour les paramètres 1 et 6 | **Corrigée par M009** (division par `MTOT` au lieu de `MDYN(I)`). Gradient du LM faux de 4 % sur `p1` dans le cas `rel_rvp`. |
| B8 | `io.f`, `ENTER_PRIOR` | `READ(5,'(a)')` pour des nombres, priors de types 2, 3 et 4 | **Corrigée par M007 (lecture seulement, voir B13).** Le programme s'arrête à la lecture des priors log-normal, linéaire et fixe, écrits comme le fait l'interface. Seul le type 1 fonctionne. |
| B9 | `COMPUTE_MAP`, `RATIO_PRIORS` (`kepellip.f:188`, `kepellip.f:556`), même schéma dans `kepu.f` | Une orbite sans vitesses de l'étoile (`MULTIPLA` faux), prior sur `m0` (coefficients `1 0`) : le terme `BCOF(0)*MSTAR` est ignoré et le prior s'évalue sur une masse nulle, sans avertissement | **Corrigée par M012** (prior refusé à la lecture, avec un message). Le prior n'agissait pas, avec ou sans vitesses de la planète. Confirmé sur des runs réels lancés avec et sans ce prior : mêmes résultats, `Map` décalée d'une constante. Exemple : 51 Eri (une orbite, astrométrie relative et vitesses de la planète, prior 1,55 ± 0,03 Msun), où `a` et `P` ajustés donnent `a³/P²` ≈ 1,36 Msun, loin de la valeur du prior. Les coefficients `1 1` donnent bien un prior sur la masse totale (vérifié par `make test`). |
| B10 | `RATIO_PRIORS` de `kepu.f` (facteur `(Q/QT)*(UPQT/UPQ)`) | Le rapport contient `(q/q')(P/P')` en plus du prior `1/(q·P)` de `COMPUTE_MAP` et du commentaire. Le jacobien du commentaire oublie que les paramètres sont `ln q` et `ln P_q` | **Corrigée par M011** (le nouveau rapport des priors est calculé pour la paramétrisation lambda, sans ce facteur). Le MCMC universel échantillonnait en `1/(q²·P²)` au lieu de `1/(q·P)`. Mesuré exactement par le jacobien numérique. Intention à confirmer avec Hervé. |
| B11 | `RATIO_PRIORS` des deux versions, boucle sur les priors de masse | Chaque prior ajoute `M'/M` et `m0'/m0`, y compris les priors de masse positive et plusieurs fois pour la même masse. `COMPUTE_MAP` ne les contient pas | **Corrigée par M013** (choix « multiplier »). Quand les masses étaient ajustées, le prior réellement échantillonné diffère de celui qu'indique `COMPUTE_MAP`, et favorise les grandes masses. Sur un run réel à deux orbites, `m0` 0,563 Msun contre 0,536 avec un prior plat sur les masses, estimé par repondération ; orbites externes presque inchangées. Prior voulu à confirmer avec Hervé. |
| B12 | `GELMAN_RUBIN` (`mcmc.f:347`), départ des chaînes | Les 8 chaînes partent du même point (solution du LM). Au premier calcul, l'écart entre chaînes est quasi nul, donc `r` ≈ 1 et `ts` est énorme : le critère d'inférence (`RINF`, `TSINF`) est atteint avant que les chaînes aient exploré | **Corrigée par M015.** Constaté le 2026-11-06 sur `rel_vr_planete` en classique (relative + vitesses de la planète, `--rapide`) : inférence dès le premier contrôle dans 5 runs sur 5, avec le code de `a83b9d1` comme après M006. Distributions trop étroites ou asymétriques, `P` entre 1,1 et 3,1 sigma de la vérité selon le run. Même mécanisme probable pour la fausse convergence de `arc_periastre` en universel (`Tests/parametrisation/README.md`). Correction à discuter avec Hervé : départs dispersés, période de chauffe, ou critère sur la seconde moitié des chaînes. |
| B13 | `USER_PRIOR_MAP` et `USER_PRIOR_RATIO` (`mcmc.f`) | Type 0 : rapport du MCMC `ln m / ln m_essai` au lieu de `m / m_essai`. Type 3 : bornes lues mais jamais appliquées. Type 4 : largeur 1e-3 dans le MAP, 1e-6 dans le MCMC et le L-M. Type 6 : dans le MAP, `EX = EXP(X)` avec `X` non initialisé (en pratique 0, donc toute masse remplacée par 1) | **Corrigée par M010.** Le type 3 (« Uniform » de l'interface) ne bornait pas la masse. Le type 6 (masses positives, ajouté automatiquement) laissait le MAP sans pénalité pour une masse négative ; le MCMC n'était pas touché. Type 0 : jamais écrit par l'interface. |
| B14 | Tous les appels à `RANDOM_NUMBER` (`mcmc.f`, `MRQRENEW` de `kepellip.f` et `kepu.f`, `io.f`), sans `RANDOM_SEED` | La graine du générateur aléatoire n'est jamais fixée. Avec gfortran, elle change alors à chaque exécution | Deux lancements identiques ne donnent pas le même résultat : on ne peut pas reproduire un run à l'identique. Le LM en dépend aussi, car `MRQRENEW` repart d'angles tirés au hasard après une solution impossible : sur GG Tau, il finit dans un autre minimum environ une fois sur cent (396 sur 400 avant M010, 398 sur 400 après). Choix à faire avec Hervé : graine fixe, graine donnée en entrée, ou graine écrite dans la sortie pour pouvoir refaire un run. Dans `make test`, l'essai GG Tau est relancé jusqu'à trois fois pour ne pas échouer au hasard. |
| B15 | `MCMCP` (`mcmc.f`), intervalle entre deux contrôles de Gelman-Rubin une fois la convergence atteinte | `NGEL = MIN(NGROUP,NINT(LONG*1d-2/NGROUP)*NGROUP)` vaut 0 tant que `LONG` < 5 000 balayages : aucun pas entre deux contrôles | **Corrigée par M016.** Quand la convergence arrivait avant 5 000 balayages, les chaînes ne bougeaient plus et les mêmes 8 positions étaient gardées à chaque contrôle jusqu'à remplir la réserve. Constaté le 2026-11-08 sur `rel_vr_planete` (`--rapide`) : 9 solutions distinctes sur 20 001, ω hors de l'intervalle à 99,7 %. Explique aussi une partie de B12. Les longs runs réels (GG Tau, 51 Eri), qui convergent après des millions de pas, ne sont pas touchés. |
| B16 | `MRQELEMENTS` (`kepellip.f`, `kepu.f`), erreur sur l'inclinaison | `PLA(I)%DINC` n'est calculée que dans la branche `RADVEL` : sans vitesses radiales, `DIDP` est rempli mais `ERR_COMB` n'est jamais appelée, et l'erreur affichée sur i vaut 0 | Non corrigée. Constaté le 2026-11-09 dans le dépôt Gaia (BH3 en Levenberg-Marquardt seul). Sans effet sur le MCMC : seul l'affichage du Levenberg-Marquardt (choix 4) est touché. |
| B17 | `MRQELEMENTS` (`kepellip.f`, `kepu.f`), deux branches | `DIDP(4)` et `DIDP(5)` au lieu de `DIDP(DKAL+4)` et `DIDP(DKAL+5)` : pour la deuxième orbite et les suivantes, l'erreur sur i est calculée avec les paramètres de la première | Non corrigée. Constaté le 2026-11-09 en relisant B16. Même portée : affichage du Levenberg-Marquardt seulement. |
| B18 | Programmes principaux (`CHARACTER*80 FILES`), `WRITE_DISTRIB` et `WRITE_DISTRIB_DOUBLING` (`CHARACTER*80 LINE`), `READ_DATAFILE` (`STATUS='UNKNOWN'`) | Les noms de fichiers sont coupés à 80 caractères. Le fichier de données au nom coupé est ouvert en `UNKNOWN`, donc créé vide s'il n'existe pas | **Corrigée par M020.** Avec un chemin de plus de 80 caractères (l'interface écrit des chemins absolus), le programme s'arrête sur `End of file` à la lecture des données, sans dire pourquoi. Le `gzip` des résultats avait la même limite. Constaté le 2026-11-10 en testant l'interface Gaia ; déjà connu d'Antoine. |
| B19 | `MRQRENEW` (`kepellip.f`, `kepu.f`), redémarrage du L-M après une solution impossible | Pendant le premier passage sans jitter des vitesses (`MRQMIN(NPAR-1,...)`), la masse de l'étoile est en `P(NPAR-1)` = `P(NEL*NPLA+2)`. `MRQRENEW` y écrit ln(jitter) et met la masse en `P(NPAR)`, hors du tableau du passage (`P` déclaré `DIMENSION(NPAR)`) | **Corrigée par M022.** Vitesses de l'étoile avec jitter et prior fixe sur sa masse : après un redémarrage, la masse de l'étoile vaut exp(ln σ_jit) (0,0195 Msun pour un jitter initial de 0,01 km/s) ; le prior fixe donne alors une pénalité de 10¹² que le L-M ne rattrape plus, et la solution finale est fausse. Constaté le 2026-11-10 sur HD 114762 (vitesses de Lick seules, χ² = 2,3·10⁶ au lieu de 119). Sans jitter, pas d'effet. |

**Fragilité constatée (corrigée par M018) :** un job tué pendant l'écriture de `dump.dat` laissait un fichier illisible à la reprise (`End of file`), car le fichier était écrit sur place.

**À savoir, sans être une anomalie du calcul :** avec une orbite et sans vitesses de l'étoile, les masses individuelles ne sont pas ajustées. La colonne `m` contient alors la masse totale, dans l'unité de la planète, et `m0` la valeur d'entrée. Une valeur comme 1079 pour une planète en Mjup vient de là (1,03 Msun).

**Déjà corrigé avant ce journal :** une colonne `m` infinie apparaît dans les runs continués avec un ancien binaire : `PLA%MUNIT` n'était pas gardé dans le fichier de sauvegarde. Le commit `83ee116` du 2025-12-12 l'a corrigé.

Correction de B5 et B6 essayée sur une copie : `kepu.f` retrouve alors le χ² du modèle indépendant dans tous les cas de test. La correction de B1 a été essayée de la même façon, sur le cas simulé relatif + vitesses radiales.

## Fonctionnalités à ajouter

| Code | Fonctionnalité | Pourquoi |
|---|---|---|
| F1 | Vitesses radiales absolues pour l'étoile **et** pour les compagnons : un seul type de donnée avec une colonne « corps » (0 pour l'étoile, k pour le compagnon k) et la vitesse systémique `v0` commune. L'actuel type 2 en est le cas « étoile seule » ; un ancien fichier sans colonne « corps » reste lu comme des vitesses de l'étoile. Le type 3 (vitesses relatives compagnon moins étoile) reste inchangé. | L'amplitude absolue d'un compagnon vaut `K_p = (2π/P)·m⋆/(m⋆+m_p)·a sin i/√(1−e²)` : elle dépend du rapport des masses, contrairement à la vitesse relative. Avec l'étoile et le compagnon mesurés en absolu, `K⋆/K_p = m_p/m⋆`. Demandé par Antoine (2026-11-08). Le manuscrit décrit déjà les deux formes (`observations.tex`), alors que le code ne gère que la forme relative pour les compagnons. |
| F2 | Interface d'analyse : pouvoir choisir l'orbite représentative entre le meilleur χ² (actuel, `BestOrbit.py`) et le MAP (déjà lu par `TransferData.py`, jamais utilisé), ou afficher les deux. | Le meilleur χ² ignore les priors ; le MAP en tient compte, ce qui compte avec un prior informatif (par exemple une masse serrée). Depuis M010, M011 et M013, le MAP décrit bien le prior échantillonné par le MCMC. Amélioration secondaire : les résultats restent les distributions. Proposé le 2026-11-08. |
| F3 | Suivre un run en temps réel. **En partie fait par M017** : fichier provisoire `xxx_partial.dat.gz` (positions récentes des chaînes avant la récolte, réserve ensuite), lisible par l'interface d'analyse, et état du run dans `xxx_partial.txt`. Reste : un journal de suivi (R, TS et χ² de chaque chaîne au fil du run, avec chaîne et ordre des positions pour tracer les traces), un bouton « rafraîchir » et une mention « non convergé » dans l'interface d'analyse, et la durée réelle des écritures à mesurer sur un gros run. | Avant M017, seul `out.txt` (pas, χ² des chaînes, critère R) renseignait pendant le run. Demandé par Antoine (2026-11-08) : voir où en est le run sans attendre la fin. |

## Modifications

## M024 — Directions du MCMC recalculées tôt après la chauffe

- **Date :** 2026-11-10
- **Fichiers :** `Code/mcmc/sub/commons.f` (constantes `NDIAG1` = 10 001 et `NDIAG2` = 100 001), `Code/mcmc/sub/mcmc.f` (`MCMCP` : condition du recalcul des directions, un bloc repéré ; et un appel à `CHANGE_DIRS` ajouté sans repère dans le bloc M015 de fin de chauffe, qui est du code ajouté par M015).
- **Changement :** les directions de proposition (`CHANGE_DIRS`, axes propres de la covariance des chaînes) sont recalculées en fin de chauffe, puis quand `LONG` franchit 10 001 et 100 001 balayages après la chauffe, en plus du recalcul tous les `NDIAG` = 500 000 balayages d'origine. La condition porte sur le franchissement du seuil (`LONG ≥ seuil` et `LONG − NGROUP < seuil`), pas sur une valeur exacte de `LONG`.
- **Pourquoi :** avant M024, les chaînes gardaient les axes de départ (directions des paramètres eux-mêmes) jusqu'à 500 000 balayages. Quand les paramètres sont fortement corrélés, les pas le long de ces axes sont minuscules et la convergence se fait attendre pendant des millions de balayages. Recalculer tôt, dès que les chaînes ont exploré un peu, aligne les pas sur les corrélations.
- **Mesures (systèmes synthétiques, `--rapide`, 2026-11-09 et 10) :**
  - `arc_court` : même distribution, ESS de l'échantillon final 3 à 6 fois plus grand qu'avec le code d'avant.
  - `arc_tres_court` (cas le plus difficile) : pas de gain ni de perte visibles (début de la récolte à 6,0 et 7,1 millions de balayages contre 7,2).
  - Variante écartée : recalcul à chaque doublement de `LONG` (fin de chauffe, 1 000, 2 000, 4 000…). Sur `arc_tres_court`, la récolte commence plus tard (9,0 et 18,6 millions de balayages) : des recalculs trop fréquents gênent l'exploration.
- **Vérification :** `rel_vr_planete` (`--rapide`) : recalculs bien faits en fin de chauffe, vers 10 000 et 100 000 balayages ; run terminé en 1 min 20 s ; tous les éléments à moins de 1,8 σ de la vérité. `make test` passe.
- **Variante écartée (2026-11-10) : directions gelées pendant la récolte** (aucun recalcul une fois `INFER` vrai, y compris celui tous les 500 000 balayages). Sur `arc_court` et `rel_vr_planete`, la récolte commence juste après la chauffe, donc les recalculs de M024 tombent pendant la récolte. 2 runs par variante et par système : mêmes distributions (écarts à la vérité identiques à 0,05 σ près), mais ESS minimal plus faible avec le gel (`rel_vr_planete` : 3 000 et 7 400 contre 16 400 et 17 500 ; `arc_court` : 2 600 et 5 000 contre 4 700 et 5 400). Les recalculs pendant la récolte aident donc, sans biais visible : M024 est gardé tel quel. Avec M024, `arc_court` se termine en 200 000 à 300 000 balayages.
- **Effet sur les résultats :** même distribution a posteriori ; convergence plus rapide, ou inchangée, selon le cas.

## M023 — Icône et logo d'ORACLE (sans repère)

- **Date :** 2026-11-10
- **Fichiers :** `Interface/Items/oracle-app.icns`, `.ico`, `.svg` (icône dessinée par Antoine, pour macOS, Windows et Linux) et `oracle-app.png` (512 px, tiré du `.icns`) et `oracle-app-macos.png` ajoutés, ainsi que le logo d'Antoine : `oracle-logo-fond-blanc.svg` (vectoriel) et `.jpg`, `oracle-logo-sans-fond.png` (encre bleu nuit sur fond transparent, pour fond clair), `oracle-logo-sombre.png` (encre claire sur fond bleu nuit) ; images de 1600 × 334 px ; `Interface/Main.py` (code d'Antoine).
- **Changement :** `Main.py` donne à l'application son nom (`setApplicationName('Oracle')`) et son icône (`setWindowIcon`) : `oracle-app-macos.png` sur macOS, `oracle-app.png` ailleurs.
- **Version macOS :** sur Mac, le dessin d'une icône occupe 824 px sur 1024, centré, avec une marge transparente (grille d'Apple, mesurée sur Zoom, Calendrier et VS Code : 80,5 % et 9,8 % de marge de chaque côté) ; le dessin d'origine remplit tout le carré et paraissait plus gros que les autres icônes du Dock. `oracle-app.icns` (toutes les tailles, par `iconutil`) et `oracle-app-macos.png` sont refaits avec le dessin réduit à 80,5 % et centré ; `.ico`, `.svg` et `oracle-app.png` gardent le dessin d'origine.
- **Pourquoi :** l'icône apparaît sur les fenêtres, dans la barre des tâches (Windows, Linux) et dans le Dock (macOS), quel que soit le moyen de lancement. Les formats `.icns`, `.ico` et `.svg` serviront aux lanceurs (`Oracle.app`, raccourci Windows, `oracle.desktop`) que créera le futur script d'installation.
- **Vérification :** interface lancée hors écran par un `Oracle.app`, sans erreur.
- **Effet sur les résultats :** aucun.

## M022 — Corrige B19 : redémarrage du L-M au premier passage sans jitter

- **Date :** 2026-11-10
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `Code/mcmc/sub/kepu.f`, routine `MRQRENEW` (deux blocs chacun).
- **Changement :** `P` est déclaré `DIMENSION(N)` (dimension du L-M en cours) au lieu de `NPAR` ; la masse de l'étoile est remise en `P(N)` au lieu de `P(NPAR)` ; le jitter des vitesses n'est remis que si `N = NPAR`, c'est-à-dire s'il fait partie du passage.
- **Pourquoi :** voir B19. La masse est toujours le dernier paramètre du passage, `P(N)`, au premier comme au second.
- **Vérification :** HD 114762, 74 vitesses de Lick (Kane et al. 2011, ApJ 735, L41), jitter, masse de l'étoile fixée à 0,83 Msun : χ² = 118,8 au lieu de 2,3·10⁶ ; P = 83,9145 ± 0,0032 j (Kane : 83,9151 ± 0,0030), e = 0,3385 ± 0,0050 (0,3354 ± 0,0048), ω = 201,97° (201,28°), m sin i = 11,0 MJ (10,98). Sans jitter, résultat inchangé. `make test` passe.
- **Remarque :** avec les vitesses de l'étoile seules, l'ω affiché par ORACLE est égal à celui de Kane et al., qui suivent la convention habituelle des vitesses radiales (ω de l'orbite de l'étoile, z vers l'arrière du ciel). Ce que cela dit de la convention d'ORACLE avec l'astrométrie (ω du compagnon ou de l'étoile, sens de z) reste à vérifier.
- **Effet sur les résultats :** seulement pour les runs avec vitesses de l'étoile et jitter dont le premier passage du L-M a redémarré ; le MCMC part ensuite de la solution du L-M.

## M021 — Interface : point décimal, commentaires du fichier de lancement, a priori de masse (sans repère)

- **Date :** 2026-11-10
- **Fichiers :** `Interface/Main.py`, `Interface/NewSimu/PriorMass.py`, `Interface/NewSimu/WindowSetNewSimu.py` (code d'Antoine, sans repère).
- **Changement :**
  1. `Main.py` impose le point décimal (`QLocale.setDefault(QLocale(QLocale.Language.C))`) avant de créer les fenêtres.
  2. Les commentaires des choix 3 et 4 du `start.sh` étaient inversés (« LM only » pour 3, « MCMC only » pour 4). Les codes écrits étaient justes.
  3. Les champs des a priori de masse (moyenne, écart-type, bornes, valeur) passent de 2 à 4 décimales.
  4. L'unité d'un a priori de masse s'écrivait deux fois (`msMsun`) : seule `ms` ou `mj` reste.
  5. `Create startup files` rend `start.sh` exécutable, comme le faisait déjà `Start the simulation`.
- **Pourquoi :** (1) sur un système réglé en français, Qt peut écrire `1,66` : Fortran lit alors deux nombres, 1 et 66, sans erreur. Constaté le 2026-11-10 avec l'interface lancée hors écran ; les `start.sh` réels d'Antoine avaient des points. (3) Un écart-type de 0,001 Msun était impossible à saisir. (4) Fortran ne lit que deux caractères pour l'unité : sans effet, mais trompeur.
- **Vérification :** `start.sh` produit hors écran avec points décimaux, a priori `0.7600 ms`.
- **Effet sur les résultats :** aucun pour les runs déjà faits (leurs `start.sh` avaient des points).

## M020 — Corrige B18 : noms de fichiers de plus de 80 caractères

- **Date :** 2026-11-10
- **Fichiers :** `Code/mcmc/main/astrom_mcmco.f`, `Code/mcmc/main/astrom_univ_mcmco.f` (`FILES` en `CHARACTER*512`), `Code/mcmc/sub/io.f` (`LINE` de `WRITE_DISTRIB` et `WRITE_DISTRIB_DOUBLING` en `CHARACTER*600`, ouverture du fichier de données dans `READ_DATAFILE`). Sans repère, dans des lignes ajoutées par M017 et M018 : `BASE` passe de 256 à 512 caractères et `LINE` de 600 à 1 600 dans `WRITE_PARTIAL` (`mcmc.f`), `LINE` de 600 à 1 100 dans `WRITE_DUMP` (`io.f`), pour contenir deux ou trois noms longs.
- **Changement :** noms de fichiers jusqu'à 512 caractères. Le fichier de données est ouvert en `STATUS='OLD'` : s'il n'existe pas, le programme s'arrête avec `Data file not found : <nom>`.
- **Pourquoi :** voir B18. Les sous-programmes reçoivent déjà les noms en `CHARACTER*(*)` : seules les déclarations des programmes principaux et des commandes `gzip`/`mv` fixaient une longueur.
- **Dump :** les noms y sont écrits et relus en `(a)` : un dump écrit avant M020 se relit sans changement.
- **Vérification :** `make test` passe. GG Tau en Levenberg-Marquardt avec un chemin de données de 167 caractères : χ² + Fmap = 25,0651, comme la référence. Fichier de données absent : arrêt avec le message.
- **Effet sur les résultats :** aucun.

## M019 — Arrêt sur le nombre de tirages indépendants (ESS) au lieu de TS

- **Date :** 2026-11-09
- **Fichiers :** `Code/mcmc/sub/commons.f` (constantes `ESSMIN` et `PESSMIN`), `Code/mcmc/sub/mcmc.f` (`GELMAN_RUBIN` : critère de convergence, appel à la fin du second remplissage ; nouvelle routine `ESS_RESERVE` en fin de fichier).
- **Changement :**
  1. **Convergence sur R seul** : `OK` demande R < `RCONV` = 1,01 pour chaque paramètre, sans `TS` > `TSCONV` (constante gardée, plus utilisée). `TS` reste affiché dans `out.txt`. Le critère d'inférence (`RINF`, `TSINF`) ne change pas.
  2. **ESS de la réserve** : quand le second remplissage de M015 remplit la réserve, `ESS_RESERVE` calcule pour chaque paramètre d'ajustement le nombre de tirages indépendants, par moyennes de lots : pour chaque chaîne, `n` = `NINFER` solutions gardées dans l'ordre, lots de taille `sqrt(n)`, temps d'autocorrélation `τ = taille × variance des moyennes de lots / variance`, ESS = somme sur les chaînes de `n/τ`. Une chaîne immobile compte pour 1. Les ESS sont écrits dans `out.txt`.
  3. Si le plus petit ESS atteint `ESSMIN` = 400, le run s'arrête comme avant (réserve pleine et 5 contrôles convergés de suite). Sinon, la probabilité de récolte est abaissée pour viser 2 × `ESSMIN` (`p ← min(p, 1) × ESS / (2 × ESSMIN)`) et la réserve est remplie à nouveau (`NMOD` = 0). Cette probabilité peut descendre sous `MINPROB` = 0,1, jusqu'à `PESSMIN` = 0,01 (complété le 2026-11-09 : sur `arc_tres_court`, la probabilité était déjà à `MINPROB` et l'ESS ne pouvait pas remonter). Espacer la récolte ne biaise pas l'échantillon, cela allonge seulement le remplissage (environ `NINFER` × 100 / p balayages, soit au plus 500 millions de balayages pour `NINFER` = 50 000). Si la probabilité est déjà à `PESSMIN`, l'échantillon est gardé avec un avertissement. La probabilité baisse au moins de moitié à chaque nouveau remplissage, donc la boucle se termine.
- **Références :**
  - R : Gelman & Rubin 1992, *Stat. Sci.* 7, 457 (statistique déjà utilisée par Hervé).
  - TS : Ford 2006, *ApJ* 642, 505 (`Tz`, aussi utilisé par RadVel avec un seuil de 1 000).
  - Seuils R < 1,01 et ESS > 400 : Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021, *Bayesian Anal.* 16, 667.
  - Moyennes de lots de taille `sqrt(n)` : Flegal & Jones 2010, *Ann. Stat.* 38, 1034.
- **Pourquoi :**
  - Avec les formules de `GELMAN_RUBIN`, `TS` ≈ `NCH`/(R² − 1) : avec 8 chaînes, `TS` > 1 000 revient à R < 1,004. C'était le vrai critère d'arrêt (GG Tau : environ 300 millions de pas, 30 % du run, entre R < 1,01 et l'arrêt).
  - `TS` mesure l'écart entre chaînes, pas la qualité de l'échantillon gardé. Celle-ci dépend de la durée du second remplissage, et l'ESS la mesure directement.
  - Test sans `TS` (`TSCONV` = 0, `arc_court`, 2 runs, 2026-11-09) : mêmes distributions qu'avec `TS` et qu'avec un remplissage 4 fois plus long (écarts < 0,05 σ). ESS minimal de 700 à 1 000 (Ω, ω), contre 1 300 avec `TS` et 14 500 avec le facteur 4.
- **Vérification (2026-11-09, `Tests/build_synthetic/ess/`, `--rapide`) :**
  - `arc_court`, 2 runs : ESS minimal 2 034 et 1 728, arrêt au premier remplissage. Distributions identiques aux runs précédents (tous les paramètres à moins de 1,4 σ de la vérité, écarts entre runs < 0,05 σ).
  - `rel_vr_planete` : ESS minimal 18 079, vérité retrouvée (moins de 1,8 σ, comme pour M015).
  - `arc_court` avec `ESSMIN` monté à 3 000 : premier ESS 2 165, donc nouveau remplissage à la probabilité 0,36 ; ESS final 7 925, même distribution.
  - Comparaison avec un calcul indépendant en Python (somme de Geyer) sur les éléments physiques : même ordre de grandeur. Le calcul du Fortran porte sur les paramètres d'ajustement : il donne environ 20 % de plus pour Ω et ω, qui y sont mélangés à d'autres quantités.
  - `make test` passe.
- **Effet sur les résultats :** arrêt plus tôt quand l'échantillon est déjà assez riche ; remplissage plus long, avec une récolte plus espacée, quand il ne l'est pas. Choix de méthode à présenter à Hervé.

## M018 — `dump.dat` écrit sous un nom temporaire puis renommé

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/io.f`, routine `WRITE_DUMP` (trois blocs repérés).
- **Changement :** le dump est écrit dans `dump.dat.tmp`, puis renommé en `dump.dat` (`mv -f`, par `CALL SYSTEM` comme le `gzip` des résultats).
- **Pourquoi :** un job tué pendant l'écriture (fin de durée sur le cluster, par exemple) laissait un `dump.dat` tronqué et la reprise échouait (`End of file`). Le renommage est atomique : le dump précédent reste intact jusqu'à ce que le nouveau soit complet.
- **Vérification :** run arrêté puis repris depuis le dump, aucun `dump.dat.tmp` restant. `make test` passe.
- **Effet sur les résultats :** aucun.

## M017 — F3 : fichier provisoire pour suivre un run en temps réel

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/commons.f` (constantes `NPARTIAL`, `TPARTIAL`, `FPARTIAL`, compteur `NHARV`), `Code/mcmc/sub/mcmc.f` (deux appels dans `MCMCP`, une ligne dans `GELMAN_RUBIN`, nouvelle routine `WRITE_PARTIAL` en fin de fichier).
- **Changement :** à chaque contrôle de Gelman-Rubin, `WRITE_PARTIAL` garde les `NPARTIAL` = 1 000 dernières positions de chaque chaîne. Toutes les `TPARTIAL` = 300 s (ou plus si l'écriture est longue : l'intervalle vaut au moins 20 fois la durée de la dernière écriture, donc au plus 5 % du temps de calcul), elle écrit `xxx_partial.dat.gz` (`xxx` = nom du fichier de résultats sans `.dat`), au même format que les résultats : les positions récentes des chaînes avant la récolte, toute la réserve ensuite. Pendant la récolte, la réserve n'est réécrite que si au moins `FPARTIAL` = 5 % de ses solutions sont nouvelles depuis la dernière écriture (compteur `NHARV`, incrémenté à chaque solution gardée dans `GELMAN_RUBIN`) : inutile de réécrire un fichier qui n'a presque pas changé, par exemple dans les longues fins de run. Le fichier est écrit sous un nom temporaire puis renommé, pour qu'un lecteur ne tombe jamais sur un fichier en cours d'écriture. `xxx_partial.txt` donne l'état : nombre de balayages, récolte commencée (1) ou non (0), nombre de solutions gardées, nombre de positions récentes. Les deux fichiers sont retirés à la fin du run ; ils restent si le run est interrompu. Seule la version parallèle (`MCMCP`, la seule utilisée) est concernée.
- **Pourquoi :** voir où en est un run sans attendre la fin (F3, demandé par Antoine). Avant, le fichier de résultats n'était écrit pendant le run qu'une fois la réserve pleine et tous les 10 millions de pas.
- **Vérification :** runs `rel_vr_planete` avec `TPARTIAL` réduit à 0,5, 2 et 5 s : fichier écrit avant la récolte (positions récentes) puis pendant (solutions gardées, réserve entière pendant le second remplissage de M015), fichiers retirés à la fin, aucun fichier temporaire restant. L'interface d'analyse (`TransferData.TransfertSimu`) lit le fichier provisoire sans changement. `make test` passe.
- **Effet sur les résultats :** aucun. Coût : au plus 5 % du temps de calcul. Une réserve de 400 000 orbites avec deux orbites fait plusieurs centaines de Mo de texte avant compression : l'intervalle s'allonge alors de lui-même. Sans cette règle, avec une écriture toutes les 0,5 s, le run test était 30 fois plus long.
- **Sortie d'Hervé inchangée :** `GELMAN_RUBIN` réécrit toujours `results.dat.gz` tous les 10 millions de pas quand la réserve est pleine. Elle s'interrompt pendant le second remplissage de M015. Les deux sorties font en partie doublon : à voir avec Hervé.
- **Reste à faire :** journal de suivi (R, TS et χ² de chaque chaîne au fil du run, pour les traces), bouton « rafraîchir » dans l'interface d'analyse, mention « non convergé » à l'affichage.

## M016 — Corrige B15 : plus de contrôles sans pas une fois la convergence atteinte

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/mcmc.f`, routine `MCMCP` (un bloc repéré).
- **Changement :** une fois R < 1,01, l'intervalle entre deux contrôles vaut `NGROUP` (100 balayages). L'ancienne formule donnait le même résultat au-delà de 5 000 balayages, mais 0 en deçà.
- **Pourquoi :** avec un intervalle nul, les chaînes ne bougeaient plus et les mêmes positions étaient gardées à chaque contrôle (B15).
- **Vérification :** voir M015 (5 runs sur 5 corrects, environ 14 300 solutions distinctes sur 20 001).
- **Effet sur les résultats :** seuls les runs qui convergent avant 5 000 balayages sont touchés, donc des cas faciles. Les runs longs ne changent pas.

## M015 — Corrige B12 : départs dispersés, chauffe, réserve remplie à nouveau après la convergence

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/commons.f` (`NBURN`, `KDISP`, `KREFILL`, variable `REFILL`), `Code/mcmc/sub/mcmc.f` (`MCMCP`, `INIT_MCMC`, `GELMAN_RUBIN`), `Code/mcmc/sub/io.f` (`WRITE_DUMP`, `READ_DUMP`).
- **Changement :**
  1. **Départs dispersés** (`INIT_MCMC`) : la chaîne 1 part de la solution du L-M ; les autres en sont décalées de `KDISP` = 2 pas `BETA` (tirage gaussien) le long de chaque direction, avec au plus 100 essais pour un point accepté par `RATIO_PRIORS`.
  2. **Chauffe** (`MCMCP`) : pendant `NBURN` = 1 000 balayages, ni Gelman-Rubin ni récolte. À la fin, moyennes et covariances cumulées sont remises à zéro et `LONG` repart de 1.
  3. **Second remplissage** (`GELMAN_RUBIN`) : la récolte commence toujours à R < 1,1 (`RINF` inchangé). Au premier passage à la convergence (R < 1,01 et `TS` > 1 000), `NMOD` est remis à 0 : la réserve est réécrite dans l'ordre, et la probabilité de récolte est recalculée pour que ce remplissage dure `KREFILL` = 1 fois le temps mis pour converger. Le run s'arrête comme avant (réserve pleine et 5 contrôles convergés de suite) : toutes les solutions gardées datent alors d'après la convergence. Pas de seconde remise à zéro si R repasse au-dessus de 1,01.
  4. **Reprise** : `REFILL` est écrit à la fin de `dump.dat` et relu. Un ancien fichier sans cette ligne reste lisible ; la remise à zéro est alors refaite.
- **Pourquoi :** les 8 chaînes partaient du même point, d'où un R ≈ 1 au premier contrôle (B12). Et avec le remplacement aléatoire d'Hervé, rien ne garantissait que les solutions récoltées avant la convergence (entre 1,1 et 1,01) soient remplacées avant l'arrêt : dans un run qui converge vite, elles formaient l'essentiel de l'échantillon. Choix d'Antoine (2026-11-08) : garder la récolte dès 1,1 pour l'aperçu, mais un résultat final tiré après la convergence. Facteur 1 plutôt que le 4 d'Hervé pour le second remplissage : après la convergence, les chaînes « oublient » leur position bien plus vite qu'elles n'ont mis à converger ; le facteur 4 aurait rendu les runs environ 5 fois plus longs, contre environ 2 fois au plus.
- **Vérification :** `rel_vr_planete` en classique, 5 runs (`--rapide`) : avant, 0 run correct sur 5 (P entre 1,1 et 3,1 σ, distributions trop étroites) ; après M015 et M016, 5 sur 5 donnent la même distribution (P à 1,3 σ de la vérité, tous les paramètres à moins de 1,8 σ, barres symétriques), environ 14 300 solutions distinctes sur 20 001. Reprise depuis un `dump.dat` écrit pendant le second remplissage : pas de seconde remise à zéro, même distribution ; même fichier sans la ligne `REFILL` : remise à zéro refaite, même distribution. `make test` passe.
- **Effet sur les résultats :** les runs qui convergent vite donnent maintenant des distributions correctes. Les runs sont plus longs, d'au plus le temps mis pour converger. Les longs runs réels devraient peu changer : leurs anciennes solutions d'avant la convergence étaient déjà presque toutes remplacées (à vérifier sur GG Tau). Choix de méthode à présenter à Hervé.

## M014 — Prior logarithmique par défaut sur la masse de l'étoile

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `Code/mcmc/sub/kepu.f` : `COMPUTE_MAP`, `RATIO_PRIORS` et `MRQMAPASTROM` (un bloc repéré `M014` dans chacune, dans les deux fichiers). Le bloc de `RATIO_PRIORS` remplace le second bloc de M013. `Tests/reference/reference.json` refait.
- **Changement :** quand `m0` est ajustée (`MULTIPLA`), son prior par défaut est maintenant logarithmique (`1/m0`) au lieu de plat. `COMPUTE_MAP` ajoute `-ln m0`. Le L-M (`MRQMAPASTROM`, où `FMAP = -2 ln prior`) ajoute `+2 ln m0` et sa dérivée `2` par rapport à `ln m0`. Le MCMC n'a plus de facteur `m0'/m0` : un prior logarithmique en `m0` est plat en `ln m0`, son paramètre.
- **Pourquoi :** `m0` est une grandeur d'échelle, comme `a`, `P`, `q` et `P_q`, qui ont toutes un prior logarithmique dans ORACLE. Le prior plat en `m0` (aucun terme dans `COMPUTE_MAP`) était la seule exception et favorisait les grandes masses. Choix d'Antoine (2026-11-08), à présenter à Hervé : c'est un changement de convention, pas la correction d'une erreur.
- **Vérification :** `make test` : le rapport du MCMC suit toujours exactement `COMPUTE_MAP` et le jacobien dans les 22 cas, le gradient du L-M suit les différences finies, et `MAP = χ² + Fmap` est vérifié. La référence a été refaite : dans les 11 cas où `m0` est ajustée, le MAP augmente de 8,105645 et le critère du L-M baisse de 16,211290, soit exactement `-ln m0` et `2 ln m0` pour `m0` = 1,02 Msun en unités internes. Les χ² ne changent pas. Les autres différences de la référence (cas universels) sont des arrondis de 1e-16 à 1e-13 en relatif, venus de M011 (calcul de `s` à partir de `λ`), très en dessous de la tolérance de 1e-10.
- **Effet sur les résultats :** les runs où `m0` est ajustée sans prior de l'utilisateur la favorisent un peu moins grande qu'avant. L'effet est faible quand `m0` est bien contrainte par les données ou par un prior. Les valeurs de MAP affichées pour ces runs sont décalées d'une constante (`-ln m0` varie peu). Rien ne change pour les runs à une orbite sans vitesses de l'étoile.

## M013 — Corrige B11 : le prior de masse multiplie les priors par défaut, dans le MCMC comme dans le MAP

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `Code/mcmc/sub/kepu.f`, routine `RATIO_PRIORS` (deux blocs repérés `M013` dans chaque fichier). Tests : `Tests/jacobian.py` (plus d'écart attendu) et `Tests/run_tests.py`.
- **Changement :** pour chaque prior de masse `u(μ)`, le rapport du MCMC ne contient plus que `u(μ')/u(μ)`. Les facteurs `M'/M` (pour chaque masse dynamique concernée) et `m0'/m0` (si `m0` est concerné) sont retirés. Un seul facteur `m0'/m0` est ajouté après la boucle quand `MULTIPLA` est vrai : c'est le jacobien de `ln m0`, puisque le MCMC échantillonne `ln m0` et que `COMPUTE_MAP` suppose un prior plat en `m0`. Avant, ce facteur était fourni, sans que ce soit son rôle, par le prior de masse positive sur `m0`.
- **Pourquoi :** `RATIO_PRIORS` et `COMPUTE_MAP` ne décrivaient pas le même prior. Le MCMC **remplaçait** le prior logarithmique par défaut par le prior de l'utilisateur le long de la masse (intention du commentaire du code), alors que le MAP le **multipliait**. Les facteurs étaient en plus appliqués aux priors de masse positive (type 6) et comptés plusieurs fois pour une même masse. Choix d'Antoine (2026-11-08), à présenter à Hervé : « multiplier ». C'est le seul choix bien défini avec plusieurs priors sur des combinaisons de masses différentes, et il se résume ainsi : « les priors de l'utilisateur s'ajoutent aux priors par défaut d'ORACLE ». Avec un prior étroit (par exemple 1,55 ± 0,03 Msun), les deux choix donnent presque le même résultat.
- **Vérification :** `make test` : le rapport des priors suit maintenant exactement `COMPUTE_MAP` et le jacobien calculé numériquement dans les 22 cas (classique et universel, avec et sans masses ajustées, avec prior sur l'étoile, avec jitter). Plus aucune anomalie connue dans les tests. Référence inchangée (elle porte sur χ² et MAP, pas sur le rapport).
- **Effet sur les résultats :** les runs où les masses sont ajustées (vitesses de l'étoile, plusieurs orbites ou Gaia) échantillonnent maintenant le prior qu'annonce le MAP, qui ne favorise plus les grandes masses. Sur un run réel à deux orbites, estimé par repondération : `m0` passe d'environ 0,563 à 0,536 Msun, les orbites externes changent peu. Les runs à une orbite sans vitesses de l'étoile ne changent pas.

## M012 — Corrige B9 : refus d'un prior sur la masse de l'étoile seule quand les masses ne sont pas ajustées

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/io.f`, routine `ENTER_PRIOR` (un bloc repéré `M012`). Interface : `Interface/NewSimu/WindowSetNewSimu.py`, nouvelle méthode `LockMassPriorsOrNot` (sans repère).
- **Changement :** à la lecture d'un prior de masse, si les masses individuelles ne sont pas ajustées (`MULTIPLA` faux : une orbite, sans vitesses de l'étoile ni Gaia) et si le prior demande `m0` séparément (`BCOF(0)` non nul, c'est-à-dire des coefficients différents pour `m0` et `m1`), le programme s'arrête avec un message qui indique de donner le même coefficient à `m0` et `m1` (masse totale, par exemple `1 1`).
- **Interface :** quand les masses individuelles ne sont pas ajustées (deux corps, sans vitesses de l'étoile ni astrométrie absolue), le coefficient du compagnon (`m2` dans la fenêtre) est grisé et recopie celui de l'étoile (`m1`) : seul un prior sur la masse totale peut être saisi. Il redevient modifiable dès qu'on coche les vitesses de l'étoile ou l'astrométrie absolue, ou qu'on passe à plus de deux corps. Vérifié en ouvrant la fenêtre sans affichage (Qt `offscreen`).
- **Pourquoi :** dans ce cas, seule la masse totale `m0 + m1` est un paramètre. `COMPUTE_MAP` et `RATIO_PRIORS` ignoraient le terme en `m0` et évaluaient le prior sur une masse nulle : il n'agissait pas, sans message (cas de 51 Eri). Les vitesses de la planète ne changent rien : ORACLE ne connaît que des vitesses relatives à l'étoile (type 3), qui ne dépendent que de la masse totale. Choix d'Antoine (2026-11-08), à présenter à Hervé, parmi trois options : refuser, appliquer le prior à la masse totale avec un avertissement, ou seulement avertir.
- **Vérification :** `make test` : le prior `1 0` sur le cas `rel` est refusé avec le message, dans les deux versions. Les priors sur la masse totale (`1 1`) et les cas où les masses sont ajustées agissent toujours comme attendu. Référence inchangée.
- **Effet sur les résultats :** aucun résultat valide ne change. Un fichier d'entrée avec un tel prior, qui donnait un run sans prior, est maintenant refusé.
- **Suite prévue :** vitesses radiales absolues des compagnons (voir la section « Fonctionnalités à ajouter »). Avec elles, les masses individuelles deviennent ajustables sans vitesses de l'étoile, et un prior sur `m0` seul redevient valable.

## M011 — Paramètres 6 et 7 de la version universelle : (cos λ, sin λ)/P_q (corrige aussi B10)

- **Date :** 2026-11-08
- **Fichiers :** `Code/mcmc/sub/kepu.f` : `COMPUTE_MAP` (un bloc), `ELEMENTS` (trois blocs), `RATIO_PRIORS` (quatre blocs), nouvelle routine `NU_TO_S` en fin de fichier. `Code/mcmc/main/astrom_univ_mcmco.f` : deux blocs (déclarations, point de départ du MCMC). Tous repérés `M011`. Tests : `Tests/driver.f`, `Tests/longitude_univ.f` et `Tests/jacobian.py` adaptés.
- **Changement :** dans le MCMC universel, les paramètres 6 et 7 étaient `ln P_q` et `s·f(e)` (variable universelle à `T0` multipliée par `f(e) = 10e²/(1+10e²)`). Ils deviennent `(cos λ, sin λ)/P_q`, avec `λ = ϖ + ν` la longitude vraie à `T0`, comme dans la version classique. `λ` est définie pour toutes les coniques, et `s` s'en déduit sans inversion numérique (`NU_TO_S`, formules de l'ellipse, de l'hyperbole et de la parabole). Le rapport des priors (`RATIO_PRIORS`) utilise le jacobien de ces nouveaux paramètres : `e·(1+e·cos ν)²/((1+e)^(3/2)·P_q³)` par orbite, sans le facteur `(q/q')(P_q/P_q')` de B10. Un point d'essai hyperbolique au-delà de l'asymptote (`1+e·cos ν ≤ 0`) est refusé. Le Levenberg-Marquardt garde ses paramètres (`ln P_q`, `s·f(e)`) : seule la conversion vers le point de départ du MCMC change, dans `astrom_univ_mcmco.f`.
- **Pourquoi :** décision d'Antoine (2026-11-06) après les tests de `Tests/parametrisation/` : `s·f(e)` est un mauvais paramètre d'ajustement (queue très lourde, MCMC coincé près du périastre et déclaré convergé à tort). Avec lambda, le MCMC retrouve les orbites quasi paraboliques (`arc_parabole`, e de 0,94 à 1,04) et converge 4 à 10 fois plus vite autour du périastre. Voir `Tests/parametrisation/README.md`.
- **Vérification :** le code est identique, aux commentaires près, à la variante validée dans `Tests/parametrisation/` (plus les corrections M006, M008 et M009 faites depuis). `make test` : référence inchangée (χ² et MAP identiques sur les 11 cas universels), et le rapport des priors universel suit maintenant exactement le jacobien calculé numériquement : il ne reste que l'écart B11, comme en version classique.
- **Effet sur les résultats :** les runs universels échantillonnent dans de nouveaux paramètres, avec le prior annoncé par `COMPUTE_MAP` (`1/(q·P_q)`, B10). Les fichiers de sortie contiennent toujours les éléments d'orbite. **Un fichier de reprise (dump) d'un run universel lancé avant M011 n'est plus utilisable** : ses paramètres 6 et 7 n'ont pas le même sens. La version classique n'est pas touchée.

## M010 — Corrige B13 : priors de masse de types 0, 3, 4 et 6

- **Date :** 2026-11-06
- **Fichiers :** `Code/mcmc/sub/mcmc.f`, routines `USER_PRIOR_MAP` (trois blocs) et `USER_PRIOR_RATIO` (deux blocs), repérés `M010`. `Code/mcmc/sub/io.f`, routine `ENTER_PRIOR` (un bloc). Interface : `Interface/NewSimu/PriorMass.py` (sans repère : l'interface est le code d'Antoine). Nouveau test : `Tests/priors.f` et la section « Priors sur les masses » de `Tests/run_tests.py`.
- **Changement :**
  - Type 0 (logarithmique, prior en 1/m) : le rapport du MCMC vaut `m / m_essai`. Il valait `ln m / ln m_essai`, qui ne correspond à aucun prior.
  - Type 3 (uniforme, « Uniform » dans l'interface) : les bornes lues par `ENTER_PRIOR` sont appliquées. Le MCMC refuse une masse hors des bornes (rapport 0) et le MAP la pénalise de 1e30. Le prior reste plat entre les bornes. Le L-M (`FMAP_PRIOR`) n'est pas modifié : il traite toujours ce prior comme une constante, car un mur ne se dérive pas.
  - Type 4 (fixe) : le MAP utilise la même largeur relative (1e-6) que le MCMC et le L-M. Il utilisait 1e-3.
  - Interface : l'entrée « Log », qui envoie le type 2 (log-normal), est renommée « Log-normal » pour ne pas être confondue avec le prior logarithmique (type 0). Le type 0 n'est pas ajouté à l'interface : sans prior, ORACLE est déjà à peu près log-uniforme sur la masse dynamique (ln a et ln P plats), et l'effet d'un prior de plus dépend de B11, à discuter avec Hervé. `ENTER_PRIOR` donne maintenant une unité (1) à tout prior : le type 0, écrit à la main, n'en lisait pas et l'unité restait non définie dans les fichiers de sortie.
  - Type 6 (masses individuelles positives, ajouté automatiquement quand les masses sont ajustées) : la ligne `EX = EXP(X)`, où `X` n'avait pas de valeur, est supprimée. Pour une masse très négative, le terme `-LOG(1+EXP(-S))` est remplacé par son équivalent `S`, pour éviter un `EXP` infini, comme le fait déjà le rapport du MCMC.
- **Pourquoi :** les types 2, 3 et 4 ne pouvaient pas tourner avant M007 (le programme s'arrêtait à leur lecture), le type 0 n'est pas proposé par l'interface, et l'erreur du type 6 donnait un terme proche de 0 pour des masses positives : aucun n'avait été vérifié en fonctionnement.
- **Vérification :** nouveau test, pour chaque type : densité attendue dans le MAP, rapport du MCMC égal à `exp(MAP(m') - MAP(m))`, L-M égal à -2 fois le MAP à une constante près, masse non modifiée par le MAP, masses hors support refusées par le MCMC et pénalisées dans le MAP. Avant la correction : échec des types 0, 3, 4 et 6, succès des types 1 et 2. Après : tous passent. Le test de lecture contient aussi un prior de type 0 écrit à la main (le type, sans ligne de paramètres). `make test` : référence inchangée (aucun cas de référence n'utilise ces types).
- **Effet sur les résultats :** un run avec un prior uniforme borne maintenant vraiment la masse. Un run avec un prior fixe affiche un MAP cohérent avec ce qui est échantillonné. Les runs à masses ajustées affichent un MAP qui pénalise les masses négatives ; leur MCMC ne change pas. Rien ne change pour les priors gaussiens et log-normaux.

## M009 — Corrige B7 : dérivée des vitesses de la planète dans le LM universel

- **Date :** 2026-11-06
- **Fichier :** `Code/mcmc/sub/kepu.f`, routine `MRQVFITP` (un bloc repéré `M009`).
- **Changement :** la dérivée partielle de la vitesse `VY` par rapport à la masse dynamique de l'orbite était divisée par `MTOT`. Elle est maintenant divisée par `MDYN(I)`, comme la ligne voisine pour `VX`.
- **Pourquoi :** `MTOT` ne vaut `MDYN(I)` que si `MULTIPLA` est vrai. Avec une seule orbite et sans vitesses de l'étoile, `MTOT` reste la masse de l'étoile, et la dérivée était fausse dans le rapport `MSTAR/MDYN`. La masse dynamique dépend de q et de P_q, d'où l'erreur sur les dérivées des paramètres 1 et 6. La même ligne dans `MRQVFIT` (vitesses de l'étoile) est juste et n'est pas modifiée : ces vitesses imposent `MULTIPLA` vrai (`io.f:214`), donc `MTOT` = `MDYN(I)` à cet endroit.
- **Vérification :** `make test` : le gradient analytique de `rel_rvp/univ` est maintenant égal aux différences finies (il s'en écartait de 4 % sur `p1`). Référence inchangée.
- **Effet sur les résultats :** aucun χ² ni MAP ne change. Seul le Levenberg-Marquardt universel avec des vitesses de la planète, une orbite et sans vitesses de l'étoile est concerné : il suivait une direction légèrement fausse, et la matrice de covariance qu'il calcule était légèrement fausse. Le MCMC n'utilise pas ces dérivées.

## M008 — Corrige B3 : écriture hors du tableau `DMUDP` dans `MRQELEMENTS` universel

- **Date :** 2026-11-06
- **Fichier :** `Code/mcmc/sub/kepu.f`, routine `MRQELEMENTS` (un bloc repéré `M008`).
- **Changement :** après le calcul de l'erreur sur la masse, le code remet à zéro les cases `DMUDP(DKAL+1-NEL)` et `DMUDP(DKAL+6-NEL)`, remplies pour l'orbite précédente. Il le faisait aussi pour la première orbite, où ces indices valent -6 et -1. La remise à zéro est maintenant faite seulement si `I.GT.1`. Elle reste nécessaire ensuite, car `DMUDP` sert au calcul de l'erreur sur `tp`.
- **Pourquoi :** en compilation normale, deux zéros étaient écrits hors du tableau, dans une zone qui peut appartenir à une autre variable. `kepellip.f` n'est pas concerné : il ne fait pas cette remise à zéro et ne réutilise pas `DMUDP`.
- **Vérification :** `make test` : 363 OK, 0 échec, référence inchangée. Le contrôle des tableaux (`-fcheck=all`) passe maintenant sur les 11 cas universels, qui s'arrêtaient tous sur cette ligne.
- **Effet sur les résultats :** aucun χ² ni MAP ne change. Les erreurs affichées après le LM universel ne dépendent plus du contenu imprévisible de la mémoire.

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
