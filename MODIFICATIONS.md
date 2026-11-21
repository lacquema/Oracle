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

Ces anomalies existent dans le code de `bcced4b`. `make test` les affiche comme `CONNU`. Une fois corrigée, une anomalie passe dans une entrée `M...` ci-dessous. Corrigées à ce jour : B5 et B6 (M005), B1, B2 et B4 (M006), B8 (M007), B3 (M008), B7 (M009), B13 (M010), B10 (M011), B9 (M012), B11 (M013), B12 (M015), B15 (M016), B18 (M020), B19 (M022), B20 (M025), B21 (M026), B22 (M027), B14 (M029), B24 (M033), B25 (M034), B23 (M035), B26 (M036), B27 (M037), B28 (M038), B29 (M040), B30 (M041), B31 (M042), B32 (M043), B33 (M044), B34 (M045), B35 (M046), B36 (M047), B37 (M048), B38 (M049), B16 et B17 (M050).

| N° | Où | Problème | Effet constaté |
|---|---|---|---|
| B23 | Levenberg-Marquardt (prior de type 3, « Linear ») puis `SEARCH_BETA` (`mcmc.f`) | Le L-M ne respecte pas les bornes d'un prior plat sur une masse et peut finir hors de ces bornes, où le prior est nul (`MAP = −1e30`). Le MCMC part alors d'un point impossible : aucun pas n'est accepté, même un pas nul, les pas (`BETA`) tombent à 0 et `SEARCH_BETA` boucle sans fin, sans message (sa seule sortie d'échec est un pas trop grand, `BETAMAX`). | **Corrigée par M035.** Trouvée le 2026-11-13 : vitesses de l'étoile seules, prior plat sur la masse de l'étoile de 0,1 à 10 Msun : le L-M finit à `M*` = 1e-9 Msun et `a` = 0,016 au, puis le run tourne 8 h sans écrire de ligne (la sortie, mise en tampon, faisait croire à un L-M qui ne finissait pas). Vu avec `sample` sur une copie compilée avec `-g`, et en écrivant `BETA` et le taux d'acceptation dans `SEARCH_BETA` : tous nuls après 1,3 million de balayages. Les cas normaux du L-M convergent en 4 à 17 itérations, GG Tau en 510. |
| B24 | `STORE_SOLUTION` (`io.f`, les deux versions, l. 1320 et 1473) | `P` déclaré `P(NPAR)` mais lu jusqu'à `P(NPAR+2)` (χ² et MAP). | **Corrigée par M033.** Bénin : le tableau passé a bien `NPAR+2` éléments (`PSAV`). Signalé par `-fcheck=all` (campagne du 2026-11-11) ; contourné dans la copie de test. |
| B25 | `MCMCP` (`mcmc.f:97`) et `WRITE_DUMP` (`io.f:637`) | `PROBINFER` est écrit dans `dump.dat` avant d'avoir reçu une valeur (il n'est fixé qu'au premier passage dans `GELMAN_RUBIN`). | **Corrigée par M034.** Bénin : valeur relue puis recalculée à la reprise. Vu avec `-finit-real=snan` (campagne du 2026-11-11). |
| B26 | `INIT_MCMC` (`mcmc.f`), contrôle du départ ajouté par M035 | `MRQCHECK` (paramètres du L-M) appelée avec les paramètres du MCMC, différents en variables universelles (paramètres 6 et 7). | **Corrigée par M036.** Trouvée à l'audit du 2026-11-14 : en universel, run refusé à tort avec un prior plat sur une masse dynamique. Régression de M035, restée un jour, jamais poussée. |
| B27 | `READ_INPUTLINE` (`io.f`), conversion (Sep, PA) → (Dec, RA) | Covariance des erreurs : `Sep` au lieu de `Sep²` devant la variance de l'angle de position. | **Corrigée par M037.** Trouvée à l'audit du 2026-11-14 : coefficient de corrélation faux (exemple : −0,163 au lieu de +0,151), donc χ² faux pour toutes les données en (Sep, PA). Jamais testé auparavant. |
| B28 | `METRO_STEP` (`mcmc.f`), probabilité d'acceptation | `MIN(NaN, 1)` vaut 1 avec gfortran : un pas d'essai dont le χ² est NaN est accepté. | **Corrigée par M038.** Trouvée à l'audit du 2026-11-14. La chaîne garde un χ² NaN et accepte ensuite tout. Aucun NaN dans les runs normaux testés ; vu seulement dans les cas dégénérés refusés depuis par M032. |
| B29 | `SEARCH_BETA` (`mcmc.f`), réglage des pas avant le MCMC | Sa seule sortie d'échec est un pas trop grand (`BETAMAX`). Si aucun pas n'est accepté (départ où le prior est nul ou le χ² NaN), les pas tombent à 0 et la boucle tourne sans fin, sans message. | **Corrigée par M040.** Cause du run de 8 h sans sortie de B23 ; M035 et M036 refusent depuis les départs interdits connus, M040 est le filet de sécurité pour les autres. Trouvée à l'audit du 2026-11-14 (A3). |
| B30 | `WRITE_DUMP` et `READ_DUMP0` (`io.f`), fichier de sauvegarde | Les erreurs d'écriture du dump ne sont pas vérifiées : un dump incomplet (disque plein, quota) remplace le bon. À la reprise, un `dump.dat` absent donne une erreur obscure de gfortran. | **Corrigée par M041.** Trouvée à l'audit du 2026-11-14 (A8). Reproduite le 2026-11-18 sur un petit disque virtuel : le dump de 9,4 Mo remplacé par un fichier de 6,9 Mo, et `Dump file created... ok`. |
| B31 | `CHANGE_DIRS` (`mcmc.f`), nouvelles directions pendant le run | L'échec de `SEARCH_BETA` (un pas au-delà de `BETAMAX`, direction que les données ne contraignent pas) est ignoré : la chaîne garde ce pas. Au départ (`INIT_MCMC`), le même échec est rattrapé par un nouvel essai. | **Corrigée par M042.** Relevée à l'audit du 2026-11-14 (avec A3). Jamais vue dans un run. |
| B32 | `READ_INPUT` (`io.f`), question sur le jitter | `JITNUM` n'est lu que s'il y a des VR de l'étoile (type 2), mais `COMPUTE_MAP` le teste dès qu'il y a des VR (type 2 ou 3). Sans VR de l'étoile, il n'est jamais fixé. | **Corrigée par M043.** Relevée à l'audit du 2026-11-14 (A6). Sans effet constaté : une variable de module part à 0 avec gfortran. Une autre valeur fausserait le MAP (voir M043). |
| B33 | `ENTER_ORBITS`, `CONVERT_UNIT`, lecture de la distance et du jitter (`io.f`) | Unités non vérifiées : une unité de masse autre que `ms`/`mj` (orbite ou prior) ou de distance autre que `pc`/`mas` laisse la valeur indéfinie, sans message. Une réponse sur le jitter autre que 0 ou 1 ajoute un paramètre jamais utilisé. | **Corrigée par M044.** Relevée à l'audit du 2026-11-14 (A7). Constaté le 2026-11-18 : avec `mx` ou un prior en `xx`, le run tourne normalement, sur des masses fausses. |
| B34 | `ENTER_ORBITS` (`io.f`), inclinaison de départ | i = 0 exactement : 0/0 (nœud indéfini) et χ² NaN dès le L-M. i = 180 exactement : `tan(i/2)` infini, le L-M oscille sans finir. | **Corrigée par M045.** Relevée à l'audit (A7). |
| B35 | `KEPU` (`kepu.f`), équation de Kepler universelle non résolue | Le code écrit `epic fail !` puis s'arrête par `stop`, c'est-à-dire avec le code 0, comme un run réussi. | **Corrigée par M046.** Trouvée le 2026-11-18 (départ à i = 180°, universel). |
| B36 | `MCMCP` (`mcmc.f`), taille des groupes de pas | `NGROUP = NGEL`, lu dans le dump à la reprise. Un dump d'avant M016 (`NGEL` = 0) donne une division entière par zéro, `NGEL` > 200 un débordement de `LL(NPAR,200)`. | **Corrigée par M047.** Relevée à l'audit (A4). Constaté le 2026-11-18 : `NGEL` = 0, arrêt en 23 lignes avec le code 0 ; `NGEL` = 300, SIGBUS. |
| B37 | `JD_TO_DATE` (`utils.f`) | Jour julien vers date : le demi-jour de l'algorithme (jours comptés depuis midi) manque, et la fraction de jour est rangée dans un entier. Jusqu'à 1,5 jour d'écart sur l'aller-retour avec `DATE_TO_JD`. | **Corrigée par M048.** Relevée à l'audit. Ne sert qu'à écrire des données simulées au format date (choix 5). |
| B38 | Messages autour du dump (`mcmc.f`, `io.f`) | Après M041, un dump non écrit était suivi de `ok`, et l'arrêt sur `stop.flag` annonçait `dump file written` même sans dump (`oracle run` le croyait). À la reprise, un `logfile.dat` absent arrêtait le run sur une erreur de gfortran. | **Corrigée par M049.** Limites relevées en écrivant M041 (2026-11-18). |
| B39 | Paramétrisation du MCMC (`kepellip.f`, `kepu.f`) avec des données sans astrométrie | Sans astrométrie (ni relative ni Gaia), Ω et i ne sont vus par aucune donnée, mais ils entrent dans les paramètres 2 à 7 (ω + Ω, Ω − ω, i). Ces paramètres ne se mélangent pas et le run n'aboutit pas, avec ou sans loi a priori sur les masses, sans aucun avertissement. M032 ne refuse que le cas sans prior et sans VR relatives, et son message laisse croire qu'un prior suffit. | Relevée le 2026-11-18 (classée d'abord comme piste F4). VR de l'étoile seules + prior : pas d'inférence en 8 h (2026-11-13). Sans prior : VR de l'étoile seules (une et deux orbites) et VR de l'étoile + VR relatives, pas d'inférence en 45 min ; plantage (0/0, NaN) en variables universelles avec les VR seules (M032). Le 2026-11-18 aussi : deux orbites, astrométrie du seul compagnon intérieur + VR de l'étoile, sans prior : l'orbite intérieure converge, **l'orbite extérieure (vue seulement par les VR) non** (r de 1,13 à 1,48 en 90 min). Le défaut porte donc sur **chaque orbite sans astrométrie**, pas seulement sur les runs sans aucune astrométrie. Correction prévue : F4 (paramétrisation propre à ce cas), puis revoir M032. Non testée par `make test`. |
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
| B14 | Tous les appels à `RANDOM_NUMBER` (`mcmc.f`, `MRQRENEW` de `kepellip.f` et `kepu.f`, `io.f`), sans `RANDOM_SEED` | La graine du générateur aléatoire n'est jamais fixée. Avec gfortran, elle change alors à chaque exécution | **Corrigée par M029** (graine tirée puis écrite, choix d'Antoine du 2026-11-11). Deux lancements identiques ne donnaient pas le même résultat : on ne peut pas reproduire un run à l'identique. Le LM en dépend aussi, car `MRQRENEW` repart d'angles tirés au hasard après une solution impossible : sur GG Tau, il finit dans un autre minimum environ une fois sur cent (396 sur 400 avant M010, 398 sur 400 après). Choix à faire avec Hervé : graine fixe, graine donnée en entrée, ou graine écrite dans la sortie pour pouvoir refaire un run. Dans `make test`, l'essai GG Tau est relancé jusqu'à trois fois pour ne pas échouer au hasard. |
| B15 | `MCMCP` (`mcmc.f`), intervalle entre deux contrôles de Gelman-Rubin une fois la convergence atteinte | `NGEL = MIN(NGROUP,NINT(LONG*1d-2/NGROUP)*NGROUP)` vaut 0 tant que `LONG` < 5 000 balayages : aucun pas entre deux contrôles | **Corrigée par M016.** Quand la convergence arrivait avant 5 000 balayages, les chaînes ne bougeaient plus et les mêmes 8 positions étaient gardées à chaque contrôle jusqu'à remplir la réserve. Constaté le 2026-11-08 sur `rel_vr_planete` (`--rapide`) : 9 solutions distinctes sur 20 001, ω hors de l'intervalle à 99,7 %. Explique aussi une partie de B12. Les longs runs réels (GG Tau, 51 Eri), qui convergent après des millions de pas, ne sont pas touchés. |
| B16 | `MRQELEMENTS` (`kepellip.f`, `kepu.f`), erreur sur l'inclinaison | `PLA(I)%DINC` n'est calculée que dans la branche `RADVEL` : sans vitesses radiales, `DIDP` est rempli mais `ERR_COMB` n'est jamais appelée, et l'erreur affichée sur i vaut 0 | **Corrigée par M050.** Constaté le 2026-11-09 dans le dépôt Gaia (BH3 en Levenberg-Marquardt seul). Sans effet sur le MCMC : seul l'affichage du Levenberg-Marquardt (choix 4) est touché. |
| B17 | `MRQELEMENTS` (`kepellip.f`, `kepu.f`), deux branches | `DIDP(4)` et `DIDP(5)` au lieu de `DIDP(DKAL+4)` et `DIDP(DKAL+5)` : pour la deuxième orbite et les suivantes, l'erreur sur i est calculée avec les paramètres de la première | **Corrigée par M050.** Constaté le 2026-11-09 en relisant B16. Même portée : affichage du Levenberg-Marquardt seulement. |
| B18 | Programmes principaux (`CHARACTER*80 FILES`), `WRITE_DISTRIB` et `WRITE_DISTRIB_DOUBLING` (`CHARACTER*80 LINE`), `READ_DATAFILE` (`STATUS='UNKNOWN'`) | Les noms de fichiers sont coupés à 80 caractères. Le fichier de données au nom coupé est ouvert en `UNKNOWN`, donc créé vide s'il n'existe pas | **Corrigée par M020.** Avec un chemin de plus de 80 caractères (l'interface écrit des chemins absolus), le programme s'arrête sur `End of file` à la lecture des données, sans dire pourquoi. Le `gzip` des résultats avait la même limite. Constaté le 2026-11-10 en testant l'interface Gaia ; déjà connu d'Antoine. |
| B19 | `MRQRENEW` (`kepellip.f`, `kepu.f`), redémarrage du L-M après une solution impossible | Pendant le premier passage sans jitter des vitesses (`MRQMIN(NPAR-1,...)`), la masse de l'étoile est en `P(NPAR-1)` = `P(NEL*NPLA+2)`. `MRQRENEW` y écrit ln(jitter) et met la masse en `P(NPAR)`, hors du tableau du passage (`P` déclaré `DIMENSION(NPAR)`) | **Corrigée par M022.** Vitesses de l'étoile avec jitter et prior fixe sur sa masse : après un redémarrage, la masse de l'étoile vaut exp(ln σ_jit) (0,0195 Msun pour un jitter initial de 0,01 km/s) ; le prior fixe donne alors une pénalité de 10¹² que le L-M ne rattrape plus, et la solution finale est fausse. Constaté le 2026-11-10 sur HD 114762 (vitesses de Lick seules, χ² = 2,3·10⁶ au lieu de 119). Sans jitter, pas d'effet. |
| B20 | `MRQRADVEL` (`mrqfit.f`), matrice `ALPHA` du L-M | Pour ln σ_jit, vraie dérivée seconde du χ² (négative quand les résidus dépassent σ) et termes croisés jitter × orbite, alors que les autres paramètres sont en Gauss-Newton | **Corrigée par M025.** Le L-M s'arrête avant le minimum, le jitter glissant vers 0 : sur les vitesses de HD 114762, χ² + Fmap = 89,75 (σ_jit = 1,2 m/s) au lieu de 80,05 (13,2 m/s). Constaté le 2026-11-11. |
| B21 | `MRQSTEP` (`mrqfit.f`), pas d'essai hors des bornes | Chaque pas hors des bornes ou NaN appelle `MRQRENEW`, qui repart de l'orbite d'entrée avec des angles tirés au hasard, au lieu d'augmenter λ | **Corrigée par M026.** Boucle sans fin si le premier pas depuis l'orbite d'entrée sort des bornes : 1,7 million de redémarrages en 60 s sur HD 114762 (départ à i = 14,8°). Constaté le 2026-11-10 et analysé le 2026-11-11. |
| B22 | `ENTER_ORBITS` (`io.f`) puis toutes les divisions par `S` = e/√(1−e²) | Excentricité de départ e = 0 : 0/0 | **Corrigée par M027.** Arrêt dans `KEPLER0` (`glup`) dès le début du L-M. Constaté le 2026-11-11. |

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
| F4 | **Corrigera B39 (anomalie, pas seulement une amélioration : Antoine, 2026-11-18).** Ajustement des vitesses radiales seules (et de Gaia seule), pour que le code soit général. Aujourd'hui refusé sans prior sur les masses (M032) ; avec un prior, le MCMC ne converge pas en 8 h. | Essai du 2026-11-13 (VR de l'étoile seules, prior gaussien sur la masse de l'étoile, 61 millions de pas) : convergent la période, `V0` et la masse (R ≤ 1,02) ; pas les paramètres 2 à 7 (R de 1,11 à 1,49), qui contiennent `Ω` ou `i` (`ϖ = ω + Ω` dans les paramètres 2 et 3, `Ω − ω` dans 6 et 7), alors que les VR ne mesurent pas `Ω` ni `i` séparément de la masse. Piste : une paramétrisation propre à ce cas (`P, K, e, ω, Tp`, comme les codes de VR), avec `i` et `Ω` hors de l'ajustement. Antoine (2026-11-13) : piste d'amélioration pour plus tard. **Complément du 2026-11-18** (systèmes synthétiques sans prior, version classique, 45 min, `Oracle-Gaia/Tests/synthetic.py`) : VR de l'étoile seules, une et deux orbites : pas de plantage ni de NaN, mais pas d'inférence (r jusqu'à 1,44) ; **VR de l'étoile + VR relatives, sans astrométrie** (cas accepté par M032 et par le site) : même chose, r de 1,12 à 1,52 sur les paramètres qui contiennent `Ω` ou `i` (Ω n'est vu par aucune donnée sans astrométrie). Avec de l'astrométrie relative (deux orbites, des deux compagnons ou du seul extérieur, avec VR de l'étoile), le run converge sans prior et retrouve toutes les masses. Diagnostic d'Antoine (2026-11-18) : le problème vient de la paramétrisation, pas utilisable sans astrométrie ; le refus (M032, Prepare, `translate.py`) n'est qu'un garde-fou en attendant, on n'en ajoute pas. **À vérifier quand F4 sera fait :** relancer ces cas (`vr_seules_libre`, `deux_vr_seules_libre`, `vr_vrel_libre`), avec et sans prior ; revoir alors la condition de refus (M032, site, `translate.py`) et la passer si possible en avertissement (principe d'Antoine : on prévient, l'utilisateur décide) ; garder le prior seulement là où les masses sont vraiment dégénérées (fonction de masse). **Incohérence de M032 notée le 2026-11-18 (non corrigée, choix d'Antoine) :** M032 mêle deux problèmes. (1) Dégénérescence des masses sans données relatives : une loi a priori la règle (Gaia seule + prior converge, BH3). (2) Ω et i invisibles sans aucune astrométrie : rien ne converge, avec ou sans prior (F4). D'où deux incohérences : VR de l'étoile seules + prior est accepté et ne converge pas, alors que le message de M032 dit d'ajouter un prior ; VR de l'étoile + VR relatives sans astrométrie est accepté et ne converge pas. Règle cohérente à reprendre avec F4 : sans aucune astrométrie (ni relative ni Gaia), prévenir que la paramétrisation actuelle ne convergera pas ; sans données relatives et sans prior, demander une loi a priori sur une masse. |
| F5 | Excentricité légèrement surestimée (biais de Lucy et Sweeney 1971) : vérifier quel prior sur `e` le code impose réellement (jacobien des paramètres 2 et 3, `e (cos ϖ, sin ϖ)/√(1−e²)`) et le comparer à un prior uniforme sur `e` ; confirmer le biais avec des graines indépendantes entre les cas ; décider s'il faut changer le prior ou seulement l'expliquer. **Expliqué dans la documentation du site (section « Reading the results ») le 2026-11-13.** | Campagne de tests du 2026-11-11 au 13 (bilan dans le dépôt privé, `Tests/campagne/README.md`) : sur 96 runs à une orbite (vraie `e` = 0,3), quantile moyen de la vraie `e` 0,43 ± 0,03 ; elle tombe sous l'intervalle à 68 % dans 32 % des runs et au-dessus dans 8 %. Environ 2 σ. Le bruit allonge plus souvent qu'il ne raccourcit le vecteur `(e cos ω, e sin ω)`. Demandé par Antoine (2026-11-13). |

## Modifications

## M059 — Corrige C5 et C6 : données simulées (choix 5) relisibles, et limites contrôlées

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`SIMUL_DATA`, un bloc ; `WRITE_ASTROM`, deux blocs).
- **Changement :** `WRITE_ASTROM` écrit en tête ce que lit `READ_DATAFILE` : le décalage JD, une ligne de tirets, le numéro du type (1), puis l'en-tête de colonnes, les tirets et les mesures (le décalage n'est plus écrit après l'en-tête). `SIMUL_DATA` s'arrête avec un message au-delà de `NMW` (100) dates, ou pour un indice de compagnon hors de 1..`NPLA`.
- **Pourquoi :** constats C5 et C6 de l'audit du 2026-11-18. Le fichier produit par le choix 5 ne se relisait pas (« Bad real number in item 1 of list input ») ; au-delà de 100 dates, `IPLW` et `JDW` étaient écrits hors de leurs bornes.
- **Vérification (2026-11-18) :** choix 5 (une orbite, 4 dates) puis choix 4 sur le fichier produit : avant, code 2 ; après, lecture et L-M normaux (code 0). 101 dates : « Error : more than 100 dates », code 1. `make test` passe.
- **Effet sur les résultats :** aucun (outil de simulation).

## M058 — Corrige C10 : point de départ de l'équation de Kepler universelle pour une orbite hyperbolique

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/kepu.f` (`kepu_p3solve`, un bloc).
- **Changement :** dénominateur de la cubique `(mu-alpha*q)/6` au lieu de `-alpha*q/6`.
- **Pourquoi :** constat C10 de l'audit du 2026-11-18. Au périastre, f(s) = q·s·c1 + μ·s³·c3 − dt ≈ q·s + (μ − αq)·s³/6 ; c'est la cubique de Swift (`drift_kepu_p3solve`, avec r0 = q), et le développement qu'utilisent déjà `fpp` et `fppp` de `kepu_new`. Sans μ, le point de départ s'écartait de la racine d'un facteur 2 à 4 dès que dt n'est pas petit.
- **Vérification (2026-11-18) :** calcul indépendant (q = 2, μ = 1) : vraie racine, ancien départ, nouveau départ. α = −0,01, dt = 1000 : 17,01 ; 63,96 ; 17,83. α = −0,1, dt = 100 : 6,97 ; 13,04 ; 7,52. α = −1, dt = 1000 : 6,51 ; 14,28 ; 12,49. `make test` passe (χ² universel = classique = Python inchangés : seul le point de départ de Newton change).
- **Effet sur les résultats :** orbites hyperboliques (e > 1, variables universelles) seulement : moins d'itérations et moins d'échecs de la résolution ; la solution elle-même ne change pas.

## M057 — Corrige C4 : indice de compagnon contrôlé à la lecture des données

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`READ_DATAFILE`, premier passage, un bloc).
- **Changement :** pour les types 1 et 3, un indice de compagnon hors de 1..`NPLA` arrête le programme (`STOP 1`) avec « Error : companion index k in the data file, but n orbit(s) ».
- **Pourquoi :** constat C4 de l'audit du 2026-11-18 : l'indice lu servait directement à `PLA(IPLA)`, d'où une écriture hors du tableau. Prepare refusait déjà ce cas ; le Fortran utilisé directement plantait.
- **Vérification (2026-11-18) :** cas `rel` de `make test` (une orbite) avec une ligne d'indice 2 : avant, « Segmentation fault », code 139 ; après, le message ci-dessus, code 1. `make test` passe.
- **Effet sur les résultats :** aucun.

## M056 — Corrige C2 : trois arrêts sur erreur sortaient avec le code 0

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/sub/kepellip.f` (`KEPLER0`) ; `Code/mcmc/sub/utils.f` (`JACOBI`) ; `Code/mcmc/sub/io.f` (`READ_DATAFILE`, au milieu du bloc M020, fermé puis rouvert : les blocs imbriqués sont interdits).
- **Changement :** `STOP 1` et un message commençant par « Error : » quand l'équation de Kepler classique n'est pas résolue (« glup »), quand `JACOBI` ne converge pas, et quand le fichier de données manque.
- **Pourquoi :** constat C2 de l'audit du 2026-11-18. gfortran rend 0 pour `STOP` et pour `STOP 'texte'` (vérifié) : `oracle run` annonçait « done » pour « glup » et pour l'échec de `JACOBI` (le fichier absent était vu grâce à « not found »). M046 avait corrigé le même défaut en universel seulement.
- **Vérification (2026-11-18) :** fichier de données absent : code 1, « Error : data file not found : absent.dat » (avant : code 0). Les deux autres cas (anomalie moyenne NaN, covariance avec NaN) ne se provoquent pas sans modifier le code : changements relus. `make test` passe.
- **Effet sur les résultats :** aucun.

## M055 — Corrige C1 : départ du MCMC, plus de boucle sans fin quand un pas grandit sans limite

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`INIT_MCMC`, trois blocs).
- **Changement :** les échecs de `SEARCH_BETA` au départ du MCMC (un pas au-delà de `BETAMAX`) sont comptés ; après 20, arrêt (`STOP 1`) avec un message qui donne les numéros des paramètres dont le pas a dépassé `BETAMAX`, c'est-à-dire que ni les données ni les lois a priori ne limitent.
- **Pourquoi :** constat C1 de l'audit du 2026-11-18 (déjà vu le 14 : A2). `INIT_MCMC` recommençait la recherche sans limite (« Beta initiation process failed. Retrying... ») : le run ne finissait jamais, ne disait pas pourquoi et remplissait le journal. M040 traitait les pas qui tendent vers 0, M042 l'échec dans `CHANGE_DIRS`, pas ce cas.
- **Vérification (2026-11-18) :** cas `rv_seul` de `make test` avec une ligne de tirets en trop (aucune donnée lue, avertissement M051), choix 3 : avant, 7 millions de lignes en 45 minutes sans fin (runs synthétiques du jour) ; après, arrêt en 0,2 s, code 1, après 20 essais, paramètre 8 (V0, que rien ne contraint sans données) nommé. `make test` passe.
- **Effet sur les résultats :** aucun pour un run qui démarre ; un run qui ne pouvait pas démarrer s'arrête avec une explication.

## M054 — Corrige C9 : dérivées secondes du terme de loi a priori dans le L-M classique ; nouveau contrôle dans `make test`

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/sub/kepellip.f` (`MRQMAPASTROM`, deux blocs) ; `Tests/driver.f` et `Tests/run_tests.py` (nouveau contrôle).
- **Changement :**
  1. `D2FMAP(IM+7,IM+7)` part de sa propre valeur, et non de `D2FMAP(IM+6,IM+6)` (déjà augmentée à la ligne précédente : les deux termes s'annulaient et (7,7) recevait l'ancienne valeur de (6,6)).
  2. Termes croisés entre deux orbites i ≠ k, `D2LMUF(IM+1,KM+6)` et `(IM+1,KM+7)` : facteur `FF2` = FF1_i·FF1_k ajouté, comme dans le terme symétrique.
  3. `make test` compare désormais, pour chaque cas et chaque version, les dérivées secondes de `MRQMAP` aux différences finies centrées de son gradient (écart relatif par colonne au plus 1e-4, contrôle « dérivées secondes du prior »).
- **Pourquoi :** constat C9 de l'audit du 2026-11-18. La dérivée seconde de la loi a priori en P7 était perdue pour chaque orbite, dans tout run classique qui passe par le L-M, et la matrice était asymétrique avec une loi a priori sur les masses de plusieurs orbites. `make test` ne contrôlait que le gradient. La version universelle (`kepu.f`) était juste (vérifiée numériquement).
- **Vérification (2026-11-18) :** programme de contrôle `Oracle-Gaia/Tests/campagne/audit_2026-11-18/d2_derivees_secondes.f`, cas `rel_2orb` avec loi a priori sur la masse totale et sur m₂ : avant, (7,7), (14,14), (8,6), (8,7), (1,13), (1,14) faux ; après, aucun écart. Nouveau contrôle de `make test` : passe sur tous les cas après M054 ; essayé sur une copie sans M054 (voir plus bas), il échoue.
- **Effet sur les résultats :** L-M classique seulement (chemin, nombre d'itérations, covariance finale et erreurs affichées sur P et t_p). Le MCMC n'utilise pas ces dérivées.

## M053 — Corrige C11 : le L-M respecte aussi la borne inférieure de l'excentricité

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `kepu.f` (`MRQCHECK`, un bloc chacun).
- **Changement :** un pas du L-M qui donne `e` sous `VPRIOR(3)%BOUND(1)` est refusé, comme un pas au-dessus de la borne supérieure.
- **Pourquoi :** constat C11 de l'audit du 2026-11-18. `MRQCHECK` ne testait que la borne supérieure, alors que le MCMC (`RATIO_PRIORS`) refuse les deux, et que le départ du MCMC est contrôlé par `RATIO_PRIORS` depuis M036 : avec une borne inférieure non nulle, le L-M pouvait finir dessous et le run s'arrêtait (« starting point of the MCMC is outside the bounds »).
- **Vérification (2026-11-18) :** cas `rel` de `make test` (e vrai 0,3), e borné à [0,5 ; 0,9], départ e = 0,6, choix 2. Avant : L-M fini à e = 0,283, arrêt, code 1. Après, classique et universel : L-M fini dans les bornes (e = 0,9, la vraie valeur étant hors de l'intervalle imposé), MCMC lancé. `make test` passe.
- **Effet sur les résultats :** aucun quand la borne inférieure vaut 0 (cas habituel).

## M052 — Corrige C8 : le jitter initial des VR donné par l'utilisateur est enfin utilisé ; jitter nul remplacé

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/main/astrom_mcmco.f` et `astrom_univ_mcmco.f` (quatre blocs chacun) ; `Code/mcmc/sub/kepellip.f` et `kepu.f` (`INIT_DATA`, un bloc chacun).
- **Changement :** le programme principal garde le jitter initial (`SJV0 = STAR%SIGJV`) avant le premier passage du L-M et repart de `LOG(SJV0)` au second passage, au lieu de `LOG(STAR%SIGJV)`. Au choix 3 (MCMC sans L-M), il remet aussi `P` dans la disposition complète avant `MRQELEMENTS` (`P(NPAR) = P(NPAR-1)`, `P(NEL*NPLA+2) = LOG(SJV0)`). Dans `INIT_DATA`, un jitter initial nul ou négatif est remplacé par `SIGV0` (0,1 km/s), avec un message.
- **Pourquoi :** constat C8 de l'audit du 2026-11-18 (`Oracle-Gaia/Notes/audit_2026-11-18.md`). Pendant le premier passage (sans jitter), la masse de l'étoile occupe la case du jitter ; `MRQELEMENTS`, qui lit toujours le jitter en `P(NEL*NPLA+2)`, prenait la masse pour le jitter, et le second passage (ou le MCMC au choix 3) partait d'un jitter égal à m₀ en unités internes (0,512 km/s par masse solaire), quelle que soit la valeur donnée. Avec Gaia (version privée), la case contient le premier paramètre Gaia, et le L-M finissait à un jitter de milliers de km/s (correction propre à Gaia : G016 dans `Oracle-Gaia`). Le plancher corrige A10 de l'audit du 14 : un jitter donné à 0 donnait `LOG(0)` une fois C8 corrigé.
- **Vérification (2026-11-18) :** cas `rel_rv_jit` de `make test`, classique et universel. Choix 3 : jitter de départ 0,0200 km/s pour 0,02 donné (avant : 0,5226, soit m₀ = 1,02 Msun), 0,1 km/s pour 0 donné (message affiché). Choix 2 : le L-M trouve 0,00803 km/s dans les quatre cas. `make test` passe.
- **Effet sur les résultats :** point de départ du second passage du L-M et du MCMC au choix 3. Le minimum du L-M peut changer quand l'ancien départ l'égarait (cas Gaia surtout).

## M051 — Avertissement : type de données annoncé sans aucune mesure lue

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`READ_DATAFILE`, après le premier passage de lecture, repère `[MOD M051]`).
- **Changement :** après le premier passage, pour chaque type de données annoncé comme présent (`ISDATA`), si aucune mesure de ce type n'a été lue, le programme l'écrit (« Warning : no data of type k read in FICHIER, although this type was declared present », puis le rappel du format). Le run continue : c'est un avertissement, pas un arrêt (Antoine : on prévient, l'utilisateur décide).
- **Pourquoi :** la lecture s'arrête sans message dès que la ligne attendue comme numéro de type n'est pas un nombre, par exemple deux lignes de tirets à la suite après l'en-tête. Le run tournait alors sans ces données : χ² initial nul, et un MCMC qui ne démarre jamais (« Beta initiation process failed » en boucle). Trouvé le 2026-11-18 sur des systèmes synthétiques (`Tests/synthetic.py` de `Oracle-Gaia`, vitesses sans astrométrie relative : une ligne de tirets en trop), après 45 minutes de calcul perdues. La règle de format est notée dans `CLAUDE.md` (« Format du fichier de données »).
- **Vérification (2026-11-18) :** `make test` passe ; fichier de vitesses de l'étoile précédé de deux lignes de tirets : l'avertissement s'affiche pour le type 2 ; même fichier corrigé : aucun avertissement.
- **Effet sur les résultats :** aucun.

## M050 — Corrige B16 et B17 : erreur sur l'inclinaison affichée par le Levenberg-Marquardt

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `kepu.f` (`MRQELEMENTS`, deux blocs chacun) ; `Tests/driver.f` (écrit i et son erreur) ; `Tests/run_tests.py` (nouveau contrôle, lecture des nombres à exposant à trois chiffres).
- **Changement :** dans les deux branches (avec et sans VR), les dérivées de i sont rangées en `DIDP(DKAL+4)` et `DIDP(DKAL+5)`, les paramètres de l'orbite traitée, au lieu de `DIDP(4)` et `DIDP(5)`. Dans la branche sans VR, `PLA(I)%DINC` est maintenant calculée (`ERR_COMB`), comme dans l'autre.
- **Pourquoi :** B16 (sans VR, erreur jamais calculée : 0 ou valeur indéfinie, jusqu'à 3e-314 dans le programme de test) et B17 (deuxième orbite et suivantes : erreur calculée avec les paramètres de la première). Les formules des dérivées étaient justes : i = 2 atan(T), T = |(P4, P5)|, di/dPk = 2 Pk / (T (1 + T²)).
- **Vérification (2026-11-18) :**
  - Nouveau contrôle de `make test`, « inclinaison et son erreur » : pour chaque cas, chaque version et chaque orbite, i et son erreur sont recalculés en Python à partir des paramètres et de la covariance diagonale du programme de test. Avant : 10 échecs sur 24 (B16 : `rel`, `rel_2orb`, `rel_corr`, `rel_seppa` ; B17 : orbite 2 de `rel_rv_2orb_jit`, 0,358° au lieu de 0,914°). Après : 24 sur 24, à 1e-9 près.
  - Affichage réel du L-M (choix 4), classique : `rel`, `+/- 0.000000` devient `+/- 0.765222` ; `rel_2orb`, `0.000000` et `0.000000` deviennent `0.787620` et `0.520691` ; `rel_rv_2orb_jit`, orbite 1 inchangée (0,664372), orbite 2 de 0,638223 à 0,514802. Inclinaisons inchangées.
  - `make test` passe (394 OK).
- **Effet sur les résultats :** aucun sur le MCMC ni sur l'ajustement ; seule l'erreur sur i affichée par le L-M change.

## M049 — Corrige B38 : messages du dump fidèles, `logfile.dat` recréé à la reprise

- **Date :** 2026-11-18
- **Fichiers :** `Code/mcmc/sub/io.f` (`WRITE_DUMP` : argument de sortie `OKD`, trois blocs ; `INIT_SIMU` : ouverture du journal) ; `Code/mcmc/sub/mcmc.f` (déclaration de `OKD` dans `MCMC` et `MCMCP`, les trois appels de `WRITE_DUMP`).
- **Changement :**
  - `WRITE_DUMP` renvoie `OKD`, vrai seulement si le dump a été écrit et renommé. Les appelants n'écrivent `ok` que dans ce cas.
  - Arrêt sur `stop.flag` (M039) : si le dump n'a pas pu être écrit, `Stopped on request, dump file could not be written : previous dump kept` et `STOP 1`. `oracle run` cherche la phrase `Stopped on request, dump file written`, absente ici : il annonce un arrêt sans sauvegarde et la reprise depuis le dump précédent, ce qui est exact. Rien à changer côté Python.
  - Reprise : `logfile.dat` ouvert en `STATUS='UNKNOWN'` (créé s'il manque, à la suite sinon) au lieu de `'OLD'`.
- **Pourquoi :** B38, limites laissées par M041.
- **Vérification (2026-11-18) :** cas `rel` classique, 1 thread, graine 1113 :
  - disque virtuel plein (comme M041) sans `stop.flag` : `Dump file created... Warning : dump file could not be written (disk full ?), previous dump kept`, sans `ok` ; dump d'avant intact ; le run continue.
  - disque plein avec `stop.flag` : avertissements, puis `Stopped on request, dump file could not be written : previous dump kept`, code 1 ; dump d'avant intact.
  - disque normal avec `stop.flag` : `ok`, `Stopped on request, dump file written`, code 0, comme avant.
  - reprise sans `logfile.dat` : `Here we go again`, le run continue, `logfile.dat` créé.
  - `make test` passe.
- **Effet sur les résultats :** aucun.

## M048 — Corrige B37 : conversion jour julien vers date

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/utils.f` (`JD_TO_DATE`, un bloc).
- **Changement :** `I = INT(JD+0.5)` et `F = JD+0.5-I`, `F` réel. Le reste de l'algorithme (Duffett-Smith et Zwart, 2011 ; le même que Meeus) est inchangé.
- **Pourquoi :** B37. L'algorithme compte les jours depuis midi (JD + 0,5) ; le code prenait `INT(JD)` et rangeait la fraction dans un entier. `JD_TO_DATE` ne sert qu'à écrire des données simulées au format jour/mois/année (choix 5, format `f5.2` pour le jour) : les dates écrites étaient décalées de 0,5 à 1,5 jour par rapport aux instants du modèle.
- **Vérification (2026-11-18) :** petit programme qui appelle les deux routines du code, extraites telles quelles.
  - Avant : 2451544,5 (1er janvier 2000 à 0 h) donne le 31 décembre 1999, puis 2451543,5 au retour ; écart maximal de 1,5 jour sur 200 000 dates tirées entre 1582 et 2100.
  - Après : 2451544,5 donne le 1,000 janvier 2000 ; 2451545,25 donne le 1,750 ; écart maximal sur l'aller-retour nul (0,0 j) pour les 200 000 dates.
  - `make test` passe.
  - De bout en bout (ajouté avec M049) : choix 5, une orbite, format jour/mois/année, cinq dates données avec des fractions de jour. Avant : `15.25 6 2010` écrit `14.00 6 2010`, `1.00 1 2000` écrit `31.00 12 1999`, `1.00 3 2024` écrit `29.00 2 2024`. Après : les cinq dates ressortent telles qu'entrées ; positions simulées identiques (même graine).
- **Effet sur les résultats :** aucun sur les ajustements. Les données simulées au format date sont maintenant aux bons instants.

## M047 — Corrige B36 : taille des groupes de pas fixe, indépendante du dump

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`MCMCP`, un bloc).
- **Changement :** `NGROUP = 100` au lieu de `NGROUP = NGEL`, et `NGEL` ramené à 100 s'il vaut autre chose après `INIT_MCMC`.
- **Pourquoi :** B36. Le code écrit toujours `NGEL` = 100 (départ et contrôles de `GELMAN_RUBIN`, M016), mais à la reprise `NGEL` vient du dump.
- **Vérification (2026-11-18) :** reprise (choix 1) pendant 30 s d'un dump du cas `rel` classique, `NGEL` mis à 0, 100 et 300 dans le dump :
  - avant : 0, arrêt après 23 lignes avec le code 0 (`IEEE_INVALID_FLAG`), sans rien calculer ; 300, plantage (SIGBUS) ; 100, normal ;
  - après : les trois reprennent, sorties identiques ligne à ligne (1 346 lignes) ; 100 identique à l'ancien code.
  - `make test` passe.
- **Effet sur les résultats :** aucun pour un dump écrit par le code actuel.

## M046 — Corrige B35 : échec de l'équation de Kepler universelle signalé comme une erreur

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/kepu.f` (`KEPU`, un bloc).
- **Changement :** après `epic fail !` et ses valeurs, le code écrit `Error : universal Kepler equation not solved (dt, q, mu, alpha above : NaN or extreme values). Check the starting orbit.` et s'arrête par `STOP 1` au lieu de `stop`.
- **Pourquoi :** B35. Le code de sortie 0 faisait passer l'échec pour un run réussi (`oracle run`, scripts). Vu avec un départ à i = 180° en universel : le L-M finit dans un minimum absurde (masse de l'étoile 1e-9 Msun), puis la masse NaN arrive dans `KEPU` au départ du MCMC. La cause en amont (départ très loin de la vraie orbite) n'est pas traitée ici.
- **Vérification (2026-11-18) :** même entrée (cas `rel_rv` universel, i de départ 180°, ramenée à 179,9° par M045) : avant, `epic fail !` et code 0 ; après, `epic fail !`, le message d'erreur et code 1. `make test` passe.
- **Effet sur les résultats :** aucun ; seul le code de sortie change.

## M045 — Corrige B34 : inclinaison de départ ramenée dans [0,1° ; 179,9°]

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`ENTER_ORBITS`, un bloc après celui de M027).
- **Changement :** une inclinaison de départ inférieure à 0,1° est remplacée par 0,1°, une inclinaison supérieure à 179,9° par 179,9°, avec un message, comme M027 pour e = 0.
- **Pourquoi :** B34. i = 0 : `tan(i/2)` = 0 et le nœud est indéfini (`CO`, `SO`, `CP`, `SP` en 0/0 dans `ELEMENTS`) : χ² NaN, arrêt par M031 sur un message qui n'indique pas la cause. i = 180 : `tan(i/2)` ≈ 1,6e16, le L-M oscille sans finir (plus de 5 000 lignes en 30 s).
- **Vérification (2026-11-18) :** cas `rel_rv` (vraie i = 62°), classique et universel, 30 s, départs à différentes inclinaisons, code M040 :
  - 0 : χ² NaN dans le L-M (code 1) ; 180 : L-M sans fin ;
  - 0,1°, 1°, 150°, 170°, 175° : bon minimum (χ² 48,76), dans les deux versions ;
  - 0,001°, 0,01°, 0,05°, 120°, 178°, 179°, 179,5°, 179,9° : mauvais minimum dans l'une ou les deux versions (χ² de 5 800 à 23 000), parfois `epic fail` en universel (voir M046). Le L-M dépend du point de départ ; ce n'est pas propre aux inclinaisons extrêmes.
  - Après M045 : départ à 0 → 0,1°, bon minimum (χ² 48,76) dans les deux versions ; départ à 180 → 179,9°, le L-M finit (mauvais minimum, comme 179,9° avant). Départ à 62° : inchangé.
  - `make test` passe.
- **Effet sur les résultats :** aucun pour une inclinaison de départ entre 0,1° et 179,9°.

## M044 — Corrige B33 : unités et réponse sur le jitter vérifiées à la lecture

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`ENTER_ORBITS`, `CONVERT_UNIT`, lecture de la distance et du jitter : quatre blocs).
- **Changement :** arrêt avec un message (`STOP 1`) si l'unité de masse d'une orbite ou d'un prior n'est ni `ms` ni `mj`, si l'unité de distance n'est ni `pc` ni `mas` (majuscules ou minuscules), ou si la réponse sur le jitter n'est ni 0 ni 1. Dans `CONVERT_UNIT`, `UPCASE` porte maintenant sur l'argument `UNT` et non sur `UNIT`, variable de l'appelant passée comme `UNT` : même effet.
- **Pourquoi :** B33. Le lanceur `oracle` du dépôt privé écrit toujours `ms`, `mj`, `pc` ou `mas` ; une entrée écrite à la main pouvait se tromper sans message.
- **Vérification (2026-11-18) :** cas `rel_rv` classique (et `rv_seul` pour le prior), 30 s, code M040 contre M044 :
  - masse en `mx` : avant, le run tourne normalement (masse fausse) ; après, `Error : unknown mass unit for planet # 1 : MX (ms or mj)`, code 1.
  - prior en `xx` : avant, le run tourne ; après, `Error : unknown mass unit in a prior : XX (ms or mj)`, code 1.
  - distance en `kpc` : avant, χ² NaN dans le L-M (code 1, M031) ; après, `Error : unknown distance unit : KPC (pc or mas)`, code 1.
  - jitter `2` : avant, rien en 30 s ; après, `Error : jitter answer must be 1 or 0, not 2`, code 1.
  - `MJ` et `PC` en majuscules : acceptés, comme avant.
  - `make test` passe.
- **Effet sur les résultats :** aucun pour une entrée valide.

## M043 — Corrige B32 : `JITNUM` fixé à 0 avant la question sur le jitter

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`READ_INPUT`, un bloc).
- **Changement :** `JITNUM = 0` juste avant la question sur le jitter, posée seulement s'il y a des VR de l'étoile.
- **Pourquoi :** B32. Avec des VR relatives seules (type 3), `COMPUTE_MAP` teste `RADVEL.AND.(JITNUM.EQ.1)` sans que `JITNUM` ait été lu. Il vaut 0 parce que gfortran met les variables de module à 0, ce que la norme ne garantit pas.
- **Vérification (2026-11-18) :** cas `rel_rvp` (astrométrie relative et VR relatives, pas de VR de l'étoile), classique, 40 s, 1 thread, graine 1113. Copies où `JITNUM` vaut 1 au départ (déclaration `JITNUM = 1` dans `commons.f`), pour simuler une valeur non nulle :
  - sans M043 : MAP de départ faux (−28,796846 au lieu de −28,796788 : terme de prior du jitter compté à tort), puis toute la suite diffère ;
  - avec M043 : sortie identique ligne à ligne (1 417 lignes) à celle du code normal ;
  - code normal avec M043 : identique au code M040.
  - `make test` passe.
- **Effet sur les résultats :** aucun avec gfortran.

## M042 — Corrige B31 : échec de la recherche des pas après un changement de directions

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`CHANGE_DIRS`, une déclaration et un bloc).
- **Changement :** pour chaque chaîne, les pas d'avant sont gardés ; si `SEARCH_BETA` échoue, elle est relancée depuis ces pas, au plus 10 fois en tout. Si elle échoue encore, la chaîne garde ses pas d'avant et le code écrit `Warning : step search failed for chain`. Sans échec, le calcul est le même qu'avant (une seule recherche, partie des mêmes pas).
- **Pourquoi :** B31. Au départ, `INIT_MCMC` relance la recherche tant qu'elle échoue ; `CHANGE_DIRS` ignorait l'échec, et la chaîne gardait un pas supérieur à 1e10 jusqu'au changement de directions suivant. La limite de 10 essais évite de boucler sans fin sur une direction qui ne se contraint pas.
- **Vérification (2026-11-18) :** copie hors du dépôt, `NINFER` = 5000, 1 thread, graine 1113.
  - Cas normaux (`rel`, `rel_rv`, `rel_2orb`, `rel_seppa` en classique et universel, `rel_rv_jit_parab` en universel), 90 s, 2 à 4 changements de directions chacun : sortie identique ligne à ligne à celle du code M040.
  - Échec forcé (copie où `FAIL` est mis à vrai après chaque recherche de `CHANGE_DIRS`), cas `rel`, 120 s : 4 changements de directions, 32 avertissements (8 chaînes × 4), le run continue (`countok` augmente jusqu’à la fin des 120 s).
  - Échec forcé aux deux premiers essais seulement : aucun avertissement, le troisième essai réussit, le run continue.
  - `make test` passe.
- **Effet sur les résultats :** aucun quand la recherche réussit, c'est-à-dire dans tous les runs vus jusqu'ici.

## M041 — Corrige B30 : un dump incomplet ne remplace plus le précédent ; dump manquant signalé à la reprise

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/io.f` (`WRITE_DUMP` : une déclaration et deux blocs ; `READ_DUMP0` : un bloc).
- **Changement :**
  - `WRITE_DUMP` écrit une dernière ligne témoin, `END OF DUMP`, ferme le fichier temporaire, puis relit ses 12 derniers octets (accès direct, sans relire tout le fichier). Si la ligne n'y est pas, le fichier temporaire est effacé, le dump précédent est gardé, le code écrit `Warning : dump file could not be written (disk full ?), previous dump kept` et le run continue. Sinon, renommage comme avant (M018).
  - Si le fichier temporaire ne peut pas être ouvert (dossier protégé, disque plein) : même avertissement, le run continue. Avant, le run s'arrêtait sur une erreur de gfortran.
  - `READ_DUMP0` (reprise, choix 1) : `dump.dat` ouvert en `STATUS='OLD'` ; s'il manque, `Error : dump file not found` et `STOP 1`, comme M020 pour le fichier de données.
- **Pourquoi :** B30 (point A8 de l'audit). Toutes les écritures du dump avaient `IOSTAT=ERROR`, mais `ERROR` n'était jamais lu. Surtout, **gfortran (16.2, macOS) ne signale pas le disque plein** : `IOSTAT` reste à 0, même avec `FLUSH`, et les écritures sont perdues sans bruit (vérifié par un petit programme : 20 000 lignes écrites sur un disque de 4 Mo, aucune erreur, fichier de 4 Mo ; Python voit bien l'erreur `No space left on device`). Vérifier `IOSTAT` ne suffit donc pas, d'où la ligne témoin : les écritures après le disque plein sont perdues, elle manque. `READ_DUMP` s'arrête avant elle : les anciens dumps (sans elle) et les nouveaux se relisent pareil. Le Python (`oracle`) ne lit pas le contenu du dump.
- **Vérification (2026-11-18) :** copie hors du dépôt, `NINFER` = 5000, cas `rel` classique, 1 thread, graine 1113. Disque virtuel de 16 Mo (`hdiutil`) contenant déjà un bon `dump.dat` de 9,4 Mo, run de 40 s qui écrit son premier dump au départ du MCMC :
  - ancien code (M040) : `dump.dat` remplacé par un fichier tronqué de 6,9 Mo, sortie `Dump file created... ok` ;
  - nouveau code : `dump.dat` intact (identique octet par octet), fichier temporaire effacé, avertissement écrit, le run continue.
  - Disque normal : dump identique octet par octet à celui de l'ancien code, plus la ligne témoin (12 octets).
  - Reprise (choix 1) pendant 30 s depuis un dump nouveau, depuis un ancien dump, et d'un dump nouveau avec l'ancien code : `Here we go again`, sorties identiques ligne à ligne entre les trois.
  - Reprise sans `dump.dat` : nouveau code, `Error : dump file not found : dump.dat`, code 1 ; ancien code, `Sequential READ or WRITE not allowed after EOF marker`, code 2.
  - `make test` passe.
- **Limites :** l'appelant écrit encore `ok` après l'avertissement (`Dump file created... Warning ...` puis ` ok`). Arrêt sur `stop.flag` (M039) : si le dump échoue, le message `Stopped on request, dump file written` est quand même écrit, et `oracle run` le croit. Constat : une reprise sans `logfile.dat` s'arrête aussi sur une erreur de gfortran (`Cannot open file 'logfile.dat'`) ; `oracle` le garde toujours avec le dump.
- **Effet sur les résultats :** aucun. Le dump gagne une ligne à la fin.

## M040 — Corrige B29 : `SEARCH_BETA` s'arrête avec un message si aucun pas n'est accepté

- **Date :** 2026-11-18
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`SEARCH_BETA`, un bloc repéré après le calcul de `FAIL`).
- **Changement :** à chaque balayage du réglage des pas, si l'un des pas (`BETA`) est tombé sous 1e-30 en valeur absolue, le code écrit `Error : no MCMC step accepted (step size -> 0)` et `forbidden starting point (priors, bounds) or chi2 not a number`, puis s'arrête (`STOP 1`).
- **Pourquoi :** B29 (point A3 de l'audit du 2026-11-14). Un pas ne tombe à 0 que si aucun essai n'est accepté, même minuscule : le point de départ est alors interdit (prior nul, hors des bornes, χ² NaN). `UPDATE_BETA` divise le pas par 100 à chaque mise à jour sans acceptation : il passe sous 1e-30 en quelques dizaines de mises à jour, d'où un arrêt en moins d'une seconde. Les pas légitimes les plus petits sont de l'ordre de 1e-10 : le seuil ne gêne pas un run normal. M035 et M036 refusent déjà les départs interdits connus ; M040 couvre les autres. `CHANGE_DIRS` appelle aussi `SEARCH_BETA` (chaînes en cours de run) : une chaîne bloquée y arrête maintenant le run au lieu de le faire tourner sans fin. Constat sans correction : `CHANGE_DIRS` ignore le `FAIL` renvoyé par `SEARCH_BETA`.
- **Vérification (2026-11-18) :** copie hors du dépôt, `NINFER` = 5000, 1 thread, graine 1113, quatre variantes : avec et sans M040, avec et sans le contrôle de départ de M036.
  - Départ interdit voulu (cas `rv_seul`, choix 3, prior plat de 2 à 10 Msun sur la masse de l'étoile, masse de départ 1,02), contrôle de M036 neutralisé : sans M040, rien en 60 s (boucle sans fin, B23 reproduite), classique et universel ; avec M040, arrêt avec le message, code 1, en moins d'une seconde, classique et universel.
  - Même cas avec le code complet : refusé par M036 (« outside the bounds »), comme avant.
  - Cas normaux (`rel`, `rel_rv`, `rel_2orb`, `rel_seppa` en classique et universel, `rel_rv_jit_parab` en universel), 90 s chacun avec et sans M040 : `beta initiation successful` partout, et la sortie écran est identique ligne à ligne entre les deux versions (de 68 à 1 072 lignes comparées).
  - `make test` passe.
- **Effet sur les résultats :** aucun pour un run qui démarre normalement. Un run qui tournait sans fin sans rien écrire s'arrête maintenant avec un message.

## M039 — Amélioration : arrêt demandé par un fichier `stop.flag`, avec sauvegarde

- **Date :** 2026-11-17
- **Fichiers :** `Code/mcmc/sub/mcmc.f` (`MCMCP`, une déclaration et un bloc à la fin de la boucle principale).
- **Changement :** à chaque tour de la boucle principale de `MCMCP`, une fois la chauffe faite, le Fortran regarde s'il existe un fichier `stop.flag` dans son dossier. S'il le trouve, il écrit le dump (`WRITE_DUMP`, comme la sauvegarde périodique), retire le fichier, écrit `Stopped on request, dump file written` et s'arrête (`STOP`). Pendant la chauffe (1 000 balayages), il ne le lit pas : une reprise saute la chauffe (`BURNT = .NOT.NEW`, M015), donc un dump écrit avant sa fin serait faux. Le fichier peut être créé à la main (`touch stop.flag`) ou par un lanceur. Le lanceur `oracle run` (à venir dans ce dépôt), sur Ctrl-C ou SIGTERM, crée le fichier, attend au plus 25 s que le Fortran s'arrête, puis le tue s'il ne l'a pas fait (pendant le Levenberg-Marquardt, il n'y a rien à sauvegarder, et il le tue tout de suite). Le Fortran est lancé en ignorant SIGINT et SIGTERM, que le terminal et les gestionnaires de tâches envoient à tous les processus : seul `oracle run` les reçoit et fait l'arrêt propre.
- **Pourquoi :** demandé par Antoine (2026-11-17). Avant, Ctrl-C ou un gestionnaire de tâches tuait le Fortran, et la reprise partait du dernier dump périodique. Avec `FREQDUMP` à 10 millions de balayages, ce dump n'existait souvent pas : le run ne pouvait pas reprendre, ou refaisait tout depuis l'ancien dump. Une simple pause (SIGSTOP) a été écartée : elle garde la mémoire et ne survit ni à un redémarrage ni à la limite de temps d'un serveur de calcul.
- **Serveurs de calcul :** Slurm envoie SIGTERM à tous les processus du job à la limite de temps (ou à `scancel`), puis SIGKILL après `KillWait`, 30 s par défaut. L'arrêt propre tient dans ce délai : moins de 3 s dans l'essai ci-dessous, sur un petit système. Un gestionnaire au délai plus court, ou un gros système dont le dump est long à écrire, peut le couper : la reprise part alors du dump précédent, comme avant. Pour plus de marge, demander le signal avant la limite (Slurm : `#SBATCH --signal=TERM@120`, envoyé 120 s avant).
- **Vérification (2026-11-17) :** système d'exemple à une orbite (astrométrie relative, 2 threads), lancé par `oracle run`, arrêt en pleine récolte (100 001 balayages, aucun dump périodique) :
  - SIGINT à tout le groupe de processus, comme Ctrl-C dans un terminal : arrêt en 2,8 s, `dump.dat` écrit, `stop.flag` retiré, état `stopped`.
  - SIGTERM à tout le groupe puis SIGKILL 30 s après, comme Slurm : même résultat.
  - `oracle run --continue` sur le premier : `Here we go again`, `countok` repart de 1 005 (1 004 à l'arrêt) et `nmod` continue ; arrêté à son tour par `timeout` (SIGTERM), il se sauvegarde de nouveau.
  - `make test` passe.
- **Effet sur les résultats :** aucun. Sans `stop.flag`, le calcul est inchangé (un `INQUIRE` par tour de boucle, soit tous les 100 balayages environ).

## M038 — Corrige B28 : un pas d'essai dont le χ² n'est pas un nombre est refusé

- **Date :** 2026-11-14
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`METRO_STEP`, un bloc repéré).
- **Changement :** la probabilité d'acceptation vaut 0 si le χ² d'essai ou le rapport des priors est NaN, ou si leur produit l'est (0 × infini). Si le χ² courant est lui-même NaN, tout essai fini reste accepté, comme avant, pour que la chaîne puisse en sortir. Sinon la formule est inchangée.
- **Pourquoi :** B28, trouvée à l'audit du 2026-11-14. `ALPHA = MIN(EXP(-0.5d0*(CHI2TRY-CHI2))*RAPQ,1.d0)` : avec gfortran, `MIN(NaN, 1)` vaut 1 (vérifié par un petit programme, de −O0 à −O3). Un pas d'essai à χ² NaN était donc accepté ; la chaîne gardait un χ² NaN, puis acceptait tous les pas suivants. C'est ce qui s'est passé dans les cas dégénérés de la campagne du 2026-11-11 (lignes `chi2 = NaN`), refusés depuis par M032. Le code est compilé sans `-ffast-math` (`Makefile`), donc les tests `X.NE.X` sont fiables (M031 en dépend aussi).
- **Vérification (2026-11-14) :** copies hors du dépôt, cas `rel` (classique, 1 thread, graine 4242).
  - Sans NaN : fichier de résultats identique avant et après M038, octet pour octet après décompression.
  - Avec un NaN injecté (copie où `CHI2SEUL` renvoie NaN une fois sur 1000) : avant M038, le fichier de résultats contient 124 solutions à χ² NaN et des solutions jusqu'à χ² = 131 (pas acceptés à tort ; médiane 25,9) ; après, aucun NaN, χ² maximal 57 (53 sans injection), et les mêmes distributions que sans injection (a = 1,998 ± 0,022 au, e = 0,282 ± 0,009, i = 59,6 ± 0,8°).
  - `make test` passe (370 OK).
- **Effet sur les résultats :** aucun pour un run sans NaN (aucun χ² d'essai NaN dans les 35 runs de la nuit du 13 au 14, 6 cas, classique et universel) ; une chaîne ne peut plus être emportée par un NaN.

## M037 — Corrige B27 : corrélation des erreurs pour les données en (séparation, angle de position)

- **Date :** 2026-11-14
- **Fichiers :** `Code/mcmc/sub/io.f` (`READ_INPUTLINE`, un bloc repéré) ; `Tests/generate_cases.py`, `Tests/run_tests.py`, `Tests/cases/rel_seppa/point.json`, `Tests/reference/reference.json`.
- **Changement :** dans la conversion (Sep, PA) → (Dec, RA), la covariance des erreurs vaut `cos·sin·V_sep − Sep²·sin·cos·V_pa + Sep·(cos² − sin²)·Cov(sep, pa)`. Le code écrivait `Sep` au lieu de `Sep²` devant `V_pa` (terme homogène à des mas au lieu de mas²). Les erreurs `dDec`, `dRA` étaient justes ; seul le coefficient de corrélation ρ était faux.
- **Pourquoi :** B27, trouvée à l'audit du 2026-11-14. Exemple (première mesure du cas `rel_seppa`) : Sep = 16,427 mas, PA = 159,214°, dSep = 1 mas, dPA = 3,488°, ρ(sep, pa) = 0,2 : le vrai ρ(Dec, RA) vaut +0,151, le code donnait −0,163. L'ellipse d'erreur de chaque mesure était mal orientée, donc le χ² faux. La conversion n'était vérifiée par aucun test (`generate_cases.py` : « pas de recalcul indépendant », et `run_tests.py` sautait ce cas).
- **Tests :** `generate_cases.py` refait maintenant la conversion en Python, à partir des valeurs écrites dans le fichier, par la matrice jacobienne (`J C Jᵀ`), et `run_tests.py` compare le χ² du Fortran à celui du modèle Python pour `rel_seppa` comme pour les autres cas. Régénérer les cas ne change que `point.json` de `rel_seppa` (données identiques). Avant la correction, ce contrôle échoue (Python 187,039, Fortran 253,836) ; après, il passe (187,039 des deux côtés, classique et universel). La référence est refaite pour ce cas seulement : 6 valeurs changent (χ², MAP et critère de `rel_seppa`, classique et universel).
- **Vérification (2026-11-14) :** `make test` : 370 OK (368 plus les deux nouveaux contrôles indépendants de `rel_seppa`). Le diff de `reference.json` ne touche que `rel_seppa`.
- **Effet sur les résultats :** **tous les ajustements faits avec des données en (Sep, PA) sont concernés** (formats 5 à 8 de l'astrométrie relative, et du type 4) : leur χ² utilisait une corrélation fausse. L'effet dépend de la forme des ellipses d'erreur. Il existe même pour des erreurs rondes : quand l'erreur radiale égale l'erreur tangentielle (dSep = Sep·dPA) sans corrélation, la vraie corrélation en (Dec, RA) est nulle, alors que le code donnait environ cos(PA)·sin(PA), jusqu'à ±0,5. Les runs concernés sont à refaire pour des résultats exacts. Les données en (Dec, RA) ne sont pas touchées (GG Tau par exemple).

## M036 — Corrige B26 : départ du MCMC contrôlé dans ses propres paramètres (régression de M035 en variables universelles)

- **Date :** 2026-11-14
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`INIT_MCMC`, un bloc repéré, entre les deux parties du bloc M035).
- **Changement :** le contrôle du point de départ du MCMC (M035) ne passe plus par `MRQCHECK` mais par `RATIO_PRIORS(P,P)` : le départ est refusé si `BAD` est vrai ou si le rapport n'est pas strictement positif (prior nul, ou NaN). Le message de M035 est inchangé.
- **Pourquoi :** B26, trouvée à l'audit du 2026-11-14. `MRQCHECK` est une routine du Levenberg-Marquardt. En variables universelles, le L-M garde l'ancienne paramétrisation (`P6 = ln Pq`, `P7 = s·f(e)`) alors que le MCMC utilise `(cos λ, sin λ)/Pq` depuis M011 (le programme principal convertit). Appelée avec les paramètres du MCMC, elle calculait `Pq = EXP(P6)`, soit environ 1 jour, donc des masses dynamiques fausses : avec un prior plat (type 3) sur une masse dynamique, le run était refusé à tort. Les vérifications de M035 ne portaient que sur un prior sur la masse de l'étoile, qui n'utilise pas ces masses. `RATIO_PRIORS` applique exactement les bornes et les priors du MCMC, dans la bonne paramétrisation des deux versions. En classique, L-M et MCMC ont les mêmes paramètres : rien ne changeait. Le bloc M035 de `MRQCHECK` (pas du L-M) est correct et reste.
- **Vérification (2026-11-14) :** copies hors du dépôt.
  - Astrométrie relative + VR, prior plat de 0,5 à 2 Msun sur la masse totale (`1 1`) : avant M036, refusé à tort en universel (« starting point of the MCMC is outside the bounds »), accepté en classique ; après, le MCMC démarre dans les deux versions.
  - VR seules, prior plat de 0,1 à 10 Msun sur la masse de l'étoile : démarre, classique et universel.
  - Départ volontairement hors des bornes (MCMC seul, prior plat de 2 à 10 Msun, masse de départ 1,02) : refusé avec le message, code 1, classique et universel.
  - Cas normaux (`rel`, `rel_rv`, `rel_2orb` en classique ; `rel_rv`, `rel_2orb`, `rel_rv_jit_parab` en universel) : démarrent.
  - `make test` passe.
- **Effet sur les résultats :** aucun en classique ; en universel, les runs avec un prior plat sur une masse dynamique ne sont plus refusés (défaut présent seulement entre M035 et M036, jamais poussé).

## M035 — Corrige B23 : prior plat respecté par le L-M, départ hors des bornes refusé

- **Date :** 2026-11-13
- **Fichiers :** `Code/mcmc/sub/kepellip.f` et `Code/mcmc/sub/kepu.f` (`MRQCHECK`, mêmes blocs dans les deux), `Code/mcmc/sub/mcmc.f` (`INIT_MCMC`).
- **Changement :**
  1. `MRQCHECK` refuse aussi un pas du Levenberg-Marquardt dont une combinaison de masses sort des bornes d'un prior plat (type 3, « Linear »). Les masses sont calculées comme dans `COMPUTE_MAP` (`MUF = Σ BCOF·MDYN + BCOF(0)·m0`). Un tel pas est traité comme les autres pas hors des bornes (λ × 10, M026).
  2. Au démarrage d'un nouveau MCMC, le point de départ passe par `MRQCHECK`. S'il est hors des bornes (période, demi-grand axe, excentricité, masses ou prior plat), le programme s'arrête avec un message et les paramètres (`STOP 1`).
- **Pourquoi :** B23. Le L-M ignorait les bornes d'un prior plat (sa dérivée y est nulle) et pouvait finir là où le prior est nul ; le MCMC n'y acceptait aucun pas et `SEARCH_BETA` bouclait sans fin, sans message. Choix d'Antoine (2026-11-13) : corriger la cause (1) et refuser un départ impossible (2), plutôt qu'une sortie de secours dans `SEARCH_BETA`. Une limite d'itérations du L-M, envisagée d'abord sur un mauvais diagnostic, a été écartée : les cas normaux convergent en 4 à 17 itérations (GG Tau : 510).
- **Vérification (2026-11-13) :** copies hors du dépôt, en classique et en universel.
  - VR de l'étoile seules, prior plat de 0,1 à 10 Msun (le cas qui tournait 8 h sans rien écrire) : le L-M finit sur la borne (`M*` = 0,1 Msun, MAP fini), puis le MCMC démarre (pas trouvés), dans les deux versions. Il ne converge pas en 15 min (R jusqu'à 1,8 en classique, 1,4 en universel), comme avec le prior gaussien : c'est F4, pas B23.
  - Départ volontairement hors des bornes (MCMC seul, prior plat de 2 à 10 Msun, masse de départ 1,02) : arrêt avec le message, code 1.
  - Cas normaux de la campagne (`rel`, `rel_rv`, `rel_rv_prior`, `rel_rv_jit`, `rel_2orb` en classique ; `rel_rv`, `rel_2orb` en universel) : L-M et démarrage du MCMC inchangés, `rel` mené jusqu'au bout.
  - `make test` passe.
- **Effet sur les résultats :** aucun quand le L-M restait dans les bornes ; avec un prior plat, le point de départ du MCMC est désormais dans les bornes.

## M034 — Corrige B25 : `PROBINFER` reçoit une valeur dès le départ

- **Date :** 2026-11-13
- **Fichier :** `Code/mcmc/sub/mcmc.f` (`INIT_MCMC`, branche d'un nouveau run, un bloc repéré).
- **Changement :** `PROBINFER = 0` au démarrage d'un nouveau run, à côté de `INFER = .FALSE.`.
- **Pourquoi :** B25. `PROBINFER` n'était fixé qu'à la convergence (`GELMAN_RUBIN`), mais `WRITE_DUMP` l'écrivait dans `dump.dat` avant : une valeur non initialisée. Elle ne sert qu'une fois `INFER` vrai, et elle est alors recalculée ; une reprise relit `INFER` faux et ne l'utilise pas non plus. Aucun effet sur le calcul, seulement un fichier propre.
- **Vérification (2026-11-13) :** deux copies compilées avec `-finit-real=snan` (réels non initialisés à NaN), cas `rel` de la campagne, dump tous les 100 balayages, arrêt au bout de 3 s (balayage 301, avant la convergence) : `PROBINFER` vaut NaN dans `dump.dat` sans M034, 0 avec. `make test` passe.
- **Effet sur les résultats :** aucun.

## M033 — Corrige B24 : dimension de `P` dans `STORE_SOLUTION`

- **Date :** 2026-11-13
- **Fichier :** `Code/mcmc/sub/io.f` (`STORE_SOLUTION`, routine interne de `WRITE_DISTRIB` et de `WRITE_DISTRIB_DOUBLING`, un bloc repéré dans chacune).
- **Changement :** `P` est déclaré `P(NPAR+2)` au lieu de `P(NPAR)`.
- **Pourquoi :** B24. La routine lit le χ² et le MAP en `P(NPAR+1)` et `P(NPAR+2)`. Les tableaux passés (`PSTART`, `PSAV(:,I)`) ont bien `NPAR+2` éléments, donc le calcul était juste, mais `-fcheck=all` arrêtait le programme à l'écriture des résultats (campagne de tests du 2026-11-11, où il fallait contourner).
- **Vérification (2026-11-13) :** copie compilée avec `-O0 -fcheck=all`, cas `rel` de la campagne (astrométrie relative, `NINFER` = 2000, 2 threads) mené jusqu'au bout : `res.dat.gz` écrit, aucune erreur de bornes, sans le contournement qu'il fallait pendant la campagne. `make test` passe.
- **Effet sur les résultats :** aucun.

## M032 — Masses non déterminées par les données : run refusé au démarrage

- **Date :** 2026-11-12
- **Fichiers :** `Code/mcmc/sub/io.f` (`INIT_SIMU`, juste après la lecture du nombre de priors sur les masses, un bloc repéré) ; `Tests/cases_definition.py`, `Tests/cases/rv_seul/input_*.txt`, `Tests/reference/reference.json`.
- **Changement :** un nouveau run sans astrométrie relative (type 1), sans vitesses radiales du compagnon (type 3) et sans prior sur les masses est refusé : message `Error : the masses are not constrained by the data : stellar RV and/or Gaia astrometry alone (no relative astrometry, no companion RV). Add a prior on the stellar mass (or fix it : type 4).`, puis `STOP 1` (comme M012). Les reprises (`dump.dat`) ne passent pas par ce contrôle.
- **Pourquoi :** trouvé pendant la campagne de tests du 2026-11-11. Avec les vitesses de l'étoile seules, les données ne contraignent qu'une combinaison des masses (fonction de masse, `m_p sin i/M^(2/3)`), pas leur échelle. Les masses dérivaient vers 0 jusqu'à `SIGMA = 0` dans `CHI2SEUL` (`MFRAC(I) = MUF(I)/SIGMA`, 0/0, trouvé avec lldb) : NaN puis arrêt (`epic fail` en variables universelles), ou run qui tourne une heure sans résultat (version classique). Même chose avec l'astrométrie de Gaia seule, dans le dépôt privé (« Gaia seule » et « Gaia + VR » sans prior). Aucun de ces runs n'aboutissait. Un prior, même sur une seule masse, fixe l'échelle. Choix d'Antoine (2026-11-12) : refuser au démarrage plutôt qu'avertir ou s'arrêter en cours de route.
- **Tests :** le cas `rv_seul` (vitesses de l'étoile seules, sans prior) était refusé à son tour. Il reçoit le prior de `rel_rv_prior`, une gaussienne sur la masse de l'étoile (1,0 ± 0,1 Msun) : `cases_definition.py`, puis `generate_cases.py`, qui ne change que les deux fichiers d'entrée de `rv_seul` (données identiques). La référence est refaite pour ce cas seulement : 4 valeurs changent, MAP et critère en classique et en universel. Le MAP baisse de 0,02 et le critère monte de 0,04, soit exactement le terme du prior au point de test, (1,02 − 1,0)²/(2 × 0,1²) = 0,02. Le χ² et les autres contrôles de ce cas sont inchangés.
- **Vérification (2026-11-12) :** copie hors du dépôt, entrées de la campagne. Vitesses de l'étoile seules sans prior : refusé, code 1. Les mêmes avec un prior sur la masse : le run démarre (Levenberg-Marquardt puis MCMC). Astrométrie relative sans prior : le run démarre. `make test` : 368 OK.
- **Effet sur les résultats :** aucun pour les runs qui aboutissaient ; les runs voués à l'échec s'arrêtent tout de suite avec la cause.

## M031 — Levenberg-Marquardt : arrêt expliqué quand le χ² n'est pas un nombre

- **Date :** 2026-11-12
- **Fichier :** `Code/mcmc/sub/mrqfit.f` (`MRQMIN`, un bloc repéré).
- **Changement :** la ligne `if (chi2.ne.chi2) stop` arrêtait le programme sans rien écrire. Elle écrit maintenant à l'écran (`out.txt`) que le χ² du Levenberg-Marquardt est NaN, au point de départ ou après un redémarrage, qu'il faut vérifier l'orbite de départ et les priors, puis le numéro d'itération, λ et les paramètres en cours (variables de l'ajustement). Le programme s'arrête avec le code 1 (`STOP 1`, comme M012) au lieu de 0.
- **Pourquoi :** trouvé pendant la campagne de tests du 2026-11-11 : un run qui s'arrête ainsi laisse un `out.txt` qui finit sur une ligne de χ² ordinaire, sans cause visible, et un code de sortie 0. Choix d'Antoine (2026-11-12) : un message clair plutôt qu'un redémarrage automatique du L-M. Un pas d'essai dont le χ² est NaN est déjà refusé (`MRQSTEP`) ; le χ² en cours ne peut donc être NaN qu'au point de départ ou après `MRQRENEW`, et l'arrêt n'a lieu que si le pas suivant est aussi refusé (un pas accepté remplace le NaN, car `CTRY.GT.CHI2` est faux).
- **Vérification (2026-11-12) :** copie hors du dépôt, vitesses seules (`rv_seul` de la campagne). χ² forcé à NaN au point de départ : le premier pas est accepté et le L-M finit normalement, comme avant. χ² forcé à NaN à la troisième itération : le message et les paramètres sont écrits, code de sortie 1. `make test` passe.
- **Effet sur les résultats :** aucun ; un arrêt qui était silencieux est expliqué.

## M030 — Prior de Heaviside (masses positives) sans débordement de l'exponentielle

- **Date :** 2026-11-12
- **Fichiers :** `Code/mcmc/sub/mcmc.f`, `USER_PRIOR_RATIO` (rapport du MCMC, type 6) et la routine des dérivées du prior pour le Levenberg-Marquardt (type 6, `F`, `DF`, `D2F`).
- **Changement :** mêmes formules, réécrites pour ne jamais calculer `EXP(-S)` avec `S` très négatif (masse très négative devant la largeur du prior) : `LOG(1+EXP(-S)) = MAX(-S,0) + LOG(1+EXP(-|S|))`, et dans les dérivées `ES/(1+ES) = 1/(1+EXP(S))`. Le MAP (`USER_PRIOR_MAP`) était déjà protégé par M010.
- **Pourquoi :** trouvé pendant la campagne de tests du 2026-11-11. Pour `S < -709`, `EXP(-S)` vaut l'infini : dans le Levenberg-Marquardt, `F` est infini et `DF`, `D2F` valent infini sur infini, donc NaN, ce qui arrête le programme sans message (`MRQMIN`, ligne `if (chi2.ne.chi2) stop`). Dans le MCMC, le rapport reste infini quand le point de départ a une masse très négative, comme avant (tout pas qui s'en éloigne est accepté) ; `make test` contient ce cas.
- **Vérification (2026-11-12) :** en Python, les nouvelles formules donnent les mêmes `F`, `DF`, `D2F` et le même rapport que les anciennes partout où les anciennes sont finies (écart relatif maximal 4e-16, pour |S| < 600), et des valeurs finies au-delà (`S = -2000` : `F` = 4000, `DF` et `D2F` de l'ordre de 1e-6). `make test` passe.
- **Effet sur les résultats :** aucun quand les masses restent raisonnables ; supprime un arrêt silencieux quand le Levenberg-Marquardt passe par une masse très négative.

## M029 — Graine du générateur aléatoire tirée puis écrite (corrige B14)

- **Date :** 2026-11-11
- **Fichiers :** `Code/mcmc/sub/utils.f` (nouvelle routine `INIT_SEED`, en fin de fichier), `Code/mcmc/main/astrom_mcmco.f` et `astrom_univ_mcmco.f` (appel juste après `INIT_SIMU`).
- **Changement :** au démarrage, la graine est lue dans la variable d'environnement `ORACLE_SEED` si elle existe (un entier) ; sinon elle est tirée au hasard, comme avant. Dans les deux cas, elle sert à initialiser le générateur (`RANDOM_SEED(PUT=...)`) et elle est écrite à l'écran (`out.txt`) et dans `logfile.dat` : `Random seed (ORACLE_SEED) = ...`. Pour refaire un run à l'identique : `ORACLE_SEED=<graine> ./start.sh` (même nombre de threads). Une reprise depuis `dump.dat` tire une nouvelle graine (l'état du générateur n'est pas sauvegardé) et l'écrit aussi.
- **Pourquoi :** B14. Choix d'Antoine (2026-11-11) parmi quatre options (graine tirée et écrite, graine donnée en entrée, graine fixe, rien) : rien ne change pour l'utilisateur, deux runs restent deux tirages indépendants, et on peut refaire un run précis. La variable d'environnement évite d'ajouter une question à l'entrée (interface et `start.sh` inchangés).
- **Vérification (2026-11-11) :** `rel_vr_planete` (`--rapide`, 2 threads OpenMP), trois runs complets : deux avec `ORACLE_SEED=12345` donnent des fichiers de résultats identiques au bit près (après décompression ; les `.gz` ne diffèrent que par leur date) ; un sans graine donne un autre échantillon et écrit sa graine (678312113) dans `out.txt` et `logfile.dat`. `make test` passe.
- **Effet sur les résultats :** aucun (un run sans `ORACLE_SEED` reste un tirage au hasard).

## M028 — Tests : orbite quasi parabolique, et outils de test sûrs près de e = 1 (sans repère)

- **Date :** 2026-11-11
- **Fichiers :** `Tests/` seulement : `cases_definition.py` (orbite `PARAB` et cas `rel_rv_jit_parab`), `generate_cases.py` (borne sur e par cas, `emax`), `model.py` et `jacobian.py` (Kepler sûr près de e = 1), `driver.f` et `run_tests.py` (extrapolation de Richardson), référence complétée.
- **Changement :**
  1. Cas `rel_rv_jit_parab` : astrométrie relative, vitesses de l'étoile et jitter, orbite quasi parabolique (a = 3 au, e = 0,97, passage au périastre au milieu des mesures ; point d'évaluation à e = 0,99, borne de l'entrée relevée à 0,995), en variables classiques et universelles. Il a sa propre graine (option `seed` de `generate_cases.py`, ajoutée le même jour) : ses données ne dépendent pas de sa place dans la liste, qui n'est pas la même dans le dépôt privé ; les cas existants gardent la graine commune et ne changent pas.
  2. `model.kepler` part de M + 0,85 e·signe(sin M) (Danby 1987) au lieu de M : le Newton pouvait diverger près de e = 1.
  3. `jacobian.universal_forward` passe, pour une orbite liée, par l'équation de Kepler elliptique au lieu d'un Newton en variable universelle parti de s = dt/q, qui ne convergeait pas toujours quand q est petit (ln J sautait de 14 entre deux points voisins dans un cas Gaia du dépôt privé).
  4. `driver.f` : la dérivée de référence par différences finies utilise l'extrapolation de Richardson, (4 D(h/2) − D(h))/3, qui retire l'erreur de troncature en h² ; tolérance de bruit d'arrondi multipliée par 3 (`run_tests.py`). Au pas h seul, cette erreur atteignait 3·10⁻⁴ (au-dessus de la tolérance de 10⁻⁴) sur une orbite à e = 0,99 du dépôt privé.
- **Pourquoi :** vérifier le code près de e = 1, en particulier les variables universelles, faites pour ces orbites, sans que les outils de test eux-mêmes se trompent.
- **Résultat :** le code passe tous les contrôles sur ce cas (χ² contre le modèle Python, classique = universel, rapport des priors, dérivées du L-M, éléments reconstruits). Résultats des anciens cas inchangés.
- **Vérification :** `make test` : 368 OK. Référence complétée pour les 2 nouvelles entrées ; le diff ne contient que des ajouts.
- **Effet sur les résultats :** aucun (tests seulement).

## M027 — Excentricité de départ nulle acceptée

- **Date :** 2026-11-11
- **Fichier :** `Code/mcmc/sub/io.f` (`ENTER_ORBITS`, un bloc repéré).
- **Changement :** une excentricité d'entrée inférieure à 0,001 (en pratique e = 0) est remplacée par 0,001, avec le message `Eccentricity < 0.001 replaced by 0.001`.
- **Pourquoi :** les paramètres 2 et 3 valent e·(cos ϖ, sin ϖ)/√(1−e²). Avec e = 0, `S` = 0 et le code divise par `S` partout (`EXQ = EXC/S`, `CW = P2/S`, `SW = P3/S` dans `ELEMENTS`, `RATIO_PRIORS`, les routines du L-M et `MRQELEMENTS`) : 0/0 donne des NaN, et le programme s'arrête dans `KEPLER0` (message `glup`, puis `IEEE_INVALID_FLAG`). Le MCMC ne tombe jamais exactement sur e = 0 ; seule l'orbite de départ, donnée par l'utilisateur, peut le faire. Corriger à l'entrée évite de toucher à toutes ces divisions. Une orbite circulaire est une valeur de départ naturelle quand on ne connaît pas l'excentricité.
- **Mesure :** vitesses de HD 114762, départ à e = 0 : arrêt sur `glup` avant M027 ; après, même solution que les autres départs (χ² + Fmap = 80,048, e = 0,337).
- **Effet sur les résultats :** aucun, sauf pour les départs à e = 0, qui fonctionnent désormais.
- **Vérification :** `make test` passe.

## M026 — Levenberg-Marquardt : un pas hors des bornes augmente λ au lieu de tout recommencer

- **Date :** 2026-11-11
- **Fichier :** `Code/mcmc/sub/mrqfit.f` (`MRQSTEP` : constante `LAMBDAX` = 10⁸ et branche des pas hors des bornes, deux blocs repérés).
- **Changement :** quand le pas d'essai sort des bornes (`MRQCHECK` : e, P, a ou masse de l'étoile hors de leurs limites, ou paramètre NaN) ou donne un χ² NaN, il est refusé comme un pas qui dégrade le χ² : λ est multiplié par 10 et le L-M repart du point courant. `MRQRENEW` (nouvelle solution : orbite d'entrée, Ω, ω et λ tirés au hasard) n'est appelé que si λ dépasse `LAMBDAX`. C'est l'idée laissée en commentaire par Hervé à cet endroit (boucle `DO WHILE(CHK)` qui augmente λ).
- **Pourquoi :** avant M026, chaque pas hors des bornes appelait `MRQRENEW`, qui repart de l'orbite d'entrée et perd tout le chemin parcouru, avec λ remis à 10⁻³. Loin de la solution, le premier pas de Gauss-Newton sort souvent des bornes (pas de 80 sur e cos ω), d'où un nouveau redémarrage, et ainsi de suite : la boucle ne s'arrête jamais. Constaté le 2026-11-11 sur les vitesses de HD 114762, départ à i = 14,8° (période de départ 86,8 j) : 1,7 million de redémarrages en 60 s sans un seul pas accepté, χ² oscillant autour de 513 pendant plus de 10 min. Un grand λ ramène le pas vers la plus grande pente, avec une longueur qui tend vers 0, donc le pas finit toujours par rester dans les bornes si le point courant y est.
- **Mesures (mêmes 7 départs que M025) :** le départ à i = 14,8° se termine en moins d'une seconde (faux minimum, χ² + Fmap = 509, mauvaise période de départ). Les autres départs finissent au même point qu'avec M025 seul. Avec M026 sans M025, la boucle sans fin disparaît aussi, mais les bons départs restent à χ² + Fmap ≈ 89,75 : les deux corrections sont indépendantes.
- **Effet sur les résultats :** seulement le L-M. Les redémarrages au hasard deviennent rares, ce qui devrait aussi réduire l'effet de B14 sur le L-M (sur GG Tau, un passage dans un autre minimum environ une fois sur cent, venu de ces redémarrages ; non remesuré).
- **Vérification :** `make test` passe, essai GG Tau compris.

## M025 — Levenberg-Marquardt : jitter des vitesses par l'information de Fisher

- **Date :** 2026-11-11
- **Fichier :** `Code/mcmc/sub/mrqfit.f` (`MRQRADVEL`, trois blocs repérés).
- **Changement :** dans la matrice `ALPHA` du L-M, la dérivée seconde par rapport à ln σ_jit est remplacée par son espérance (information de Fisher) : `DV*DV` devient `SIGV2P` = σᵢ² + σ_jit² dans le terme diagonal, et les termes croisés jitter × autres paramètres (`2*DV*FS*FS*DVF`), nuls en espérance, sont retirés. `BETA` (le gradient), `MRQMAP` et le χ² ne changent pas : le minimum cherché est le même.
- **Pourquoi :** pour les paramètres orbitaux, `ALPHA` ne garde que les dérivées premières (Gauss-Newton), ce qui la rend toujours positive. Pour le jitter, elle prenait la vraie dérivée seconde, négative dès que les résidus dépassent σ (orbite encore fausse). Avec un terme diagonal négatif, augmenter λ ne ramène pas le pas vers une descente : le L-M s'arrêtait avant le minimum, le jitter glissant vers 0. L'information de Fisher est l'équivalent de Gauss-Newton pour un paramètre de variance ; avec le terme de normalisation de `MRQMAP`, le total vaut 2σ_jit⁴/(σᵢ²+σ_jit²)² par mesure, toujours positif.
- **Mesures (vitesses de Lick de HD 114762, 74 mesures, L-M seul, 8 points de départ, 2026-11-11) :**
  - avant : les départs proches de la solution s'arrêtent à χ² + Fmap = 89,750 à 89,761 (valeurs toutes différentes), avec σ_jit = 1,2 ± 1,6 m/s ;
  - après : les quatre mêmes départs donnent tous χ² + Fmap = **80,048** (identique à 9 chiffres), σ_jit = 13,2 ± 3,0 m/s. P, e, ω et la masse bougent de moins de leur erreur ;
  - un départ (i = 14,8°, période de départ 86,8 j) tournait sans fin avant M025 (voir M026) ; après, il se termine ;
  - les départs avec une mauvaise période (80 à 87 j au lieu de 84 j) finissent comme avant dans un faux minimum (χ² + Fmap ≈ 440 à 460, σ_jit ≈ 0,3 km/s qui absorbe les résidus) : le L-M reste une méthode locale.
- **Effet sur les résultats :** seulement le L-M avec jitter des vitesses (point de départ du MCMC et choix 4). Le MCMC n'utilise pas `ALPHA`. Les barres d'erreur du L-M (inverse de `ALPHA`) sur le jitter et les paramètres corrélés changent.
- **Vérification :** `make test` passe (les dérivées testées par différences finies sont celles de `BETA`, inchangées).

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

- **Vérification du nouveau remplissage (2026-11-11), avec des copies où le seuil est artificiellement haut :**
  - `arc_court` (`--rapide`), `ESSMIN` = 12 000 : ESS 5 760 au premier remplissage (p = 7,2, toutes les positions gardées), puis 11 370 à p = 0,24 et 15 710 à p = 0,114 ; arrêt normal, tous les éléments à moins de 1,4 σ de la vérité. Espacer la récolte augmente bien l'ESS, presque en proportion.
  - `rel_vr_planete` (`--rapide`), `ESSMIN` = 50 000 (inatteignable avec 20 000 solutions) et `PESSMIN` = 0,05 pour abréger : p passe de 250 à 0,18 puis à 0,05, sous `MINPROB` ; à ce plancher, le run s'arrête avec l'avertissement prévu (ESS 18 340), tous les éléments à moins de 1,8 σ de la vérité. L'ESS plafonne vers 18 000 dès p = 0,18 : les solutions sont alors déjà presque indépendantes.
  - Coût : un remplissage demande environ `NINFER` × 100 / p balayages (6,4 millions au total pour le second essai). Avec `PESSMIN` = 0,01 et `NINFER` = 50 000, le dernier remplissage peut atteindre 500 millions de balayages : à garder en tête pour les cas où l'ESS reste faible.
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
