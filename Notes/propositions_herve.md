# Trois anomalies d'ORACLE à trancher ensemble (B9, B11, B12)

Antoine Lacquement, novembre 2026. Branche `corrections`.

Le journal `MODIFICATIONS.md` liste les anomalies trouvées et les corrections déjà faites (M005 à M011). Chaque ligne modifiée y est décrite et repérée dans le code par `[MOD Mxxx]`, avec l'ancienne ligne recopiée. `make test` les vérifie.

Pour les anomalies ci-dessous, la correction dépend de ce que le code est censé faire. En attendant ton avis, j'ai pris une décision pour chacune, notée « Décision prise » : elle est facile à changer si tu vois les choses autrement. Pour chacune : le constat, les options, la décision et la question qui reste.

## B9 — Prior sur la masse de l'étoile seule quand les masses ne sont pas ajustées

**Constat.** Avec une seule orbite et sans vitesses de l'étoile, `MULTIPLA` est faux : seule la masse dynamique `m0 + m1` est un paramètre. Si l'utilisateur met un prior sur `m0` seul (coefficients `1 0`), `COMPUTE_MAP` et `RATIO_PRIORS` ignorent le terme `BCOF(0)*MSTAR` (`kepellip.f:188` et `:556`, même schéma dans `kepu.f`). Le prior est évalué sur une masse nulle et n'agit pas, sans aucun message.

Cas réel : 51 Eri (relative et vitesses de la planète, prior 1,55 ± 0,03 Msun sur l'étoile). Résultats identiques avec et sans ce prior. Les `a` et `P` ajustés donnent `a³/P²` ≈ 1,36 Msun. Avec les coefficients `1 1` (masse totale), le prior agit bien.

**Options.**
1. Refuser ce prior à la lecture, avec un message qui explique d'utiliser les coefficients `1 1` (masse totale).
2. L'interpréter comme un prior sur la masse totale, avec un avertissement. C'est presque juste quand le compagnon est léger, mais c'est une approximation cachée.
3. Garder le comportement actuel, mais avertir.

**Décision prise : option 1 (M012, 2026-11-08).** Le programme s'arrête avec un message qui demande le même coefficient pour `m0` et `m1`. Aucun résultat valide ne change, et un prior ne peut plus être sans effet en silence.

Les vitesses de la planète ne suffisent pas à séparer les masses : ORACLE ne connaît que des vitesses **relatives** à l'étoile (type 3), dont l'amplitude `(2π/P)·a sin i/√(1−e²)` ne dépend que de la masse totale. Une vitesse **absolue** du compagnon, elle, dépend de `m⋆/(m⋆+m_p)`. D'où la fonctionnalité F1 du journal, à faire ensuite : des vitesses absolues pour l'étoile et pour les compagnons, avec une vitesse systémique commune. Avec elles, un prior sur `m0` seul redeviendra valable.

**Question.** Es-tu d'accord avec ce refus, et avec l'ajout des vitesses absolues des compagnons ?

## B11 — Rapport des priors de masse : remplacer ou multiplier le prior par défaut ?

**Constat.** Dans `RATIO_PRIORS` (les deux versions), chaque prior de masse `u(μ)` multiplie le rapport par `u(μ')/u(μ)`, puis par `M'/M` pour chaque masse dynamique concernée et par `m0'/m0` si `m0` est concerné. Le commentaire du code justifie ces facteurs : le prior `u(μ)` **remplace** le prior logarithmique par défaut le long de la masse. `COMPUTE_MAP`, lui, **multiplie** le prior par défaut par `u(μ)`, sans ces facteurs. Les deux routines ne décrivent donc pas le même prior.

S'y ajoutent deux effets :
- les facteurs `M'/M` sont aussi appliqués pour les priors de type 6 (masses positives), qui ne devraient être qu'une contrainte ;
- quand plusieurs priors portent sur la même masse, le facteur est compté plusieurs fois.

`make test` mesure l'écart exactement (contre un jacobien calculé numériquement), et l'écart est entièrement expliqué par ces facteurs. Effet sur un run réel à deux orbites, estimé par repondération : `m0` vaut 0,563 Msun contre 0,536 avec le prior que donne `COMPUTE_MAP`. Les orbites externes changent peu. Le prior échantillonné favorise donc les grandes masses.

**Options.**
- **A. « Remplacer »** (l'intention du commentaire) : garder un seul facteur `M'/M` par masse et par prior utilisateur, pas pour le type 6, et ajouter le terme `+ln M` correspondant dans `COMPUTE_MAP` pour que le MAP décrive le même prior.
- **B. « Multiplier »** (ce que fait `COMPUTE_MAP`) : supprimer les facteurs `M'/M` et `m0'/m0` de `RATIO_PRIORS`. Le rapport est alors cohérent avec le MAP, sans autre changement.

**Décision prise : B, « multiplier » (M013, 2026-11-08).** Les facteurs `M'/M` et `m0'/m0` sont retirés du rapport. Un seul `m0'/m0` reste, explicitement : c'est le jacobien de `ln m0` pour un prior plat en `m0`, comme dans `COMPUTE_MAP`. Le rapport du MCMC suit maintenant exactement `COMPUTE_MAP` et le jacobien dans les 22 cas de `make test`.

Raisons : c'est le seul choix bien défini quand plusieurs priors portent sur des combinaisons de masses différentes (« remplacer » le prior par défaut n'y a plus de sens clair), et il se dit simplement : « les priors de l'utilisateur s'ajoutent aux priors par défaut ». Avec un prior étroit, A et B donnent presque le même résultat ; la différence compte pour les priors larges.

**Question.** Voulais-tu A ? Si oui, comment le définir avec plusieurs priors sur des combinaisons différentes ? Et le prior par défaut sur `m0` : il était plat (rien dans `COMPUTE_MAP`). **Décision prise : logarithmique (M014)**, comme pour `a`, `P`, `q` et `P_q`, puisque `m0` est aussi une grandeur d'échelle. Es-tu d'accord ?

## B12 — Fausse convergence de Gelman-Rubin au départ

**Constat.** Les 8 chaînes partent toutes du même point, la solution du LM (`INIT_MCMC`, `CHN(C) = CH0`). Moyennes et variances sont cumulées depuis le premier pas. Au premier contrôle, l'écart entre chaînes `B` est presque nul, donc `R` ≈ 1, et `TS = LONG*NCH*MIN(VARP/B,1)` est énorme. Le critère d'inférence (`RINF`, `TSINF`) est atteint avant que les chaînes aient exploré.

Cas mesuré : `rel_vr_planete` en version classique (relative et vitesses de la planète). Inférence dès le premier contrôle dans 5 runs sur 5. Les distributions sont trop étroites ou asymétriques, et `P` est entre 1,1 et 3,1 sigma de la vérité selon le run. C'est sans doute le même mécanisme qui a fait déclarer convergé un run universel coincé à χ² ≈ 72 000 sur un arc au périastre (avant M011). Les longs runs réels (GG Tau, 51 Eri) ne semblent pas touchés : l'inférence y arrive après des dizaines de millions de pas, quand les chaînes se sont séparées.

**Options** (combinables).
1. **Départs dispersés** : décaler chaque chaîne du point du LM de quelques pas `BETA` (déjà calculés par `SEARCH_BETA`) le long de chaque direction. C'est l'hypothèse de base de Gelman-Rubin : des départs plus dispersés que la distribution cherchée.
2. **Période de chauffe** : remettre à zéro les moyennes et covariances cumulées à un instant donné (par exemple au premier contrôle, ou à chaque doublement de `LONG`), pour que `R` ne porte que sur la suite.
3. **Seconde moitié seulement** : calculer `R` sur la seconde moitié des chaînes, comme le recommandent Gelman et Rubin. Avec des sommes cumulées, cela revient à l'option 2 avec des remises à zéro aux doublements.

**Décision prise : 1 et 2, plus un second remplissage de la réserve (M015, 2026-11-08).**
- Départs dispersés de 2 pas `BETA` autour du point du LM (la chaîne 1 reste au point du LM).
- Chauffe de 1 000 balayages sans Gelman-Rubin, puis remise à zéro des moyennes et covariances.
- La récolte commence toujours à R < 1,1, ce qui donne un aperçu pendant le run. Mais au premier passage sous 1,01, la réserve est réécrite dans l'ordre, en un temps égal à celui mis pour converger. À l'arrêt, toutes les solutions datent d'après la convergence.

En vérifiant, une seconde anomalie est apparue (**B15**, corrigée par M016) : une fois R < 1,01, l'intervalle entre deux contrôles valait `MIN(NGROUP, NINT(LONG/100/NGROUP)·NGROUP)`, soit 0 avant 5 000 balayages. Les chaînes ne bougeaient plus et les mêmes 8 positions étaient gardées jusqu'à remplir la réserve (9 solutions distinctes sur 20 001 dans un run). Les longs runs réels ne sont pas touchés.

Résultat sur `rel_vr_planete` (5 runs) : 0 run correct sur 5 avant, 5 sur 5 après, avec la même distribution à chaque fois.

Pourquoi le second remplissage : avec le remplacement aléatoire, ajouter `NMAX` solutions après la convergence en laisse environ 37 % d'avant (e⁻¹). Un run qui converge vite s'arrête dès la réserve pleine, avec surtout des solutions d'avant la convergence. Une alternative envisagée : une « boîte glissante » (on remplace toujours la plus ancienne solution), qui donne le même résultat final et un aperçu continu. Le second remplissage en est la version à un seul passage, avec deux lignes dans ton code.

**Questions.** Les départs identiques étaient-ils voulus, par exemple pour éviter de partir dans un mauvais minimum ? Quelle dispersion jugerais-tu raisonnable ? Pourquoi le remplacement aléatoire dans la réserve plutôt que des positions à intervalle régulier : était-ce pour la mémoire, ou pour une autre raison ? D'accord avec un second remplissage en un temps égal au temps de convergence, plutôt que 4 fois ? Faut-il garder l'écriture de `results.dat.gz` tous les 10 millions de pas, maintenant que le fichier provisoire `xxx_partial.dat.gz` (M017) suit le run dès le début ?

## Critère d'arrêt : TS remplacé par l'ESS de l'échantillon gardé (M019)

**Constat.** Avec les formules de `GELMAN_RUBIN`, `TS = LONG·NCH·min(VARP/B, 1)` vaut environ `NCH/(R² − 1)`. Avec 8 chaînes, `TS` > 1 000 revient à R < 1,004 : c'était le vrai critère d'arrêt, plus strict que `RCONV`. Sur GG Tau, environ 300 millions de pas (30 % du run) se sont écoulés entre R < 1,01 et l'arrêt. Et `TS` mesure l'écart entre chaînes, pas la qualité de l'échantillon gardé.

**Mesures (2026-11-09, `arc_court`, `--rapide`).**
- Mêmes distributions, à 0,05 σ près, avec `TS`, sans `TS` (2 runs), et avec un second remplissage 4 fois plus long.
- Tirages indépendants (ESS) dans l'échantillon final : minimum 700 à 1 000 sans `TS`, 1 300 avec, 14 500 avec le facteur 4. C'est la durée du second remplissage qui fixe l'ESS, pas `TS`.

**Décision prise (M019).**
- La convergence demande R < 1,01 seulement (Gelman & Rubin 1992).
- À la fin du second remplissage, l'ESS de chaque paramètre est calculé sur la réserve, par moyennes de lots (Flegal & Jones 2010).
- Si le plus petit ESS atteint 400 (seuil de Vehtari et al. 2021, avec R̂ < 1,01), le run s'arrête. Sinon, la réserve est remplie à nouveau avec une récolte plus espacée.

**Questions.**
- Pourquoi `TS` (Ford 2006) plutôt que R seul : as-tu vu des cas où R < 1,01 ne suffisait pas ?
- Le seuil de 400 te paraît-il suffisant pour nos usages (intervalles à 68 et 95 %) ?

## Directions du MCMC recalculées tôt (M024)

**Constat.** Les directions de proposition ne sont recalculées que tous les `NDIAG` = 500 000 balayages. Avant le premier recalcul, les chaînes avancent le long des axes des paramètres : quand ceux-ci sont très corrélés, les pas restent minuscules et la convergence peut demander des millions de balayages.

**Décision prise (M024).** Recalcul aussi en fin de chauffe, puis à 10 000 et 100 000 balayages après la chauffe ; le recalcul tous les 500 000 reste.

**Mesures.** `arc_court` : même distribution, ESS 3 à 6 fois plus grand. `arc_tres_court` : ni gain ni perte. Recalculer à chaque doublement de la longueur s'est révélé moins bon (récolte jusqu'à 2,6 fois plus tardive).

**Questions.**
- Pourquoi 500 000 balayages : as-tu vu des cas où un recalcul précoce, sur une covariance encore mal estimée, posait problème ?
- Les recalculs peuvent tomber pendant la récolte quand elle commence tôt (sur `rel_vr_planete`, elle commence juste après la chauffe). Le code d'origine le faisait déjà tous les 500 000 balayages. Comme le nombre de recalculs pendant la récolte reste fini, l'adaptation s'arrête et la loi échantillonnée n'est pas faussée à la longue ; mais geler les directions une fois la récolte commencée serait plus propre. Qu'en penses-tu ?

## B14, pour mémoire

La graine du générateur aléatoire n'est jamais fixée : deux lancements identiques ne donnent pas le même résultat, et le LM lui-même en dépend (redémarrages au hasard de `MRQRENEW`). Préfères-tu une graine fixe, une graine donnée en entrée, ou une graine tirée puis écrite dans la sortie pour pouvoir refaire un run ?

## Autres conventions à discuter (pas des erreurs)

Relevées le 2026-11-08 en relisant les priors. Le code les applique de façon cohérente ; la question est seulement de savoir si ce sont les bons choix.

1. **Prior sur e en version universelle.** `COMPUTE_MAP` n'a pas de terme en `e` : le prior est uniforme en `e` jusqu'à la borne donnée par l'utilisateur. Si cette borne dépasse 1, la part a priori d'orbites hyperboliques dépend de ce choix arbitraire (un tiers pour une borne à 1,5, la moitié pour 2). Testé le 2026-11-08 sur deux arcs autour d'un périastre à e = 0,98 (2 ans, puis 1 an avec 8 points) : bornes 1,5 et 3, mêmes résultats (36 % puis environ 50 % d'hyperboles, vérité retrouvée). Les données y contraignent assez e pour que la borne ne joue pas. Elle ne compterait que pour des données qui contraignent très peu e.
2. **Seuils de Gelman-Rubin** (`commons.f`) : la récolte commence à `RINF = 1,1` et le MCMC s'arrête à `RCONV = 1,01`, ce qui correspond aux recommandations actuelles (R̂ < 1,01). Depuis M015, les solutions récoltées entre 1,1 et 1,01 ne restent plus dans l'échantillon final. Comme `TS` ≈ 8/(R² − 1), `TSINF` = 100 revient à R < 1,04 : c'est ce seuil, et non 1,1, qui lance la récolte. Il ne sert plus qu'à l'aperçu. Depuis M019, l'arrêt ne dépend plus de `TS` (voir plus haut).
3. **Prior conjoint logarithmique en a (ou q) et en P (ou P_q).** Il implique un prior particulier sur la masse dynamique, auquel s'ajoute un prior de masse de l'utilisateur (choix « multiplier », M013). Une alternative courante : logarithmique en a et en masse.
4. **Prior du jitter** `1/(j + j0)` avec `j0 = 0,1 km/s` fixé dans le code (`SIGV0`). Convention standard (Gregory 2005), mais `j0` devrait sans doute dépendre des données, et ne conviendra pas pour un jitter Gaia en mas.
