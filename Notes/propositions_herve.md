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

**Je penche pour B.** C'est le comportement le plus simple à expliquer à l'utilisateur : « le prior que vous donnez s'ajoute aux priors par défaut ». La correction est petite, et `make test` la vérifie directement. A se défend aussi, mais elle est plus délicate quand plusieurs priors portent sur des combinaisons de masses différentes : « remplacer » le prior par défaut n'y est plus bien défini.

**Question.** Quel prior voulais-tu : A ou B ? Si A, comment le définir avec plusieurs priors sur des combinaisons différentes ?

## B12 — Fausse convergence de Gelman-Rubin au départ

**Constat.** Les 8 chaînes partent toutes du même point, la solution du LM (`INIT_MCMC`, `CHN(C) = CH0`). Moyennes et variances sont cumulées depuis le premier pas. Au premier contrôle, l'écart entre chaînes `B` est presque nul, donc `R` ≈ 1, et `TS = LONG*NCH*MIN(VARP/B,1)` est énorme. Le critère d'inférence (`RINF`, `TSINF`) est atteint avant que les chaînes aient exploré.

Cas mesuré : `rel_vr_planete` en version classique (relative et vitesses de la planète). Inférence dès le premier contrôle dans 5 runs sur 5. Les distributions sont trop étroites ou asymétriques, et `P` est entre 1,1 et 3,1 sigma de la vérité selon le run. C'est sans doute le même mécanisme qui a fait déclarer convergé un run universel coincé à χ² ≈ 72 000 sur un arc au périastre (avant M011). Les longs runs réels (GG Tau, 51 Eri) ne semblent pas touchés : l'inférence y arrive après des dizaines de millions de pas, quand les chaînes se sont séparées.

**Options** (combinables).
1. **Départs dispersés** : décaler chaque chaîne du point du LM de quelques pas `BETA` (déjà calculés par `SEARCH_BETA`) le long de chaque direction. C'est l'hypothèse de base de Gelman-Rubin : des départs plus dispersés que la distribution cherchée.
2. **Période de chauffe** : remettre à zéro les moyennes et covariances cumulées à un instant donné (par exemple au premier contrôle, ou à chaque doublement de `LONG`), pour que `R` ne porte que sur la suite.
3. **Seconde moitié seulement** : calculer `R` sur la seconde moitié des chaînes, comme le recommandent Gelman et Rubin. Avec des sommes cumulées, cela revient à l'option 2 avec des remises à zéro aux doublements.

**Je propose 1 et 2 ensemble.** 1 rend le test valable. 2 empêche les premiers pas, encore corrélés au point de départ, de peser sur `R`. Je vérifierais sur `rel_vr_planete` (5 runs) que l'inférence arrive plus tard et que la vérité est retrouvée, puis sur GG Tau que rien ne change.

**Questions.** Les départs identiques étaient-ils voulus, par exemple pour éviter de partir dans un mauvais minimum ? Quelle dispersion jugerais-tu raisonnable ?

## B14, pour mémoire

La graine du générateur aléatoire n'est jamais fixée : deux lancements identiques ne donnent pas le même résultat, et le LM lui-même en dépend (redémarrages au hasard de `MRQRENEW`). Préfères-tu une graine fixe, une graine donnée en entrée, ou une graine tirée puis écrite dans la sortie pour pouvoir refaire un run ?
