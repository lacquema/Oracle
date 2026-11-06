# Choix de la paramétrisation des variables universelles

Comparaison, sur des systèmes synthétiques, de trois façons d'ajuster une orbite :

| Version | Paramètres 1 à 7 | Code |
|---|---|---|
| Classique | `ln a`, `e/√(1-e²)·(cos, sin)ϖ`, `(1-e²)^¼ sin(i/2)·(cos, sin)Ω`, `(cos, sin)λ / P` | `kepellip.f`, tel quel |
| Universelle actuelle | `ln q`, `e·(cos, sin)ϖ`, `sin(i/2)·(cos, sin)Ω`, `ln P_q`, `s·f(e)` | `kepu.f`, corrigé de B5 et B6 (et B10) |
| Lambda | comme l'universelle, mais paramètres 6 et 7 = `(cos, sin)λ / P_q` | variante de `kepu.f`, corrigée de B5, B6 et B10 |

`λ` est la longitude vraie à `T0`. Elle est définie pour toutes les coniques, et `s` s'en déduit sans inversion numérique (`NU_TO_S`). `P_q = 2π √(q³/μ)`. Idée d'Antoine (2026-11-05), en réponse au constat qu'avec `s·f(e)` le paramètre 7 n'est pas un bon paramètre d'ajustement.

**Aucune de ces corrections n'est dans le dépôt.** Les scripts ci-dessous créent des copies dans `Tests/build_synthetic/variantes/` (non suivi par git).

## Refaire les copies et les vérifier

```
python3 Tests/parametrisation/prepare.py --verifie
```

- `patch_B5B6.py DOSSIER [B10]` : corrige B5 et B6 dans `kepu.f` (et B10 avec le second argument).
- `patch_lambda.py SUB MAIN` : transforme une copie déjà corrigée de B5, B6 et B10 en variante lambda (`kepu.f` et `astrom_univ_mcmco.f`).
- `driver_lambda.f`, `longitude_lambda.f`, `check_lambda.py` : programme de test pour la variante. Il vérifie sur les 11 cas de `Tests/cases` que le χ² du MCMC est égal à celui du LM et du modèle Python, et que le rapport des priors suit le jacobien calculé numériquement (écart résiduel inférieur à 1e-9 après l'écart connu B11).

## Relancer un système synthétique avec une variante

```
V=Tests/build_synthetic/variantes
ORACLE_SYNTH=$PWD/Tests/build_synthetic/essai python3 Tests/synthetic.py lance arc_periastre --rapide --threads 3
ORACLE_SUB=$PWD/$V/lambda/sub ORACLE_MAIN=$PWD/$V/lambda/main ORACLE_SYNTH=$PWD/Tests/build_synthetic/essai_l \
    python3 Tests/synthetic.py lance arc_periastre --univ --rapide --threads 3
ORACLE_SYNTH=$PWD/Tests/build_synthetic/essai_l python3 Tests/synthetic.py analyse arc_periastre --univ
```

Le nombre de pas avant l'inférence se lit dans `out.txt` (ligne `long =` juste avant `Inference now`). Un seul run par version : le hasard du départ compte (gfortran tire une graine différente à chaque exécution).

## Résultats (2026-11-05 et 06, `--rapide`, 2 ou 3 threads)

Systèmes définis dans `Tests/synthetic.py`. Une planète de 3 Mjup autour d'une étoile de 1,5 Msun à 40 pc, astrométrie relative et quelques vitesses radiales de la planète.

| Système | Arc couvert | Classique | Universelle actuelle | Lambda |
|---|---|---|---|---|
| `arc_court` (P = 26 ans, e = 0,55) | 10 ans, 39 % | vérité retrouvée, inférence à 0,9 M pas | retrouvée, 1,1 à 1,2 M pas | retrouvée, 1,1 M pas |
| `arc_tres_court` (P = 47 ans, e = 0,75) | 6 ans, 13 %, autour de l'apoastre | retrouvée, 7,8 M pas | jamais d'inférence en 69 à 78 M pas (9 h) | jamais d'inférence en 86 M pas (9 h) |
| `arc_periastre` (P = 206 ans, e = 0,95, q = 2 au) | 4 ans autour du périastre | retrouvée, 1,1 M pas | LM en échec, puis MCMC coincé à χ² ≈ 72 000 **et déclaré convergé** | retrouvée, 0,3 M pas |
| `arc_parabole` (P = 816 ans, e = 0,98, q = 2 au) | 2 ans autour du périastre | 6,8 M pas, `e` coupé à 0,99 | non lancé | retrouvée, 0,7 M pas ; `e` de 0,94 à 1,04 (95 %), 36 % d'orbites hyperboliques |

Sur les quatre arcs, les versions qui convergent donnent des distributions quasi identiques quand `e` reste loin de 1.

### Ce qu'on en retient

- **Loin du périastre, l'universel n'aide pas.** Sur `arc_tres_court`, les données fixent la distance à l'apoastre `a(1+e)` (6 %) et `a` (9 %), mais `q` à 46 % et `P_q` à 84 %, avec une corrélation de -0,97 entre `ln q` et `e`. Les paramètres universels sont alors les plus mal déterminés. `s·f(e)` a en plus une distribution à queue très lourde (aplatissement 52).
- **Autour du périastre, lambda est meilleure** : 4 à 10 fois moins de pas, et elle seule décrit `e ≥ 1`. La classique ne peut pas dépasser `e = 1` et coupe la distribution.
- **La paramétrisation actuelle `s·f(e)` est à remplacer.** Elle échoue au périastre et peut se déclarer convergée sur une mauvaise orbite : Gelman-Rubin ne voit rien si toutes les chaînes sont coincées au même endroit. Toujours regarder le χ².
- Par pas, l'universel est environ 3 fois plus lent que le classique (fonctions de Stumpff). Non optimisé.
- Décision d'Antoine (2026-11-06) : utiliser lambda. À faire dans le vrai code, après avis d'Hervé et avec les corrections B5, B6 et B10.

### Limites

- Un seul run par version et par système.
- `arc_tres_court` : les versions universelles n'ont pas été menées jusqu'au bout, et la raison exacte du mauvais mélange des chaînes n'est pas établie (pas d'état des chaînes enregistré).
- Les orbites vraies sont toujours elliptiques : le modèle Python (`Tests/model.py`) ne génère pas d'orbites hyperboliques.
