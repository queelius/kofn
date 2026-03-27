# Papermill State

## Meta
- **Stage**: experiments-designed
- **Format**: LaTeX (article class)
- **Target venue**: Technometrics or Lifetime Data Analysis
- **Repository**: /home/spinoza/github/rlang/kofn
- **Initialized**: 2026-03-25
- **Last updated**: 2026-03-26

## Authors
1. **Alexander Towell**:Southern Illinois University Edwardsville
   - Email: lex@metafunctor.com
   - ORCID: 0000-0001-6443-9897
   - Role: Sole author
   - Affiliation: PhD candidate, Computer Science, SIUE

## Thesis

Individual component lifetime estimation from k-out-of-n system data
is fundamentally limited by **permutation symmetry** in the system
density: f_sys(t; lambda_1, ..., lambda_m) is invariant under
relabeling of components, creating m! equivalent modes on the
likelihood surface.

We show that:
1. System-level data alone (Scheme 0) cannot identify individual
   component parameters at any k.
2. Domain knowledge (ordering constraints, heterogeneous fleet
   structures) provides negligible improvement (~1x).
3. **Component-level information** breaks the symmetry:
   - Periodic inspection (temporal intervals): 7-29x improvement
   - Partial autopsy (binary failed/survived for r of m components):
     5x improvement even with r = ceil(m/2)
   - Masked failed sets (which k components failed): 4x improvement
4. Partial autopsy is the most *plausible* mechanism (standard
   maintenance practice) and generalizes both masked cause-of-failure
   (series, k=1) and full autopsy (r=m).
5. Periodic inspection races ahead of autopsy as m increases because
   it provides O(m) temporal constraints vs O(k) binary constraints.

**One-sentence thesis**: Breaking the permutation symmetry inherent in
k-out-of-n system densities requires component-level information, and
even minimal component inspection (checking 2 of 4 components after
failure) delivers dramatic improvement over system-only data.

## Key Results (from package vignettes and experiments)

| Mechanism | mean_err | Improvement | Plausibility |
|-----------|----------|-------------|-------------|
| Scheme 0 (baseline) | 0.230 | 1x | Always |
| Ordering constraints | 0.200 | 1x | High |
| Heterogeneous k | 0.200 | 1x | Moderate |
| Partial autopsy (r=2/4) | 0.050 | 5x | Very high |
| Periodic inspection | 0.044 | 5x | High |
| Masked failed sets | 0.061 | 4x | Moderate |

Scheme 1 (periodic inspection) across k-spectrum (m=4):
| k | Scheme 0 | Scheme 1 | Improvement |
|---|----------|----------|-------------|
| 2 | 0.200 | 0.027 | 7x |
| 3 | 0.239 | 0.030 | 8x |
| 4 | 0.489 | 0.017 | 29x |

Scheme 1 fails for k=1 (series) due to composite likelihood
breakdown:surviving components' intervals should be conditioned
on survival past T_sys.

## Novel Contributions

1. **Symmetry-breaking framing**:reframing k-out-of-n estimation
   as a symmetry problem rather than just a censoring problem
   (Novelty: HIGH, per prior art search)
2. **Generalized masking for k-out-of-n**:candidate failed set
   likelihood marginalizing over C(|C|, k) subsets, extending
   masked cause-of-failure from series to general k
   (Novelty: MODERATE-HIGH)
3. **Partial autopsy**:inspecting r < m components without C1
   guarantee, generalizing Flehinger et al. 1998 two-stage model
   (Novelty: HIGH)
4. **Comparative analysis across k-spectrum**:systematic comparison
   of 6 symmetry-breaking mechanisms with plausibility assessment
   (Novelty: HIGH)
5. **Quantitative improvement factors**:5-29x magnitudes,
   minimal component inspection sufficiency
   (Novelty: MODERATE)

## Prior Art (key references)

- Nowik (1990):identifiability conditions, parallel non-identifiability
- Meilijson (1981):autopsy model foundation
- Komarova (2017):k-out-of-n nonparametric identification (autopsy)
- Rodrigues, Pereira & Polpo (2019):masked data for coherent systems
- Sarhan & El-Bassiouny (2003):parallel system masked data (2 comp)
- Flehinger, Reiser & Yashchin (1998):two-stage diagnosis (closest to partial autopsy)
- Miyakawa (1984), Lin/Usher/Guess (1993):foundational series masked MLE
- Craiu & Duchesne (2004):EM for masked competing risks
- Qiang & Pena (2025):shrinkage, component vs system data
- Samaniego (2007), Coolen et al. (2014):signatures for component inference

## Software

The `kofn` R package (https://github.com/queelius/kofn) implements
all methods described in the paper:
- `kofn()` constructor, `loglik()`, `fit()`, `rdata()`:S3 model API
- `loglik_masked()`, `rdata_masked()`:generalized masking
- `loglik_scheme1()`, `fit_scheme1()`:periodic inspection
- `observe_*()` functor family:composable observation schemes
- `coherent_system()`, `kofn_system()`:system infrastructure
- 307 tests, 93% coverage, R CMD check clean

## Outline

See `.papermill/outline.md` for the full section-by-section outline.

8 sections, ~20 pages for Technometrics:
1. Introduction (2pp): problem, thesis, roadmap
2. Background (3pp): coherent systems, censoring, observation schemes, likelihood
3. The symmetry problem (2.5pp): invariance proposition, practical consequences, taxonomy
4. Symmetry-breaking mechanisms (5pp): periodic inspection, generalized masking, partial autopsy, domain knowledge
5. Monte Carlo comparison (4pp): R=100 head-to-head, scaling with m, sensitivity, Weibull
6. Discussion (2pp): decision tree, scope, connections
7. Software (1pp): kofn package
8. Conclusion (0.5pp)

4 figures, 3 tables planned.

## Experiments

See `.papermill/experiments.md` for the full plan with hypotheses,
variables, and success criteria.

7 experiments, ~3 hours total:
- [ ] Exp 1: Main comparison, 6 mechanisms, R=100 (Fig 2, Tab 3)
- [ ] Exp 2: k-spectrum, Scheme 0 vs periodic across k=2,3,4
- [ ] Exp 3: Scaling with m=3,4,5,7, periodic vs autopsy (Fig 3)
- [ ] Exp 4: Sensitivity to masking noise p_mask (Fig 4a)
- [ ] Exp 5: Sensitivity to inspection delta (Fig 4b)
- [ ] Exp 6: Sensitivity to autopsy coverage r (Fig 4c)
- [ ] Exp 7: Weibull extension

## Decisions

- Composite likelihood for Scheme 1 is documented, not corrected
  (full joint likelihood requires integrating over constrained region)
- Scheme 1 excluded for k=1 (series) due to composite breakdown
- Focus on exponential for main results; Weibull as extension
