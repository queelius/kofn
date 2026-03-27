# Experiment Plan

## Hypotheses

**H1 (main)**: Component-level information (inspection, autopsy, masking)
provides statistically significant improvement over system-only data for
individual rate estimation, while domain knowledge (ordering, heterogeneous
fleet) does not.

**H2 (scaling)**: Periodic inspection's advantage over partial autopsy grows
with m, because inspection provides O(m) temporal constraints while autopsy
provides O(k) binary constraints.

**H3 (k-spectrum)**: Periodic inspection improvement over Scheme 0 increases
with k (from series to parallel), because parallel systems have the most
severe symmetry and the most room for improvement.

**H4 (minimal sufficiency)**: Inspecting r = 2 out of m = 4 components
(partial autopsy) captures most of the improvement from full autopsy (r = m),
because the first few inspected components provide the most information about
which components failed.

**H5 (robustness)**: The ranking of mechanisms is robust to (a) masking
noise (p_mask > 0), (b) inspection granularity (delta > 0), and (c) the
component distribution family (exponential vs Weibull).


## Variables

### Independent (controlled)

| Variable | Values | Rationale |
|----------|--------|-----------|
| m (components) | 3, 4, 5, 7 | Tests scaling; m=7 is near computational limit |
| k (failures) | 2, ceil(m/2), m | Series-adjacent, intermediate, parallel |
| n (sample size) | 300 | Large enough for stable MLE, small enough for fast runs |
| R (replicates) | 100 | Standard for Monte Carlo in reliability; gives SE of MAE ~ MAE/10 |
| Family | exponential, Weibull(shape=1.5) | Memoryless baseline + IFR extension |
| delta (inspection) | 0.1, 0.5, 1.0, 2.0 | Fine to coarse; median component lifetime ~ 1.5 |
| p_mask (masking noise) | 0, 0.1, 0.3, 0.5, 1.0 | Exact to fully masked |
| r (autopsy coverage) | 1, 2, ..., m | Full spectrum from minimal to complete |

### Dependent (measured)

| Metric | Definition | Primary? |
|--------|-----------|----------|
| MAE | mean(abs(sort(est) - sort(true))) | Yes |
| sum_bias | mean(sum(est)) - sum(true) | Secondary |
| sum_RMSE | sqrt(mean((sum(est) - sum(true))^2)) | Secondary |
| convergence_rate | fraction of R replicates that converge | Secondary |
| wall_time | total seconds for R replicates | Reporting only |

### Fixed

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| True rates (m=4) | (0.4, 0.6, 0.8, 1.0) | Spread factor 2.5x; not too extreme |
| True rates (m=3) | (0.4, 0.7, 1.0) | Same spread |
| True rates (m=5) | (0.4, 0.6, 0.8, 1.0, 1.2) | Same spacing |
| True rates (m=7) | (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) | Uniform spacing |
| n_starts (optimizer) | 5 | Balance speed vs mode-finding |
| Weibull shape | 1.5 (all components) | Mild IFR; shape is NOT the target parameter |
| Weibull scales (m=4) | (2.5, 1.7, 1.25, 1.0) | Chosen so mean lifetime ~ 1/rate for comparability |
| seed | 2026 (base) + replicate index | Reproducible |


## Experiments

### Experiment 1: Main Comparison (Section 5.2, Figure 2, Table 3)

**Tests**: H1
**Design**: m=4, k=2, n=300, R=100, exponential, rates=(0.4, 0.6, 0.8, 1.0)

6 mechanisms:
1. Scheme 0 (baseline)
2. Ordering constraints (cumulative increment reparameterization)
3. Heterogeneous k (pool n/2 from k=2 and n/2 from k=3)
4. Partial autopsy (r=2 of m=4)
5. Periodic inspection (delta=0.5)
6. Masked failed sets (exact, p_mask=0)

**Output**: Box plots of MAE (Fig 2), summary table (Tab 3).

**Success criterion for H1**: Component-level mechanisms (4,5,6) have
median MAE < 50% of Scheme 0 median MAE. Domain knowledge mechanisms
(2,3) have median MAE > 80% of Scheme 0 median MAE.

**Estimated time**: ~30 min (dominated by partial autopsy, 2^2 enumeration
per obs per eval, 100 reps x 5 starts).


### Experiment 2: k-Spectrum (Section 5.2 extension)

**Tests**: H3
**Design**: m=4, k in {2, 3, 4}, n=300, R=100, exponential

For each k, compare Scheme 0 vs Periodic inspection (delta=0.5).

**Output**: Table of MAE by k and scheme; improvement ratio.

**Success criterion for H3**: Improvement ratio (Scheme 0 MAE / Scheme 1 MAE)
increases monotonically with k.

**Estimated time**: ~15 min.


### Experiment 3: Scaling with m (Section 5.3, Figure 3)

**Tests**: H2
**Design**: m in {3, 4, 5}, k = ceil(m/2), n=300, R=50, exponential

3 mechanisms at each m:
- Scheme 0 (baseline)
- Periodic inspection (delta=0.5)
- Full autopsy (masked, p_mask=0)
- Partial autopsy (r = ceil(m/2))

Skip m=7 for partial autopsy (2^(7-4)=8 enumeration per obs, too slow
for R=50). Include m=7 for Scheme 0 and periodic only.

**Output**: Line plot MAE vs m (Fig 3).

**Success criterion for H2**: At m=5, periodic MAE / autopsy MAE < 1.0
(periodic wins). At m=3, the ratio > 0.8 (roughly equal).

**Estimated time**: ~45 min.


### Experiment 4: Sensitivity to Masking Noise (Section 5.4)

**Tests**: H5(a)
**Design**: m=4, k=2, n=300, R=50, exponential, p_mask in {0, 0.1, 0.3, 0.5, 1.0}

Compare loglik_masked with varying masking noise.

**Output**: MAE vs p_mask curve (part of Fig 4).

**Success criterion**: MAE degrades gracefully (monotonically) with p_mask.
At p_mask=1.0 (full masking), MAE should approach Scheme 0 baseline.

**Estimated time**: ~15 min.


### Experiment 5: Sensitivity to Inspection Granularity (Section 5.4)

**Tests**: H5(b)
**Design**: m=4, k=2, n=300, R=50, exponential, delta in {0.1, 0.5, 1.0, 2.0}

Compare fit_scheme1 with varying delta.

**Output**: MAE vs delta curve (part of Fig 4).

**Success criterion**: MAE increases monotonically with delta (coarser
inspection = less information). At small delta, MAE approaches Scheme 2
(complete data) performance.

**Estimated time**: ~20 min.


### Experiment 6: Sensitivity to Autopsy Coverage (Section 5.4)

**Tests**: H4
**Design**: m=4, k=2, n=300, R=50, exponential, r in {0, 1, 2, 3, 4}

r=0 is Scheme 0, r=4 is full autopsy.

**Output**: MAE vs r curve (part of Fig 4).

**Success criterion for H4**: Over 50% of the improvement from r=0 to r=4
is captured by r=2. Diminishing returns as r increases.

**Estimated time**: ~30 min (r=1 and r=2 involve 2^3 and 2^2 enumeration).


### Experiment 7: Weibull Extension (Section 5.5)

**Tests**: H5(c)
**Design**: m=4, k=2, n=300, R=50, Weibull(shape=1.5), scales derived from
exponential rates for comparability.

3 mechanisms: Scheme 0, periodic (delta=0.5), masked (p_mask=0).

**Output**: Table comparing exponential vs Weibull rankings.

**Success criterion**: Same ranking (periodic > masked > Scheme 0) holds
for Weibull. Improvement magnitudes may differ.

**Estimated time**: ~30 min (general density engine for Weibull is slower).


## Reproducibility

- All experiments use `set.seed(2026 + replicate_index)` for reproducibility
- Results saved to `inst/precomputed/paper/` as RDS files
- Each experiment has a standalone R script in `inst/experiments/`
- The `kofn` package version used is recorded in each output file
- Total estimated wall time: ~3 hours on a single core


## Computational Constraints

- General density engine: O(m * 2^(m-1)) per time point. At m=7, this is
  448 evaluations per observation per loglik call. With n=300 and optimizer
  needing ~100 loglik evals, that is ~13M evaluations per fit. Feasible
  but slow (~60s per fit). With R=50 and 5 starts, that is ~4 hours for
  m=7 alone. Restrict m=7 to Scheme 0 and periodic only.

- Partial autopsy enumeration: O(2^(m-r)) per observation. At m=7, r=4,
  this is 2^3=8 subsets per obs, manageable. At m=7, r=2, this is
  2^5=32 subsets per obs, still feasible. At m=7, r=1, this is 2^6=64,
  getting slow with R=50.

- Total computation budget: ~3 hours. Parallelize across experiments
  where possible (experiments are independent).
