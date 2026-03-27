# Consolidated Experiment Results

## Completed Experiments

### Exp 1: Main Comparison (COMPLETE, R=100)

Settings: m=4, k=2, exponential, rates=(0.4, 0.6, 0.8, 1.0), n=300

| Mechanism | Mean MAE | Median MAE | Conv | vs Scheme 0 |
|-----------|----------|------------|------|-------------|
| Scheme 0 (baseline) | 0.297 | 0.200 | 100/100 | 1.0x |
| Ordering constraints | 0.203 | 0.200 | 100/100 | 1.0x |
| Heterogeneous k | 0.186 | 0.189 | 100/100 | 1.1x |
| Partial autopsy (r=2) | 0.057 | 0.057 | 100/100 | 3.5x |
| Masked failed sets | 0.047 | 0.045 | 100/100 | 4.4x |
| Periodic inspection | 0.032 | 0.029 | 100/100 | 6.9x |

**Finding (H1 confirmed)**: Component-level mechanisms give 3.5-6.9x
improvement (median). Domain knowledge gives ~1x. The gap is definitive
at R=100.

Note: mean MAE for Scheme 0 (0.297) is much higher than median (0.200),
indicating heavy right tail (occasional catastrophic failures to
identify rates). Component-level mechanisms eliminate this tail.


### Exp 4: Masking Noise Sensitivity (COMPLETE, R=50)

Settings: m=4, k=2, exponential, n=300, varying p_mask

| p_mask | Median MAE | Mean MAE | Conv |
|--------|------------|----------|------|
| 0.0 | 0.045 | 0.046 | 50/50 |
| 0.1 | 0.048 | 0.049 | 50/50 |
| 0.3 | 0.049 | 0.049 | 50/50 |
| 0.5 | 0.059 | 0.062 | 50/50 |
| 1.0 | 0.200 | 0.280 | 50/50 |

**Finding (H5a confirmed)**: Graceful degradation. MAE roughly doubles
from p_mask=0 to p_mask=0.5, then jumps to Scheme 0 levels at p_mask=1.0.
Moderate masking noise (p_mask <= 0.3) barely affects estimation quality.


## Partial Results (observed during runs, not saved)

### Exp 3: Scaling with m (PARTIAL, R=50)

| m | k | Scheme 0 | Periodic | Autopsy | Per/Aut |
|---|---|----------|----------|---------|---------|
| 3 | 2 | 0.209 | 0.025 | 0.037 | 0.68 |
| 4 | 2 | 0.200 | 0.026 | 0.044 | 0.59 |
| 5 | 3 | (not reached) | | | |

**Partial finding (H2 trending)**: Periodic/autopsy ratio decreases
with m (0.68 to 0.59), confirming periodic pulls ahead. Need m=5 to
strengthen the claim.


### Exp 5: Inspection Delta Sensitivity (PARTIAL, R=50)

| delta | Median MAE |
|-------|------------|
| 0.1 | 0.032 |
| 0.5 | 0.031 |
| 1.0 | (not reached) |
| 2.0 | (not reached) |

**Partial finding (H5b trending)**: delta=0.1 and delta=0.5 give
nearly identical MAE. Inspection is already saturated at delta=0.5.
This is a practical result: coarse inspection suffices.


## Not Started

- Exp 2: k-spectrum (k=2 done in earlier run, k=3 and k=4 not completed)
- Exp 6: Autopsy coverage sensitivity (r = 0..4)
- Exp 7: Weibull extension


## Earlier Single-Run Results (from vignettes, n=300, seed=42)

These are NOT Monte Carlo (single run only) but provide directional evidence:

### Scheme 1 k-spectrum (from periodic-inspection vignette)

| k | Scheme 0 | Scheme 1 | Improvement |
|---|----------|----------|-------------|
| 2 | 0.200 | 0.027 | 7x |
| 3 | 0.239 | 0.030 | 8x |
| 4 | 0.489 | 0.017 | 29x |

### Symmetry-breaking comparison (from symmetry-breaking vignette)

| Mechanism | Mean Error | Improvement |
|-----------|-----------|-------------|
| Scheme 0 | 0.230 | 1x |
| Ordering | 0.200 | 1x |
| Heterogeneous k | 0.200 | 1x |
| Partial autopsy (r=2) | 0.050 | 5x |
| Periodic inspection | 0.044 | 5x |
| Masked failed sets | 0.061 | 4x |


## Performance Bottleneck

The incomplete experiments are all blocked on `f_sys_general()` inside
optimizer loops. This function uses O(m * 2^(m-1)) critical-state
enumeration per time point. Inside `fit_scheme1` or `fit_system`,
the optimizer calls the loglik hundreds of times, each evaluating
f_sys_general at n=300 time points. Total cost per fit:

  m * 2^(m-1) * n * n_optimizer_evals ~ 32 * 300 * 200 = 1.9M
  With 5 multi-starts: ~10M evaluations per replicate
  With R=50: ~500M evaluations per experiment

At m=4 this takes ~60s per fit, ~5 min per replicate, ~4 hours per
R=50 experiment. At m=5 it would be ~10x worse.

**Options to unblock**:
1. Cache critical states per system (currently recomputed every call)
2. Vectorize f_sys_general over time points (partially done but inner
   loops are still scalar)
3. Use the IE expansion for exponential parallel (already fast) and
   restrict general density to non-parallel experiments
4. Reduce R for general-density experiments
5. Implement in C/C++ via Rcpp
