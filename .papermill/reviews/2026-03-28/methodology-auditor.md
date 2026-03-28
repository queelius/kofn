# Methodology Auditor Report

**Date**: 2026-03-28
**Paper**: "A Binary Threshold for Component Identification in k-out-of-m Systems"

## Monte Carlo Experimental Design

### Configuration
- m=4, k=2, exponential, lambda=(0.4, 0.6, 0.8, 1.0), n=300
- Main comparison: R=100
- Sensitivity analyses: R=50
- Weibull extension: R=50

### Issue 1: R=50 Is Too Small for the Claims Made (MAJOR)

**Problem**: The sensitivity analyses use only R=50 replicates. Key claims rest on these R=50 experiments:
- The "91% of improvement at r=1" comes from the R=50 autopsy coverage experiment
- The "delta-insensitivity" comes from the R=50 inspection granularity experiment
- The Weibull extension uses R=50

With R=50 and the variability observed (e.g., SD=0.023 for partial autopsy r=2 in the R=100 experiment), the standard error of the median MAE is roughly SD/sqrt(R) ~ 0.003. This means differences of 0.003 or smaller are within noise. The r=3 to r=4 anomaly (median MAE increases from 0.042 to 0.049) is 0.007, which is at the edge of noise -- but the paper dismisses it as noise without a formal test.

**Suggestion**: Increase R to at least 200 for the sensitivity analyses. For a Technometrics paper, R=1000 would be standard for a pure Monte Carlo study. Alternatively, provide confidence intervals or standard errors for all summary statistics.

### Issue 2: Abstract Misrepresents R (MAJOR)

**Problem**: The abstract says "A Monte Carlo comparison (R=100 replicates, m=4 exponential components, n=300 observations) yields two central findings." Both "central findings" (91% threshold and delta-insensitivity) actually come from R=50 experiments, not R=100. The R=100 experiment is the main comparison (Table 3) which does not directly support either central finding.

**Suggestion**: Either run the sensitivity analyses at R=100 or correct the abstract to accurately represent which experiments support which claims.

### Issue 3: Single Configuration (MAJOR)

**Problem**: All main results use m=4, k=2, exponential, lambda=(0.4, 0.6, 0.8, 1.0). The paper acknowledges additional experiments are "available in the kofn package vignettes" (line 861) but defers them. For a methodology paper, the generalizability of the binary threshold claim requires:
- Variation in m (3, 5, 7, 10)
- Variation in k (not just k=2)
- Variation in the degree of heterogeneity (lambda ratios)
- Different families (Weibull is shown briefly but only for R=50)
- Different n (100, 500, 1000)

The precomputed results include exp2 (k-spectrum) and exp3 (scaling), but these are not presented in the paper. They should be.

**Suggestion**: Include the k-spectrum and scaling results as tables or figures. These experiments are already done -- they just need to be presented.

### Issue 4: MAE Is the Only Metric (MODERATE)

**Problem**: The paper uses only MAE of sorted estimates. This single metric obscures important information:
- **Bias**: Is the estimator biased? The sorted-MAE cannot distinguish bias from variance.
- **RMSE**: Standard metric for Monte Carlo studies.
- **Coverage probability**: Are confidence intervals from the observed information well-calibrated? This is especially important given the composite likelihood concern.
- **Component-wise error**: Which components are hardest to estimate? The aggregate MAE hides this.
- **Convergence quality**: "Conv. 100/100" is reported but what does convergence mean? Gradient norm threshold? Log-likelihood stabilization?

**Suggestion**: Add RMSE, bias, and ideally coverage probability. Component-wise error decomposition (a table showing MAE by component rank) would be informative.

### Issue 5: No Error Bars or Statistical Tests (MODERATE)

**Problem**: Tables 3 and the sensitivity tables report point estimates (median, mean, SD, IQR) but no standard errors or confidence intervals for these summary statistics. When claiming "identical accuracy" for delta=0.1 vs. delta=2.0, a formal test (or at least a confidence interval for the difference) is needed.

**Suggestion**: Add bootstrap confidence intervals for the reported medians, or use a Wilcoxon test to formally compare mechanisms.

### Issue 6: The r=3 to r=4 Anomaly (MINOR)

**Problem**: Median MAE increases from r=3 (0.042) to r=4 (0.049). The paper says "the mean MAE continues to decrease (0.043 -> 0.046)" but 0.043 -> 0.046 is an INCREASE, not a decrease. This is a factual error.

**Verified against precomputed data**: r=3 mean MAE = 0.04341, r=4 mean MAE = 0.04622. Both median and mean increase from r=3 to r=4.

**Suggestion**: Fix the factual error. The increase is likely Monte Carlo noise (R=50), but the statement "continues to decrease" is false. Either note that both median and mean increase slightly (within sampling error), or increase R to resolve the anomaly.

### Issue 7: Multi-Start Optimization (MINOR)

**Problem**: Only 5 random restarts are used. For a 4-dimensional exponential model with 24 modes (m!=4!=24), 5 restarts may not adequately explore the mode space. The convergence rate of 100/100 is reassuring but does not guarantee the global optimum was found in each replicate.

**Suggestion**: Document the sensitivity of results to the number of restarts. Consider using more restarts (10-20) for the component-level mechanisms where modes may be fewer but sharper.

### Issue 8: No Seed Reporting (MINOR)

**Problem**: Random seeds for reproducibility are not discussed. The precomputed RDS files enable reproducibility at the data level, but the paper should mention this.

## Summary

| Finding | Severity | Confidence |
|---------|----------|------------|
| Abstract misrepresents R (100 vs 50) | Major | High |
| R=50 too small for sensitivity claims | Major | High |
| Single configuration limits generalizability | Major | High |
| MAE as sole metric | Moderate | High |
| No error bars or statistical tests | Moderate | High |
| Factual error: "continues to decrease" | Minor | High |
| 5 restarts may be insufficient | Minor | Moderate |
| No seed reporting | Minor | Low |
