# Results Writer Output: Sections 5-6

## Section 5: Monte Carlo Comparison

### 5.1 Experimental Design
- Baseline: m=4, k=2, exponential, rates=(0.4,0.6,0.8,1.0), n=300, R=100
- MAE on sorted estimates as primary metric (eliminates label ambiguity)
- Multi-start optimization: 5 restarts, L-BFGS-B primary, Nelder-Mead fallback

### 5.2 Main Comparison (Table 3)
From exp1_main_comparison.rds (R=100):
| Mechanism | Median MAE | Mean MAE | SD | Improvement |
|-----------|-----------|---------|------|-------------|
| Scheme 0 | 0.200 | 0.297 | 0.187 | 1.0x |
| Ordering | 0.200 | 0.203 | 0.031 | 1.0x |
| Heterogeneous | 0.189 | 0.186 | 0.066 | 1.1x |
| Partial autopsy | 0.057 | 0.057 | 0.023 | 3.5x |
| Masking | 0.045 | 0.047 | 0.020 | 4.4x |
| Periodic | 0.029 | 0.032 | 0.014 | 6.9x |

H1 confirmed: 3.5-6.9x for component-level, ~1x for domain knowledge.

### 5.3 k-Spectrum
From exp2_k_spectrum.rds (R=100 per k):
- k=2: 7.0x, k=3: 6.3x, k=4: 4.4x
- H3 REJECTED: improvement DECREASES with k
- Explanation: Scheme 0 baseline improves with k; periodic MAE stable

### 5.4 Scaling with m
From exp3_scaling.rds (R=50):
- Per/Aut ratio: 0.58-0.68, periodic consistently better
- H2 confirmed: periodic O(m) constraints vs autopsy O(k)

### 5.5 Sensitivity Analyses
- Masking noise (exp4): graceful degradation, robust to p_mask <= 0.3
- Delta sensitivity (exp5): COMPLETELY INSENSITIVE, 0.031-0.033 across 20x range
- Autopsy coverage (exp6): phase transition at r=1, captures 69% of full benefit

### 5.6 Weibull Extension
From exp7_weibull.rds (R=50):
- Same ranking, larger magnitudes: periodic 9.1x, masked 4.8x

## Section 6: Discussion
- Decision tree for practitioners
- Scope: series systems use maskedcauses, composite likelihood limitation
- Connections: label switching, competing risks identifiability, optimal design
