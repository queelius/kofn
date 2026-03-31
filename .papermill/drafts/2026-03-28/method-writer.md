# Method Writer Output: Section 4, Section 7, Appendix B

## Section 4: Symmetry-Breaking Mechanisms

### 4.1 Periodic Inspection
- Composite log-likelihood: log f_sys(t) + sum_j log[F_j(b_j) - F_j(a_j)]
- Documented as composite (not joint) likelihood, citing Varin et al. 2011
- Series breakdown (k=1) explicitly noted as a limitation
- Mechanism: different rates -> different interval patterns

### 4.2 Generalized Masking
- Definition 2: generalized masked likelihood summing over binom(|C|,k) subsets
- Reduces to standard masked cause-of-failure for k=1
- Reduces to system density for k=m (parallel)
- Masking model: p_mask controls noise level
- Mechanism: faster components appear in candidate sets more often

### 4.3 Partial Autopsy
- Definition 3: marginalize over uninspected components
- Spectrum from r=0 (Scheme 0) to r=m (full autopsy)
- Cost: O(binom(m-r, k_unk) * k) per observation
- Key insight: first inspection breaks most of the symmetry

### 4.4 Domain Knowledge
- Ordering constraints via cumulative increment reparameterization
- Heterogeneous fleet: pooling likelihoods from different k values
- Both shown to be ineffective (results in Section 5)

### Table 2: Summary of all 6 mechanisms

## Section 7: Software
- kofn package architecture
- Key functions listed
- Quality metrics: 307 tests, 93% coverage
- Ecosystem connections

## Appendix B: Computational Details
- IE expansion formula with O(2^m) terms
- Closed-form integral for interval-censored observations
- Partial autopsy enumeration algorithm (Algorithm 1 pseudocode)
