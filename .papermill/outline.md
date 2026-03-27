# Paper Outline: Breaking Permutation Symmetry in k-out-of-n Estimation

**Target venue**: Technometrics (~20 pages, methodology + computation + application)
**Narrative arc**: Problem → Barrier → Framework → Solutions → Evidence → Recommendations

---

## 1. Introduction (~2 pages)

**Purpose**: Motivate the problem and state the thesis.

**Key arguments**:
- k-out-of-n systems are ubiquitous (redundant electronics, load-sharing structures, network routing)
- Estimating individual component parameters from system data is a foundational reliability problem
- The standard framing (censoring) is incomplete. The real barrier is *permutation symmetry*: the system density is invariant under relabeling of component parameters
- This paper provides a unified framework for understanding what information breaks the symmetry, and compares six mechanisms quantitatively

**Opens with**: A concrete example. A 4-component parallel power supply: observing system failure times tells you the system is unreliable, but not which component to replace. You need component-level information.

**Ends with**: Paper roadmap and one-sentence thesis.

**No figures or tables.**

---

## 2. Background and Notation (~3 pages)

**Purpose**: Establish the mathematical framework readers need.

### 2.1 Coherent systems and k-out-of-n structure

- Minimal path sets, minimal cut sets, structure function phi
- k-out-of-n as special case: k failures for system failure
- System lifetime T_sys as the k-th order statistic of component lifetimes
- Series (k=1), parallel (k=m), intermediate k

### 2.2 Censoring structure

- At T_sys: 1 component exact, k-1 left-censored, m-k right-censored
- The censoring *pattern* depends on k but the *symmetry* does not

### 2.3 Observation schemes

- Scheme 0: system lifetime only
- Periodic inspection: component failure times interval-censored
- Autopsy: component status (failed/survived) observed post-failure
- Connection to masked cause-of-failure literature for series (k=1)

### 2.4 Likelihood framework

- System density f_sys(t; theta) via critical-state enumeration
- Exponential: inclusion-exclusion expansion for parallel
- Weibull: numerical density via general coherent system engine

**Table 1**: Notation summary.

---

## 3. The Permutation Symmetry Problem (~2.5 pages)

**Purpose**: The core theoretical section. Establish the barrier formally.

### 3.1 Invariance of the system density

- **Proposition 1**: For a k-out-of-n system with i.i.d. components or exchangeable parameters, f_sys(t; lambda_pi) = f_sys(t; lambda) for any permutation pi.
- Proof sketch: path sets are invariant under relabeling
- Consequence: m! equivalent modes on the likelihood surface
- Connection to Nowik (1990) identifiability conditions

### 3.2 Practical consequences

- MLE convergence to symmetric point (all rates equal)
- Fisher information matrix singular at the symmetric point
- Identifiable functions (sum of rates for series) vs non-identifiable individual parameters

### 3.3 What breaks symmetry?

- Taxonomy: information about *which* component (identity) vs *when* each component failed (temporal)
- Domain knowledge (ordering, fleet structure) vs data-derived information
- Preview of the six mechanisms

**Figure 1**: Likelihood surface for m=2 parallel system showing permutation symmetry (contour plot with two equivalent modes at (lambda_1, lambda_2) and (lambda_2, lambda_1)).

---

## 4. Symmetry-Breaking Mechanisms (~5 pages)

**Purpose**: Present each mechanism with its likelihood and data requirements.

### 4.1 Periodic inspection (composite likelihood)

- Data: T_sys exact, each T_j in interval [a_j, b_j)
- Composite log-likelihood: log f_sys(t) + sum_j log[F_j(b_j) - F_j(a_j)]
- Why it works: different rates produce different inspection patterns
- Limitation: composite (not joint) likelihood; anti-conservative SEs
- Limitation: breaks down for k=1 (surviving components' intervals unconditional)
- Extension from parallel-only to general k via general density engine

### 4.2 Generalized masking (candidate failed sets)

- Data: T_sys, candidate set C with |C| >= k, true failed set F subset of C
- **Definition 1**: The generalized masked likelihood:
  L_i(theta) = sum_{F in C(C_i, k)} sum_{j in F} f_j(t_i) prod_{l in F\j} F_l(t_i) prod_{l not in F} S_l(t_i)
- For k=1 (series): reduces to sum_{j in C} w_j(t), the standard maskedcauses formula
- For k=m (parallel): F = {1,...,m} forced, reduces to f_sys(t)
- C1/C2/C3 conditions reinterpreted for general k
- Why it works: faster-failing components appear in candidate sets more often

### 4.3 Partial autopsy

- Data: T_sys, binary status (failed/survived) for r of m inspected components
- **Definition 2**: The partial autopsy likelihood:
  L_i(theta) = sum over consistent failed sets (marginalize over uninspected components)
- Relationship to generalized masking: partial autopsy is a generalization where some component statuses are unknown
- Spectrum: r=0 is Scheme 0, r=m is full masking (exact failed set)
- Why it works: same as masking, but with uncertainty about uninspected components
- Computational cost: O(2^(m-r) * C(m-r, k-k_known) * k) per observation

### 4.4 Domain knowledge mechanisms

- **Ordering constraints**: reparameterize as cumulative increments, optimize in the ascending cone
- **Heterogeneous fleet structures**: pool likelihoods from systems at different k values
- Why they fail: ordering eliminates equivalent modes but the remaining mode is still broad; heterogeneous k provides complementary censoring but the density is still permutation-invariant within each k

**Table 2**: Summary of mechanisms with data requirements, likelihood formula, and computational cost.

---

## 5. Monte Carlo Comparison (~4 pages)

**Purpose**: Head-to-head comparison on a common baseline with proper Monte Carlo.

### 5.1 Experimental design

- Baseline: m=4, exponential, rates (0.4, 0.6, 0.8, 1.0), n=300, R=100 replicates
- Primary metric: mean absolute error of sorted rate estimates
- Secondary: sum-of-rates bias, convergence rate, computation time
- k values: 2 (primary), also 3 and 4 for k-spectrum analysis

### 5.2 Main comparison (k=2, m=4)

- All 6 mechanisms on the same baseline
- Box plots of sorted rate errors across R=100 replicates

**Figure 2**: Box plots of mean absolute error for each mechanism (main result figure).

**Table 3**: Summary statistics (median MAE, IQR, convergence rate) for each mechanism.

### 5.3 Periodic inspection vs autopsy: scaling with m

- m = 3, 4, 5, 7 at k = ceil(m/2)
- Periodic inspection uses delta = 0.5
- Partial autopsy uses r = ceil(m/2)
- Does periodic inspection pull ahead as m grows?

**Figure 3**: MAE vs m for periodic inspection, full autopsy, and partial autopsy.

### 5.4 Sensitivity analyses

- Masking noise: p_mask = 0, 0.1, 0.3, 0.5, 1.0
- Inspection granularity: delta = 0.1, 0.5, 1.0, 2.0
- Autopsy coverage: r = 1, 2, ..., m

**Figure 4**: Sensitivity curves (MAE vs p_mask, MAE vs delta, MAE vs r).

### 5.5 Weibull extension

- Same baseline but with Weibull components (shape = 1.5)
- Verify that the symmetry-breaking ranking holds for non-exponential

---

## 6. Discussion (~2 pages)

**Purpose**: Practical recommendations and scope.

### 6.1 When to use which mechanism

- Decision tree based on available data:
  1. Can you inspect components during operation? -> Periodic inspection
  2. Can you inspect some components after failure? -> Partial autopsy
  3. Do you know which components failed? -> Generalized masking
  4. Do you know the component reliability ordering? -> Ordering constraints (marginal help)
  5. None of the above -> Scheme 0 (accept symmetry limitation)

### 6.2 Scope boundaries

- Series (k=1): use maskedcauses (candidate cause, not candidate failed set)
- k >= 2: kofn with inspection or masking
- Composite likelihood limitation for Scheme 1 (not the full joint)
- Computational limits: general density engine is O(2^m), partial autopsy is O(2^(m-r))

### 6.3 Connections to related problems

- Label switching in Bayesian mixture models
- Identifiability in competing risks
- Optimal experimental design for reliability systems

---

## 7. Software (~1 page)

**Purpose**: Brief description of the kofn R package.

- Package architecture: kofn model object, closure-returning generics, observation functors
- Key functions: kofn(), loglik_masked(), rdata_masked(), fit_scheme1(), observe_*()
- Availability: CRAN (pending), GitHub, r-universe
- Reproducibility: all experiments in the paper are reproducible via package vignettes

---

## 8. Conclusion (~0.5 pages)

- Restate thesis in light of evidence
- The practical takeaway: even minimal component inspection (2 of 4) breaks symmetry
- Open question: optimal inspection design (which components to inspect, how many)

---

## Appendices

### A. Proofs

- Proposition 1 (permutation invariance)
- Derivation of generalized masked likelihood

### B. Computational details

- Inclusion-exclusion expansion for exponential parallel
- Critical-state enumeration for general coherent systems
- Partial autopsy enumeration algorithm

---

## Figures and Tables Summary

| ID | Type | Content | Section |
|----|------|---------|---------|
| Fig 1 | Contour | Likelihood surface showing permutation symmetry (m=2) | 3 |
| Fig 2 | Box plot | Main comparison: MAE across 6 mechanisms (R=100) | 5.2 |
| Fig 3 | Line plot | Periodic vs autopsy scaling with m | 5.3 |
| Fig 4 | Line plots | Sensitivity: p_mask, delta, r | 5.4 |
| Tab 1 | Notation | Symbol definitions | 2 |
| Tab 2 | Summary | Mechanisms with data/likelihood/cost | 4 |
| Tab 3 | Results | Main comparison statistics | 5.2 |

**Estimated total**: ~20 pages (appropriate for Technometrics)
