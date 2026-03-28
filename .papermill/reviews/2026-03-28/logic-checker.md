# Logic Checker Report

**Date**: 2026-03-28
**Paper**: "A Binary Threshold for Component Identification in k-out-of-m Systems"

## Proposition 1 (Permutation Invariance)

### Assessment: SOUND, with minor presentation issues

**Claim**: For any k-out-of-m system with independent component lifetimes, f_sys(t; lambda_pi) = f_sys(t; lambda) for all t > 0 and any permutation pi.

**Proof strategy**: The proof argues that (1) the collection of minimal path sets of a k-out-of-m system is invariant under component relabeling (since all subsets of size m-k+1 appear), (2) therefore the critical states for component pi(j) under lambda_pi biject with the critical states for j under lambda, and (3) the substitution j' = pi^{-1}(j) converts the sum under lambda_pi into the sum under lambda.

**Verdict**: The proof is correct. The key insight -- that k-out-of-m systems have all subsets of a given size as path/cut sets, making the collection permutation-invariant -- is valid. The substitution argument in Appendix A is correctly carried out.

### Issues

1. **Minor (logical clarity)**: The proof in Section 3.1 (inline) and Appendix A (full) overlap substantially. The inline proof is a sketch that says essentially the same thing as the appendix. Consider either making the inline proof shorter (just the key insight about path set invariance) or removing the appendix.

2. **Minor (scope of Remark 1)**: Remark 1 correctly notes the result extends to any coherent system whose path set collection is permutation-invariant. This is accurate -- bridge systems with asymmetric component roles would not satisfy the condition.

## Corollary 1 (m! Equivalent Modes)

### Assessment: SOUND

The corollary follows directly from Proposition 1: if f_sys(t; lambda) = f_sys(t; lambda_pi) for all t, then the log-likelihood sum is identical at lambda and lambda_pi, so all permuted vectors achieve the same objective value.

## Claim: Symmetric Point is a "Saddle Point"

### Assessment: UNSUPPORTED

**Location**: Section 3.2, line 378
**Quoted text**: "gradient-based optimizers typically converge to the symmetric saddle point hat{lambda}_1 = ... = hat{lambda}_m"

**Problem**: The paper calls the equal-rates point a "saddle point" without proof. The equal-rates point could be:
- A local maximum (if the Hessian is negative definite)
- A saddle point (if the Hessian has mixed signs)
- A local minimum in some directions

The Monte Carlo results show convergence to this point, but convergence of gradient-based optimizers is consistent with both a local maximum and a saddle point. The claim is stated as fact but is actually an unproven assertion about the geometry of the likelihood surface.

**Suggestion**: Either prove the saddle-point claim (by analyzing the Hessian at the symmetric point) or weaken the language to "the symmetric point" or "the equal-rates point" without asserting saddle-point geometry.

## Claim: Fisher Information Has Rank 1

### Assessment: UNSUPPORTED

**Location**: Section 3.2, lines 392-395
**Quoted text**: "At the symmetric point lambda_1 = ... = lambda_m, the observed Fisher information matrix has rank 1 (in the exponential case): only the sum of lambda_j is locally identifiable."

**Problem**: This is a specific mathematical claim about the rank of a matrix. No proof, derivation, or citation is provided. For a methodology journal like Technometrics, this claim either needs a proof (even a sketch) or a reference.

**Suggestion**: Add a brief derivation or appendix showing the FIM computation at the symmetric point, or cite a reference that establishes this result.

## Likelihood Definitions (Equations 5, 6, 7)

### Assessment: CORRECT with one important caveat

**Equation 5 (Scheme 1 composite likelihood)**: The composite structure is correctly identified as an approximation. The paper appropriately acknowledges this is not a true joint likelihood and that standard errors may be anti-conservative. This is transparent and well-handled.

**Equation 6 (Generalized masking)**: Correctly marginalizes over C(|C_i|, k) candidate failed subsets. The reduction to the standard masked formula at k=1 (series) is verified.

**Equation 7 (Partial autopsy)**: Correctly marginalizes over C(m-r, k_unk) assignments of unknown failures. The boundary cases (r=0 reduces to Scheme 0, r=m reduces to full autopsy) are correctly stated.

**Caveat**: The generalized masking likelihood (Equation 6) implicitly assumes the masking mechanism is non-informative (the probability of observing candidate set C_i does not depend on the parameters theta). The masking model described (each non-failed component independently included with probability p_mask) satisfies this assumption, but this should be stated explicitly as it is a standard condition in the masked data literature.

## Logical Chain: Binary Threshold Claim

### Assessment: EMPIRICALLY SUPPORTED but OVERCLAIMED as "proven"

The abstract says "We prove that a binary identification threshold governs this problem." The paper does NOT prove this -- it demonstrates it empirically via Monte Carlo. Proposition 1 proves the symmetry exists, but the "binary threshold" (that r=1 captures 91% of improvement) is an empirical observation from a specific experimental configuration (m=4, k=2, exponential, n=300, R=50).

**Suggestion**: Change "prove" to "demonstrate" or "show empirically" in the abstract. The symmetry result is proven; the binary threshold is observed.

## Summary

| Finding | Severity | Confidence |
|---------|----------|------------|
| "Saddle point" claim unsupported | Major | High |
| FIM rank-1 claim unsupported | Major | High |
| "Prove" binary threshold overclaim | Major | High |
| Non-informative masking assumption unstated | Minor | High |
| Redundant proof (inline + appendix) | Minor | High |
