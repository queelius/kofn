# Prose Auditor Report

**Date**: 2026-03-28
**Paper**: "A Binary Threshold for Component Identification in k-out-of-m Systems"

## Overall Writing Quality: GOOD

The paper is clearly written with a logical progression from problem statement through theory, methods, experiments, and discussion. The narrative arc -- "symmetry is the barrier, component data breaks it, and surprisingly a single observation suffices" -- is compelling. The prose is generally crisp and avoids unnecessary jargon.

## Structural Issues

### Issue 1: Undefined Hypotheses (CRITICAL)

**Location**: Sections 5.2, 5.3.1, 5.3.2, 5.3.3, and Section 8
**Problem**: The paper references Hypotheses H1, H3, H4, H5a, and H5b by name, confirming or rejecting them:
- "The results confirm Hypothesis H1 decisively" (line 730)
- "This confirms Hypothesis H5a" (line 783)
- "Hypothesis H5b predicted monotonically increasing MAE" (line 804)
- "This confirms Hypothesis H4" (line 837)
- "The rejection of Hypothesis H3" (line 1042)

But these hypotheses are **never defined** in the paper. They appear to come from an experimental plan (.papermill/experiments.md) that was used during drafting but never integrated into the paper text.

A reader encountering "confirms Hypothesis H1" has no idea what H1 states.

**Suggestion**: Either (a) add a numbered hypothesis list in Section 5.1 before the experiments, or (b) remove the hypothesis labels and state the claims directly: "The results confirm that component-level mechanisms provide substantial improvement over system-only data" instead of "The results confirm Hypothesis H1."

### Issue 2: No Figures in an 18-Page Paper (MAJOR)

**Location**: Entire paper
**Problem**: The paper has zero figures. For a Technometrics article, this is highly unusual. Figures would dramatically improve communication of:
- The likelihood surface with its m! modes (a 2D slice)
- The MAE comparison across mechanisms (a bar chart or boxplot)
- The sensitivity curves (MAE vs. p_mask, MAE vs. delta, MAE vs. r)
- The phase transition at r=1 (a line plot showing the sharp drop)

The state file notes "4 figures, 3 tables planned" but only 3 tables materialized. The sensitivity analysis results are presented as inline tabular environments rather than proper tables, which makes them hard to find and reference.

**Suggestion**: Add at least 3 figures: (1) boxplot comparing the 6 mechanisms, (2) MAE vs. r showing the phase transition, (3) MAE vs. delta showing insensitivity. These are the paper's selling points and deserve visual emphasis.

### Issue 3: k-Spectrum and Scaling Results Missing (MAJOR)

**Location**: Section 5.4 and Section 8
**Problem**: The conclusion discusses "The rejection of Hypothesis H3 (improvement ratio increasing with k)" but the k-spectrum results are never presented. They are deferred to "kofn package vignettes." Similarly, the scaling experiment (exp3) is not shown. The precomputed data for both experiments exists.

If the conclusion discusses a rejected hypothesis, the evidence must appear in the paper.

**Suggestion**: Add a table or figure showing the k-spectrum results (exp2) and scaling results (exp3). Even a brief table would suffice.

### Issue 4: Sensitivity Tables Are Not Proper Floats (MINOR)

**Location**: Sections 5.3.1, 5.3.2, 5.3.3
**Problem**: The three sensitivity tables are embedded as `center` environments rather than `table` environments. They have no captions, no labels, and cannot be cross-referenced. This is inconsistent with the proper Table 1-3 formatting used elsewhere.

**Suggestion**: Convert these to proper `table` floats with captions and labels (e.g., Table 4, 5, 6).

## Narrative Issues

### Issue 5: "Immediately Eliminates Half" Is Inaccurate (MINOR)

**Location**: Section 5.3.3, line 847
**Quoted text**: "knowing whether one named component failed or survived at T_sys immediately eliminates half of the permutation equivalences"
**Problem**: The paper then says "reduces the ambiguity from 24 permutations to 6." But 24 to 6 is a 75% reduction, not "half." The word "half" is misleading.

**Suggestion**: Say "eliminates most of the permutation equivalences" or give the accurate fraction.

### Issue 6: Inconsistent Use of k-out-of-m vs. k-out-of-n (MINOR)

**Location**: Throughout
**Problem**: The paper uses k-out-of-m consistently, but the standard literature uses k-out-of-n. This is not wrong (m and n are both common), but the CLAUDE.md file uses k-out-of-n while the paper uses k-out-of-m. The state file mixes both. For consistency with the broader reliability literature, k-out-of-n is more standard, though k-out-of-m avoids collision with the sample size n.

**Suggestion**: The paper's choice of m is defensible (to avoid confusion with sample size n). Add a brief note explaining this notational choice.

### Issue 7: Remark 2 Is Buried (MODERATE)

**Location**: Lines 476-484
**Problem**: Remark 2 states that the periodic inspection analysis is restricted to k >= 2 because the composite likelihood breaks down for series systems. This is an important limitation that deserves more prominence. A reader interested in series systems might not notice this remark.

**Suggestion**: Elevate this to the opening of Section 4.1 or add a note in Section 5.1 when describing the experimental design.

## Notation Consistency

- Notation is well-defined in Table 1 and used consistently throughout
- The use of theta for the general parameter vector and lambda for the exponential rate is clean
- The bold conventions (btheta, blambda, bx) are consistent

## Section-by-Section Flow

1. **Introduction**: Clear problem statement, well-motivated by the power supply example, clean roadmap
2. **Background**: Efficient coverage of necessary material. Definition 1 is well-stated
3. **Symmetry problem**: Proposition 1 and its consequences are clearly presented
4. **Mechanisms**: Each mechanism is defined with a formal likelihood and an intuitive "Why it works" explanation. Good structure
5. **Monte Carlo**: Results are clearly reported but suffer from the hypothesis-labeling problem and missing figures
6. **Discussion**: Good practical recommendations (the decision tree). Scope and limitations are honestly discussed
7. **Software**: Appropriate for a methodology paper
8. **Conclusion**: Clean but references the unreported H3 result

## Summary

| Finding | Severity | Confidence |
|---------|----------|------------|
| Undefined hypotheses H1/H3/H4/H5a/H5b | Critical | High |
| No figures | Major | High |
| k-spectrum results discussed but not shown | Major | High |
| Remark 2 (series breakdown) insufficiently prominent | Moderate | High |
| Inline sensitivity tables not proper floats | Minor | High |
| "Eliminates half" inaccurate (actually 75%) | Minor | High |
| k-out-of-m vs k-out-of-n convention note missing | Minor | Low |
