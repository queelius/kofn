# Writing Plan: Breaking Permutation Symmetry in k-out-of-n Estimation

## Phase Assignment

| Section | Specialist | Dependencies | Key Content |
|---------|-----------|-------------|-------------|
| 2. Background and Notation | formal-writer | none | Coherent systems, censoring structure, observation schemes, likelihood framework, Table 1 |
| 3. The Permutation Symmetry Problem | formal-writer | Section 2 | Proposition 1 (invariance), practical consequences, taxonomy, Figure 1 description |
| 4. Symmetry-Breaking Mechanisms | method-writer | Sections 2-3 | Periodic inspection, generalized masking, partial autopsy, domain knowledge, Table 2 |
| 5. Monte Carlo Comparison | results-writer | Sections 3-4 | All 7 experiments, Figures 2-4, Table 3, hypothesis testing |
| 6. Discussion | results-writer | Section 5 | Decision tree, scope, connections to related problems |
| 7. Software | method-writer | none | kofn package architecture |
| Related Work (embedded in Sections 2-3) | literature-writer | none | Prior art positioning, Nowik, Meilijson, masked data literature |
| 1. Introduction | orchestrator | all sections | Motivating example, thesis, roadmap |
| Abstract | orchestrator | all sections | Self-contained summary |
| 8. Conclusion | orchestrator | all sections | Summary, implications, future work |
| Appendix A. Proofs | formal-writer | Sections 2-3 | Proposition 1 proof, masked likelihood derivation |
| Appendix B. Computational Details | method-writer | Section 4 | IE expansion, critical-state enumeration, partial autopsy algorithm |

## Parallelism Plan

Round 1 (parallel): formal-writer (Secs 2-3, App A), method-writer (Secs 4, 7, App B), literature-writer (related work)
Round 2 (parallel, after Round 1): results-writer (Secs 5-6)
Round 3 (sequential): orchestrator writes Introduction, Conclusion, Abstract

## Experiment Data Summary

All 7 experiments completed with RDS files in inst/precomputed/paper/.
Key numbers verified:
- Exp 1: R=100, 6 mechanisms, periodic median MAE=0.029, scheme0=0.200
- Exp 2: k-spectrum R=100, improvement ratio DECREASES: 7.0x, 6.3x, 4.4x
- Exp 3: Scaling m=3,4,5, periodic/autopsy ratio 0.58-0.67
- Exp 4: Masking noise p_mask, graceful degradation confirmed
- Exp 5: Delta completely insensitive (0.031-0.033 across 0.1 to 2.0)
- Exp 6: Autopsy coverage, r=1 captures 69% improvement (0.200->0.062)
- Exp 7: Weibull same ranking, larger magnitudes (9.1x periodic, 4.8x masked)
