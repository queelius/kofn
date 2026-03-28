# Citation Verifier Report

**Date**: 2026-03-28
**Paper**: "A Binary Threshold for Component Identification in k-out-of-m Systems"

## Bibliography Integrity

### Entries: 35 total in paper.bib
### Actually cited in text: 19
### Uncited entries: 16

## Critical Finding: 16 Uncited Bibliography Entries

The following entries appear in paper.bib but are NEVER cited in the manuscript:

1. **BhattacharyaSamaniego2010** - "Estimating Component Characteristics from System Failure-Time Data" (Naval Research Logistics). DIRECTLY relevant to the paper's topic. Should be cited and compared against.
2. **CramerNavarro2015** - "Progressive Type-II Censoring and Coherent Systems" (Naval Research Logistics). Relevant to censoring in coherent systems.
3. **DembinskaJasinski2021** - "MLE Based on Discrete Component Lifetimes of a k-out-of-n System" (TEST). Directly relevant to k-out-of-n estimation.
4. **DempsterLairdRubin1977** - The foundational EM paper. Should be cited when discussing EM-based approaches.
5. **FlehingerReiserYashchin1998** - "Parametric Modeling for Survival with Competing Risks and Masked Failure Causes" (Lifetime Data Analysis). The state file identifies this as the closest prior work to the partial autopsy mechanism. Not citing it is a significant omission.
6. **KunduBasu2000** - "Analysis of Incomplete Data in Presence of Competing Risks" (JSPI). Relevant to the competing risks discussion.
7. **MeekerEscobar1998** - The first edition of the Meeker/Escobar text. The 2nd edition (2022) is cited instead, which is fine, but the 1998 entry should be removed from the .bib.
8. **NgNavarroBalakrishnan2017** - "SEM Algorithm for System Lifetime Data with Known Signature" (Computational Statistics). Directly relevant to signature-based component estimation.
9. **QiangPena2025** - "Shrinkage Estimation for Component Reliability from System and Component Test Data" (Technometrics). A RECENT paper in the TARGET VENUE on the same topic. Not citing this is a serious omission.
10. **RodriguesPereira2019** - "Masked Data Analysis in Coherent Systems" (JSCS). Directly relevant.
11. **Samaniego2007** - "System Signatures and their Applications" (Springer). Foundational for the signature approach to component estimation. Should be discussed.
12. **SarhanElBassiouny2003** - "Estimation of Components Reliability in a Parallel System Using Masked System Life Data." Directly relevant to parallel system estimation.
13. **Sun2006** - "Statistical Analysis of Interval-Censored Failure Time Data." Relevant to the periodic inspection mechanism.
14. **Usher1996** - "Weibull Component Reliability-Prediction in the Presence of Masked Data." Relevant to the Weibull extension.
15. **Weibull1951** - The foundational Weibull distribution paper. Could be cited when introducing Weibull components.
16. **YangNgBalakrishnan2019** - "EM Algorithm for System-Based Lifetime Data with Unknown System Structure." Directly relevant.

## Citation Accuracy

The 19 cited references appear to be used accurately in context:
- Tsiatis (1975) and Cox (1959) are correctly cited for competing risks identifiability
- Miyakawa (1984) and Lin et al. (1993) are correctly cited for masked cause-of-failure
- Nowik (1990) is correctly cited for identifiability in coherent systems
- Komarova (2017) is correctly cited for nonparametric k-out-of-n identification
- Meilijson (1981) is correctly cited for autopsy-based identifiability
- Varin et al. (2011) is correctly cited for composite likelihood
- Stephens (2000) and Jasra et al. (2005) are correctly cited for label switching

## Missing References Not in Bibliography

1. **Navarro & Rychlik (2007)** or related work on exchangeable components in coherent systems
2. **Lindsay (1988)** - foundational composite likelihood paper, would strengthen Section 4.1
3. **Navarro (2022)** - modern textbook on system reliability theory
4. Any recent work (2020-2025) on Bayesian component estimation from system data

## Formatting Issues

- All .bib entries use consistent formatting
- DOI is provided only for Tsiatis (1975). Should be added for all entries where available.
- The numDeriv entry has "note" field with version info that is out of date (2016.8-1.1)
- CoolenCoolen-MareschalCoolen-Schrijner2014: The hyphen in the author field may cause BibTeX issues (appears to work in compilation, but the key is awkward)

## Summary

| Finding | Severity | Confidence |
|---------|----------|------------|
| 16 uncited bibliography entries | Critical | High |
| Qiang & Pena (2025) -- same venue, same topic, not cited | Major | High |
| Flehinger et al. (1998) -- closest prior art for partial autopsy, not cited | Major | High |
| Samaniego (2007), Bhattacharya & Samaniego (2010) -- signature approach not discussed | Major | High |
| DOIs missing for most entries | Minor | High |
| numDeriv version outdated | Minor | Low |
