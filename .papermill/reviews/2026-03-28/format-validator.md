# Format Validator Report

**Date**: 2026-03-28
**Paper**: "A Binary Threshold for Component Identification in k-out-of-m Systems"
**Target venue**: Technometrics

## Build Verification

- **Compilation**: pdflatex + bibtex compiles cleanly with zero warnings or errors
- **Output**: 18 pages, 394KB PDF
- **All labels resolve**: No undefined references
- **All citations resolve**: No undefined citations

## Venue Format Compliance

### Technometrics Requirements
Technometrics uses a specific LaTeX style (typically ASA format via `asa.cls` or the Technometrics template). The paper currently uses:
- `\documentclass[12pt]{article}` with 1-inch margins
- `plainnat` bibliography style
- Standard `natbib` citation commands

**Issue (Major)**: The paper does not use the Technometrics/ASA LaTeX template. For submission, it would need to be reformatted. This is normal for a pre-submission draft, but should be noted.

### Line Numbers
Line numbers are set up via `\usepackage{lineno}` but commented out (`% \linenumbers`). For review submission, these should be enabled.

## Label Resolution

### Referenced labels (all resolve): 9 unique labels referenced
- prop:symmetry, sec:periodic, sec:masking, sec:autopsy, sec:monte-carlo
- sec:autopsy-sensitivity, sec:main-comparison, sec:ie-expansion
- rem:series-breakdown, tab:notation, tab:mechanisms, tab:main-results
- alg:autopsy
- Various equations

### Orphan labels (defined but never referenced): 25 labels
Many section and definition labels are defined but never cross-referenced:
- def:autopsy-ll, def:fsys, def:masked-ll
- eq:ie, eq:mae, eq:symmetry, eq:tsys
- sec:autopsy-algorithm, sec:censoring, sec:coherent, sec:connections
- sec:consequences, sec:delta-sensitivity, sec:design, sec:domain
- sec:intro, sec:invariance, sec:likelihood, sec:masking-noise
- sec:obs-schemes, sec:proof-symmetry, sec:recommendations
- sec:scope, sec:sensitivity, sec:taxonomy

This is not a compilation issue but suggests the cross-referencing could be more thorough (e.g., "as defined in Definition 2" rather than just stating the definition).

## Table Formatting

### Proper table floats (3):
1. Table 1 (Notation summary) -- well-formatted with booktabs
2. Table 2 (Mechanisms summary) -- well-formatted
3. Table 3 (Main comparison) -- well-formatted

### Inline tables (not proper floats, 3):
- Section 5.3.1: Masking noise sensitivity (no caption, no label)
- Section 5.3.2: Inspection granularity sensitivity (no caption, no label)
- Section 5.3.3: Autopsy coverage sensitivity (no caption, no label)

These use `center` + `tabular` environments and cannot be cross-referenced. They should be converted to `table` floats.

## Figure Count: 0

No figures in the paper. For an 18-page methodology paper in Technometrics, at least 3-4 figures are expected.

## Algorithm Formatting

Algorithm 1 (Appendix C) is properly formatted using the `algorithm` + `algpseudocode` environments. It is correctly referenced from the text.

## Typography

- microtype is loaded for improved typesetting
- booktabs is used for all table rules (good)
- cleveref with capitalize option for cross-references (good)
- hyperref loaded for PDF links

## Page Budget

At 18 pages with 12pt article class, the paper is on the long side for Technometrics (which typically targets 15-20 pages in their format). Reformatting to the ASA template would change the page count. The lack of figures means the current 18 pages is text-heavy.

## Summary

| Finding | Severity | Confidence |
|---------|----------|------------|
| Not using Technometrics/ASA template | Major | High |
| No figures | Major | High |
| 3 sensitivity tables not proper floats | Minor | High |
| 25 orphan labels | Minor | High |
| Line numbers commented out | Minor | High |
| DOIs missing from most bib entries | Minor | High |
