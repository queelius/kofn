# kofn, CRAN Submission Audit

**Audit date:** 2026-04-13
**Package version:** 0.2.0
**Status:** NOT READY. Blocking issues found.

---

## Executive Summary

`kofn` is very close to CRAN-ready but has blocking issues that prevent
submission today:

1. **Three vignettes fail to build** due to stale references to the old
   `multistart_mle()` API (`$par`, `$convergence`). The solver was refactored
   to `compositional.mle`, but vignette code was not updated alongside tests
   and package code.
2. **Stress test runtime is 7m37s**, which will exceed CRAN's 10-minute
   per-check budget once vignette builds and examples are added.
3. **Stale NEWS.md entry** claims an export (`multistart_mle`) that no
   longer exists.
4. **NOTE** from R CMD check about `coef` not being imported from `stats`.
5. **Hidden directory** (`.playwright-mcp/`) leaks into the tarball.
6. **Title case** violation per CRAN convention.

All are mechanical fixes. No design or policy issues. Estimated time to
fix: 30 to 60 minutes.

**Final R CMD check status (without vignette build):** 0 errors,
2 WARNINGs (both vignette-index, will resolve when BLOCKER-1 fixed),
5 NOTEs.

---

## Severity Scale

- **BLOCKER**: CRAN auto-rejects or a human reviewer will bounce it back.
- **WARN**: Reviewer is likely to ask for changes.
- **NOTE**: Cosmetic or borderline; may or may not be raised.
- **OK**: Passes cleanly.

---

## Findings

### BLOCKER-1: Vignette build failures (API drift)

**Severity:** Blocker. CRAN requires vignettes to build cleanly.

**Affected vignettes:**
- `vignettes/exponential-parallel.Rmd` at lines 200 to 215 (chunk `fit-model`)
- `vignettes/kofn-systems.Rmd` at lines 279 to 310 (chunk `monte-carlo-comparison`),
  also lines 267-268, 352, 361, 433, 435
- `vignettes/periodic-inspection.Rmd` at lines 203 to 247 (chunk `comparison-single`),
  also lines 124, 126, 191, 193, 276, 289, 370

**Root cause:** These vignettes reference the old `multistart_mle()` result
shape (`res$par`, `res$convergence == 0`). The solver refactor (commit
`3872b36`, "Replace multistart_mle with compositional.mle solvers") updated
package code and tests but missed the vignettes, which now dispatch on a
`fisher_mle`/`mle_numerical` object that uses `coef(res)` and `res$converged`.

**Errors observed:**
```
exponential-parallel.Rmd: Error in log(), non-numeric argument (BIC.default)
kofn-systems.Rmd:        Error in boxplot.default(), invalid first argument
periodic-inspection.Rmd: Error in round(), non-numeric argument
```

**Fix:** Replace across the three vignettes:
- `res$par` becomes `coef(res)`  (or `unname(coef(res))` for arithmetic)
- `res$convergence == 0` becomes `res$converged` (note: logical, not integer)
- `mle0$par`, `mle1$par`, `mle_s0$par`, `mle_s1$par` become `coef(...)`

**Quick audit after fix:** `R CMD build /home/spinoza/github/rlang/kofn` (with
vignettes) must succeed. The three chunks above should produce HTML output.

**Not broken but flag for sanity check:**
- `vignettes/symmetry-breaking.Rmd:79` uses `res_ord$par`. This is a
  `stats::optim()` result, so `$par` is correct and should stay.
- `vignettes/weibull-em.Rmd:372-378` uses `sc$par`. This is a scenario-config
  list element, not a solver output. Leave it. But line 379 (`res$converged`)
  is correct: solver output uses the new API.

---

### BLOCKER-2: Stale NEWS.md entry

**Severity:** Blocker. Claims an export that no longer exists.

**Location:** `NEWS.md:13`

**Current:**
```
* Exported `multistart_mle()` for users building custom likelihoods.
```

**Issue:** `multistart_mle` was removed in favor of delegating to
`compositional.mle` (NAMESPACE no longer exports it). NEWS entry is factually
wrong and will confuse readers browsing the changelog after install.

**Fix:** Replace with something like:
```
* Migrated optimization to the `compositional.mle` package. The removed
  internal helper `multistart_mle()` is replaced by a private `solve_mle()`
  wrapper around `compositional.mle::lbfgsb()` with Nelder-Mead fallback.
```

---

### WARN-1: `coef` not imported from stats

**Severity:** NOTE from R CMD check, but CRAN maintainers reliably flag this
one on incoming submissions.

**Location:** `R/wei_kofn.R:517`

```r
pp <- parse_params(coef(result), m, "weibull")
```

**Check output:**
```
* checking R code for possible problems ... NOTE
mle_solver : <anonymous>: no visible global function definition for 'coef'
Undefined global functions or variables:
  coef
Consider adding
  importFrom("stats", "coef")
to your NAMESPACE file.
```

**Fix (pick one):**
1. Add `@importFrom stats coef` to a roxygen block in `R/wei_kofn.R`, then
   `devtools::document()`. Preferred since other stats imports already use
   this pattern (`@importFrom stats runif` in `R/observe.R`,
   `@importFrom stats pgamma` in `R/truncated_moments.R`).
2. Or qualify the call: `stats::coef(result)`.

---

### WARN-2: Hidden `.playwright-mcp/` in source tree

**Severity:** R CMD check NOTE. CRAN will query what it is.

**Check output:**
```
* checking for hidden files and directories ... NOTE
Found the following hidden files and directories:
  .playwright-mcp
These were most likely included in error.
```

**Cause:** The Playwright MCP browser-automation tool wrote session console
logs to `.playwright-mcp/` during the Technometrics submission. The directory
is untracked by git but sits in the source tree and was picked up by
`R CMD build`.

**Fix:**
1. Remove the directory: `rm -rf .playwright-mcp/`
2. Add defense-in-depth by appending to `.Rbuildignore`:
   ```
   ^\.playwright-mcp$
   ```
3. Also add to `.gitignore` so it doesn't get committed in future sessions.

---

### WARN-3: Title case violation

**Severity:** CRAN maintainers reliably ask for this change.

**Location:** `DESCRIPTION:2`

**Current:**
```
Title: Component Lifetime Estimation from k-out-of-n System Data
```

**CRAN requires:**
```
Title: Component Lifetime Estimation from k-Out-of-n System Data
```

(Or drop the phrase and reword; e.g., `Component Lifetime Estimation from
k-of-n System Failure Data`.)

Note: CRAN title case treats `k` as lowercase (variable name) but capitalizes
`Out` and `Of` in compound words is an awkward case. An accepted alternative
that dodges the compound problem:
```
Title: Maximum Likelihood Estimation for k-out-of-n System Data
```
(CRAN has accepted `k-out-of-n` lowercase when it is clearly a mathematical
term, but you will need to argue it in the cover message.)

---

### BLOCKER-3: Stress tests not wrapped in `skip_on_cran()`

**Severity:** Blocker. Confirmed test runtime exceeds CRAN budget.

**Locations:**
- `tests/testthat/test-stress.R`, 26 `test_that` blocks. Header comment says
  "more expensive than the unit tests (~2-5 min total)".
- `tests/testthat/test-stress2.R`, 38 `test_that` blocks.

**Confirmed runtime:** `R CMD check` reports `checking tests ... [457s/457s] OK`,
i.e., **7 minutes 37 seconds** on this machine. CRAN's 10-minute per-check
budget covers tests AND vignette rebuild AND examples AND check overhead. On
slower CRAN builders (Debian, Windows, macOS), test runtime alone will likely
exceed 10 minutes and the package will be rejected.

**Fix:** Add at the top of each stress test file:
```r
skip_on_cran()
```
or wrap each `test_that()` with `skip_on_cran()` inside.

Keep them runnable via `devtools::test()` locally and on r-universe CI.

---

### NOTE-2: VignetteBuilder field but no prebuilt index

**Severity:** NOTE from R CMD check. Only appears because I built with
`--no-build-vignettes`. Will go away once BLOCKER-1 is fixed and a normal
`R CMD build` completes.

---

### NOTE-3: `inst/precomputed/paper/` in installed tree

**Severity:** Cosmetic. CRAN tolerates `inst/` files for reproducibility but
reviewers sometimes ask.

**Details:**
- 40 KB of `.rds` files for paper experiment results
- The scripts that generate them (`inst/experiments/`) are in `.Rbuildignore`,
  so the installed package has the outputs but not the generating code

**Options:**
1. **Leave as-is** and add a `inst/precomputed/paper/README.md` explaining that
   these are frozen outputs for paper reproducibility, with a link to the
   GitHub `inst/experiments/` scripts (or to the paper repo). Recommended.
2. Move them out of the package entirely (into the paper repo at
   `~/github/papers/binary-threshold-component-identification-k-out-of-m-systems/`)
   and delete `inst/precomputed/paper/` from the package.

Option 2 is cleaner now that the paper has its own repo; the RDS files are
only consumed by paper figure generation, not by the package API.

---

## OK Findings (passes cleanly)

- All 5 dependencies are on CRAN (compositional.mle, likelihood.model,
  numDeriv, generics, maskedcauses).
- DESCRIPTION Authors@R well-formed with ORCID.
- LICENSE file uses standard MIT template.
- URL and BugReports fields populated.
- R >= 4.1.0 dependency stated (covers native pipe and lambda syntax if used).
- No `setwd()`, `Sys.setenv()`, or `.GlobalEnv` mutation in package code.
- No `:::` internal-function access to other packages.
- No `install.packages()` or `update.packages()` in package code.
- All `T`/`F` matches are math notation in comments; code uses `TRUE`/`FALSE`.
- 69 of 70 `.Rd` files have `\value` (reexports.Rd is correctly internal).
- No `\dontrun{}` misuse (5 files use `\dontrun`/`\donttest`; likely legitimate for long-running fits).
- Tarball size 844K, installed 1.2M, well under 5 MB.
- `cat()`/`print()` calls in code are confined to S3 `print.*` methods.
- `message()` in `wei_kofn.R:410` is gated behind `if (verbose)` (default off).
- `.Rbuildignore` correctly excludes `paper/`, `.claude/`, `.papermill/`,
  `inst/notes/`, `inst/experiments/`.

---

## Recommended Fix Order

```
1. Fix vignette API drift (BLOCKER-1)              [30 min]
   exponential-parallel.Rmd, kofn-systems.Rmd, periodic-inspection.Rmd
   Replace $par with coef(), $convergence == 0 with $converged

2. Rewrite NEWS.md 0.2.0 entry (BLOCKER-2)         [5 min]
   Remove the "Exported multistart_mle" line

3. Add importFrom(stats, coef) (WARN-1)            [2 min]
   @importFrom stats coef in any R file
   devtools::document()

4. Remove .playwright-mcp/ and add to ignore (WARN-2) [2 min]
   rm -rf .playwright-mcp/
   echo '^\\.playwright-mcp$' >> .Rbuildignore
   echo '.playwright-mcp/' >> .gitignore

5. Fix Title case (WARN-3)                         [2 min]
   Edit DESCRIPTION Title field

6. Add skip_on_cran() to stress tests (BLOCKER-3)  [5 min]
   One line at top of each stress test file
   Confirmed needed: tests run 457s currently

7. Decide fate of inst/precomputed/paper/ (NOTE-3) [user judgment]
   Add README or move to paper repo

8. Re-run R CMD check --as-cran with vignettes     [5 min]
   Target: 0 errors, 0 warnings, at most 1 NOTE (New submission)
```

Total: about 1 hour of mechanical work, excluding the paper-repo decision.

---

## Submission Note

After fixes, bump version to **0.2.1** (or prepare **0.3.0** if these fixes
are bundled with anything else). CRAN prefers to see a version change on
re-submission.

---

*Audit produced by the cran-audit skill.*
