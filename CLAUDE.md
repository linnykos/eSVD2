# CLAUDE.md — eSVD2 (CRAN submission)

## Workflow Instructions
1. **Always enter plan mode** before starting any non-trivial task.
2. **Use superpower skills** where relevant: `/code-review` for code changes.
3. **Always invoke `/r-style-guide`** before writing or editing any `.R` / `.Rmd` file in this repo.
4. **After every prompt**, run `/project-state`: refresh the current-state sections of the individual contributor's `CLAUDE_[name].md` in place, and append a dated entry to their `HISTORY_[name].md`. Write only non-obvious things; skip anything already in the code or git history.

## Project Context (High Level)
**Paper/Project**: `eSVD2` is the R package implementing eSVD-DE, a cohort-wide
differential-expression method for single-cell RNA-seq. The package is published
(Lin, Qiu & Roeder, *BMC Bioinformatics* 2024, 25:113,
[doi:10.1186/s12859-024-05724-7](https://doi.org/10.1186/s12859-024-05724-7)) and
currently installs from GitHub only.

**Authors**: Kevin Z. Lin (aut, cre — method, R/C++ code, analyses), Yixuan Qiu
(R/C++ code, especially the `src/` optimization backend), Kathryn Roeder (method).

**Current goal (as of 2026-08-27)**: convert this repository into a *formal R
package suitable for CRAN submission*. The emphasis is **correctness first**
(hidden failure modes, silent `NaN`/`Inf`, untested code paths, dependency
policy), not performance. The working document driving this effort is
`additional_context/CRAN_READINESS.md` — read it before touching package
structure, `DESCRIPTION`, `NAMESPACE`, `src/`, or `tests/`. Its companion is
`additional_context/UNIT_TEST_PLAN.md`, the proposed test suite (~200 tests, each
with an ID and an explicit oracle); read it before writing any test or
regenerating anything in `tests/assets/`. It is a proposal awaiting review, not
an agreed plan.

**Method in one paragraph** (so a session need not re-read the paper):
eSVD-DE models a cells × genes count matrix `A` hierarchically as
`A_ji | λ_ji ~ Poisson(ℓ_ji · λ_ji)` with `λ_ji ~ Gamma(mean = μ_ji, var = γ_j·μ_ji)`,
where `μ_ji = exp(Y_j·ᵀ X_i· + Z_j,(cc)·C_i,(cc))` is a low-rank "predictable"
expression and `ℓ_ji = exp(Z_j,-(cc)ᵀ C_i,-(cc))` is the covariate-adjusted
sequencing depth. The pipeline is: `format_covariates` → `initialize_esvd`
(glmnet Poisson ridge per gene, then SVD of `log1p(A) − ZᵀC`) →
`opt_esvd` (alternating constrained Newton in C++, two phases: case–control
coefficient held fixed, then freed) → `reparameterization_esvd_covariates`
(orthogonalize `X` against `C`, then against itself, for identifiability) →
`estimate_nuisance` (per-gene overdispersion via `gamma_rate` /
`log_gamma_rate` in C++) → `compute_posterior` (Gamma–Poisson posterior mean and
variance per cell × gene) → `compute_test_statistic` (average posteriors within
individual, summarize case and control groups as Gaussian mixtures, Welch
two-sample statistic) → `compute_pvalue` (Welch df → t-to-z transform →
`locfdr` empirical null → BH). `compute_test_per_gene` is a memory-lean
alternative that fuses the last three steps in a per-gene loop.

## Repository Layout
Paths below are relative to this file's directory and are the same for everyone.

| Path | What it is |
|---|---|
| `R/` | Package R source. Pipeline order: `format_covariates.R` → `initialization.R` → `optimization.R`(+`optimization_helper.R`) → `reparameterization.R` → `nuisance.R` → `posterior.R` → `compute_test_statistic.R` → `compute_pvalue.R`(+`multtest.R`). `eSVD.R` is the end-to-end wrapper; `compute_test_per_gene.R` is the low-memory fused path. |
| `src/` | C++ backend (Rcpp + RcppEigen + BH). `data_loader.*` (dense/sparse iteration), `distribution.h`/`family*.cpp` (7 exponential families), `constrained_newton.cpp` (line-searched Newton), `optimization.cpp` (`opt_x`/`opt_yz`), `gamma_rate.cpp` (overdispersion MLE via Boost Newton–Raphson). |
| `inst/include/` | `eSVD2.h`, exposed for `LinkingTo`. |
| `man/`, `NAMESPACE` | roxygen2-generated — never hand-edit. |
| `data/` | Shipped reference gene lists (`gandal_df`, `housekeeping_df`, `sfari_df`, `velmeshev_gene_df`). |
| `tests/testthat/` | testthat suite. Fixtures live in `tests/assets/` and are loaded with `load("../assets/...")`. |
| `vignettes/` | `eSVD2.Rmd` (self-contained toy runs), `asd.Rmd` / `asd-preprocess.Rmd` (need external downloads; heavy chunks are `eval = FALSE`). |
| `additional_context/` | Reference PDFs + `summary.md` index + `CRAN_READINESS.md` (audit) + `UNIT_TEST_PLAN.md` (proposed test suite). **Excluded from the CRAN tarball** via `.Rbuildignore`. |
| `oldcode/` | Superseded code, kept for reference; `.Rbuildignore`d. |

**`.gitignore` vs `.Rbuildignore`.** These are deliberately different.
`.gitignore` keeps build artifacts and session state out of *git history*.
`.Rbuildignore` keeps collaboration material out of the *CRAN tarball* —
`additional_context/`, `CLAUDE*.md`, `HISTORY_*.md`, `brainstorming_*.md`,
`.claude/`, `.githooks/`, `oldcode/`, `docs/` are all tracked on GitHub but must
never ship to CRAN. **When you add a new top-level file or folder, decide
explicitly which of the two it belongs in.**

## External Locations
Folders this project depends on that live **outside** the project root. A path is
external unless it can be written relative to this file's directory without `..`.

**No per-machine filesystem path appears in this file**, because such a path is
only meaningful on one person's one machine. This file names each location and
says what it is for; the path — and which machine it is on — is recorded per
person in `CLAUDE_[name].md` under *External Locations (per-machine paths)*.
Hostnames and URIs that are the same for everyone are fine here.

| Location name | Purpose | Copy semantics |
|---|---|---|
| `EXAMPLES_REPO` | Clone of <https://github.com/linnykos/eSVD2_examples> — all analyses reported in the paper. Not required to develop the package. | per-person copy |
| `PAPER_DATA` | Downloaded public datasets used by the vignettes and the paper (Adams GSE136831, Habermann GSE135893, Smillie, Velmeshev). Large; never tracked in git. | per-person copy |
| `WAS2CODE_REPO` | Clone of the Was2CODE project (Tati collaboration). Source of `R/esvd_helper.R`, the cohort-filtering wrapper being imported into `eSVD2` as §2.17 of `UNIT_TEST_PLAN.md`. Not required to develop the package once the file is imported. | per-person copy |

Git remote (same for everyone): `https://github.com/linnykos/eSVD2.git`.
pkgdown site (same for everyone): <https://linnykos.github.io/eSVD2/>.

Refer to these locations by name in prose, code comments, and session notes.
**Copy semantics matter before you write:** a write to shared storage lands in
every collaborator's view immediately; a write to a per-person copy does not, and
the copies drift.

To resolve a name to a real path, read the current user's `CLAUDE_[name].md`. If
that person has no row for the location, ask them — do not guess, and do not
reuse another collaborator's path.

## Who Is Using This Session?
**Detect the current user** by running: `echo $USER`. This table maps each login to
that person's **first-name** context file; it is the source of truth that
`/brainstorm` and `/project-state` use to resolve `$USER` to the right filename
(so the login `kevinlin` maps to `CLAUDE_kevin.md`, never `CLAUDE_kevinlin.md`).

| Username (login) | Current-state file (first name) | History archive |
|---|---|---|
| `kevinlin` | `CLAUDE_kevin.md` | `HISTORY_kevin.md` |

**File ownership.** Each row above names one person's files, and **only that
person's session writes them.** Once `$USER` resolves to a first name, that is the
only suffix you may create, edit, append to, rename, or delete — every other
collaborator's `CLAUDE_[name].md`, `HISTORY_[name].md`, and
`brainstorming_[name].md` is read-only. Read them for context when useful; never
modify them, not even to fix a typo. This repository is shared, so an edit lands
in the owner's working copy immediately and can overwrite state they wrote from
their own machine. If a collaborator's file looks wrong or stale, say so instead
of editing it. This master `CLAUDE.md` is the exception: it is shared and any
collaborator may update it. The single override is the user, in the current turn,
*directing you to write* that exact file — confirm once, then write it.

**Session startup — run this before any other work.** All four cases below are
normal; none is an error to report back to the user.

1. Run `echo $USER` and look for a matching row.
2. **Row exists and the file exists** → read that person's `CLAUDE_[name].md`
   immediately. Do **not** read `HISTORY_[name].md` at startup; it is the
   append-only session log, consulted only on demand.
3. **Row exists but `CLAUDE_[name].md` does not** → this is a collaborator's first
   session. Initialize `CLAUDE_[name].md` and `HISTORY_[name].md` from their
   templates via `/project-state`, then continue.
4. **No matching row** → ask the user their first name, add a row for them to the
   table above, then initialize their files as in case 3. Never guess a first
   name from the login.

Adding a row is a project-level change, so record it the same way as any other
(`/project-state`).

## Reference Papers
`additional_context/summary.md` indexes the PDFs in `additional_context/`. Read it
before opening any PDF. Open a PDF only when you need exact equations or verbatim
quotes that the summary does not carry.

## CRAN Conventions for This Repo
- New R files drafted by Claude are named with a `_claude` suffix so the human can
  review before integrating.
- `man/` and `NAMESPACE` are generated by roxygen2 (`devtools::document()`).
  Edit the roxygen comments, never the `.Rd`.
- Before claiming any packaging change works, run
  `R CMD build` + `R CMD check --as-cran` and quote the output. `/project-state`
  entries about check status must cite the actual result, not an expectation.

## Post-Prompt Update Instructions
After completing each user prompt, run `/project-state`. It will:
- **Refresh in place** the current-state sections of `CLAUDE_[name].md` (Project
  Status, Key Methodological Details, Open Questions / Next Steps).
- **Append a dated entry at the bottom** of `HISTORY_[name].md` recording new
  decisions, resolved/open questions, non-obvious code rationale, and empirical
  findings.

Do NOT record: things already in the code, git history, or reproducible from code.
