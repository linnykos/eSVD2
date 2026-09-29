# Brainstorming: eSVD2 on the way to CRAN
**Date:** 2026-09-29
**Author:** Kevin (with Claude Code)
**Goal:** Index of brainstorms for this project, and short go/no-go notes. Long brainstorms live in `additional_context/` and are linked from here.

## 2026-09-29: the nuisance rate and the inflated p-values in `devel`

Full document: `additional_context/OVERDISPERSION_BRAINSTORM.md` (eleven
ideas, with tables). Scripts: `additional_context/overdispersion_brainstorm/`.

| # | Idea | Go/No-Go |
|---|---|---|
| 1 | Unit-free cap, rate <= 10 x the gene's median library size | **Go. Chosen by Kevin 2026-09-29; not yet implemented in `R/`.** Report: `additional_context/overdispersion_brainstorm/overdispersion_cap_claude.html` |
| 2 | Legacy cap, rate <= the gene's largest library size | Go; 78% to 100% of genes at the cap |
| 3 | Empirical-Bayes shrinkage of the log rate | Same as 1 at 600 cells; twice the cap's false discoveries at 3000 cells with widely spread rates. Only with a cap |
| 4 | Lower end of the profile-likelihood interval | As 3; no-go at 50% |
| 5 | One rate for all genes | No-go: fails under `generate_null()` |
| 6 | Record capped and boundary genes | Go as a complement; no-go alone |
| 7 | Median in `bool_stabilize_underdispersion` | No-go alone |
| 8 | Welch variance between individuals only | No-go |
| 9 | Rescale the statistic within strata of the rate | No-go |
| 10 | A variance that accounts for the fit | Not run; half-day check described |
| 11 | One real data set | Not run; needs a `PAPER_DATA` path |
