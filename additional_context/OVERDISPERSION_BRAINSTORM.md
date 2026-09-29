# Brainstorming: the overdispersion estimate and the inflated p-values in `devel`

**Date:** 2026-09-29
**Author:** Kevin (with Claude Code)
**Goal:** Decide how `eSVD2` 1.1.0 should estimate (or bound) the per-gene
nuisance rate before CRAN, given that the uncapped `gamma_rate` raises false
discoveries (`version_comparison/version_comparison_claude.Rmd`, Q10 and Q-LFC-1
in `CLAUDE_kevin.md`). This document is for choosing a route. **Nothing in the
package was changed.**

> **Decision (Kevin, 2026-09-29): Idea 1, the cap at 10 times the gene's
> median library size. It is not to be implemented in `R/` until Kevin has
> read the report** `overdispersion_brainstorm/overdispersion_cap_claude.html`,
> which explains the change with formulas, plots and the wiki's citations.
> Section 3.6 below records what that report added, including two
> corrections to this document.

Notation: `β_j` is the Gamma **rate** the package stores in `nuisance_vec`
(large = little overdispersion), `s_ji` is the fitted covariate-adjusted library
size (`library_mat`), `μ_ji` is the fitted mean without the library
(`mean_mat`). Paper citations use the page names of the overdispersion wiki
(`amywatt/git/overdispersion_wiki/wiki/pages/<name>.md`), e.g. `lause-2021`.

---

## 1. Bottom line

1. **The divergence is not a bug in the optimizer.** A gene's rate runs to 1e7
   exactly when its Poisson likelihood is at least as high as any
   negative-binomial one, which happens when its residuals are slightly
   *under*dispersed around the fit (Pearson statistic 0.95 to 0.99). This is
   the classical boundary case of the dispersion MLE (`dai-2012`,
   `piegorsch-1990`, `lause-2021`). No better solver removes it.
2. **`master` was a capped MLE, and the cap did nearly all the work.**
   `pmin(devel rate, max_i s_ji)` reproduces `master`'s rates to a median of
   0.1% and its Welch statistics to within 0.025. In the model-generated
   regimes **78% to 100% of genes sit at that cap** (39% under
   `generate_null()`). The published behavior was close to "every gene gets the
   same unit-free rate".
3. **Better estimation does not fix the test when some genes are nearly
   Poisson in truth.** Feeding the pipeline the *true* rates gives 5.7 false
   discoveries per data set (pooled FDP 0.26) when the true rates span 0.5 to
   200, against about 1 for any capped estimate; 55 of the 57 come from genes
   whose true rate is above 30 times the library size. With true rates below
   about 30 the truth works and is the most powerful choice (section 3.6).
   The null scale of the statistic grows with the rate, so what calibrates
   the test is a **limit on how large a rate may be relative to the others**.
4. **At 600 cells any bound in the right range works about equally well**,
   and recovers `master`'s calibration. **At 3000 cells the hard caps still
   hold, and the two likelihood-based bounds (empirical-Bayes shrinkage, the
   profile lower bound) weaken**: about twice the caps' false discoveries when
   the true rates are spread widely (section 3.4). Power is unchanged or better everywhere except
   `low_count` and `wide_rate`, where the MLE's extra discoveries come with a
   worse ranking of the genes. Differences between the bounded candidates are
   within simulation noise (SE of the false-discovery count is about 0.5 per
   cell of the tables).
5. **The two downstream fixes tried (a between-individual Welch variance, a
   rate-stratified rescaling) both failed** their go/no-go.

**Recommendation.** Idea 1 (a unit-free cap, `β_j ≤ 10 · median_i s_ji`,
applied in R after `gamma_rate`) for the release, with Idea 6 (record and
report which genes were capped or at the boundary) alongside it. Idea 2 (the
legacy cap) is the alternative if reproducing the 2024 paper's numbers exactly
matters more than letting 90% of genes keep their own estimate. Idea 3
(empirical-Bayes shrinkage) performs the same at 600 cells and is the most
defensible in a methods sense, but its protection fades as cells are added
(section 3.4), so on its own it is not a safe default. Idea 11 (one real data
set) should be run before committing to either.

---

## 2. What is going on

All numbers are from the dry-runs of section 3.

### 2.1 `master` is `min(MLE, max library size)`

| Regime | Genes with devel rate above `max_i s_ji` |
|---|---|
| baseline / null / weak_de | 0.91 / 0.91 / 0.92 |
| near_poisson | 1.00 |
| wide_rate | 0.78 |
| low_count | 0.92 |
| generate_null | 0.39 |

This is why `master`'s rates had Spearman about 0 with the truth: for nine
genes in ten the "estimate" is the gene's largest library size.

### 2.2 The "two- to threefold upward bias" is mostly a units mismatch

The generator draws `λ ~ Gamma(mean = exp(nat), rate = β)` with a library size
averaging 1. The fit puts the gene intercept and the covariate effects into
`s_ji`, and `β` shares units with `s` (they enter only through `s + β`). The
same overdispersion is therefore `β_true · s` in the fit's units.

| Regime | Median `β̂/β_true`, raw | Median, unit-free (`β̂ / median_i s_ji`) | Spearman raw → unit-free |
|---|---|---|---|
| baseline | 2.50 | 1.13 | 0.59 → 0.77 |
| null | 2.43 | 1.14 | 0.60 → 0.77 |
| minimal_design | 2.34 | 1.04 | 0.61 → 0.79 |
| wide_rate | 2.04 | 1.03 | 0.85 → 0.88 |
| low_count | 0.32 | 1.39 | 0.49 → 0.64 |

So away from the boundary the MLE is a good estimator (4% to 16% high, 39% at
low counts). The remaining excess is the direction expected when the fit
absorbs part of the residual variance (`nuisance-parameter-bias`). **The
statement in `CLAUDE_kevin.md` and in the comparison report that devel
"overestimates two- to threefold" should be corrected.**

### 2.3 Divergence is the boundary of the likelihood

| Regime | Genes at the boundary | Mean Pearson statistic, boundary genes | Interior genes |
|---|---|---|---|
| baseline / null / weak_de | 0.5% / 0.6% / 0.4% | 0.98 / 0.99 / 0.97 | 1.25 |
| low_count (mean count 0.1 to 0.6) | 2.5% | 0.95 | 1.18 |
| wide_rate (true rate 0.5 to 200) | 16% | 0.95 | 1.35 |
| near_poisson (true rate 50 to 1000) | 48% | 0.95 | 1.04 |

"At the boundary" means the Poisson log-likelihood is at least that of the
returned rate. Every such gene has a rate above 1e4 and vice versa. In
`near_poisson` even the interior genes gain a median of only 0.26
log-likelihood units over Poisson, so whether a gene diverges there is a coin
flip unrelated to its true rate. `lause-2021` reports the same figure (49.9% of
`theta.ml` estimates diverging at constant `θ`), and `underdispersion` makes
the reading explicit: a boundary estimate means *uninformative*, not *Poisson*.

Real UMI data are mostly low-count, so the 2.5% of `low_count` is a better
guide to real data than the 0.5% of `baseline`. Neither is a measurement on
real data (Idea 11).

### 2.4 Why a large rate inflates the statistic

The Welch denominator is the mixture variance `mean_i(V̄_i) + var_i(m_i)`, and
the first term (the per-cell posterior variance `(A + μβ)/(s + β)²`, not divided
by the number of cells) is 95% to 98% of it. That term behaves like `μ/β` for
large `β`, so the denominator vanishes while the numerator (a difference of
fitted means, which contain the fitted case-control coefficient) does not.

Two consequences follow.

* **The Gaussianized statistic is far from N(0, 1) and its scale depends on
  the rate.** SD of the null statistic, by fifth of the rate within a data set:

  | Rates | Regime | 1 (smallest) | 2 | 3 | 4 | 5 (largest) |
  |---|---|---|---|---|---|---|
  | devel MLE | null | 0.26 | 0.29 | 0.31 | 0.35 | **0.69** |
  | devel MLE | near_poisson | 0.25 | 0.26 | 1.20 | 2.21 | **2.37** |
  | master | null | 0.26 | 0.29 | 0.32 | 0.34 | 0.38 |
  | cap at 50·s | null | 0.26 | 0.29 | 0.31 | 0.35 | 0.43 |
  | EB shrinkage | null | 0.27 | 0.29 | 0.31 | 0.35 | 0.38 |

  The empirical null fits one scale (about 0.3) to all genes. Genes in the top
  fifth of devel's rates have twice that scale, and they are the false
  discoveries.

* **`bool_stabilize_underdispersion` discards the absolute scale.** Whenever
  the geometric mean of the rates exceeds 1, all rates are divided by it.
  That is every eSVD-model regime under `devel`, but not `generate_null()`
  (geometric mean 0.7), nor `low_count` and `generate_null()` under `master`
  (0.5 and 0.4). *An earlier version of this document said "every regime";
  that was not checked and was wrong.* Only the *relative* rates reach the posterior. In
  `near_poisson` the geometric mean is about 1.8e4, which leaves the diverged
  genes near 550 and the others near 0.02 to 0.05: two populations, one
  collapsed onto the fit and one equal to the raw counts. Replacing the mean
  by the median does not help (Idea 7).

### 2.5 The true rates are not the target when some are very large

*Corrected after the `trend` regime was run: the first version of this section
drew its conclusion from `wide_rate` alone and stated it for all rates.*

| Rates | wide_rate: false / true discoveries | AUC, wide_rate | AUC, low_count |
|---|---|---|---|
| devel MLE | 37.8 / 20.1 | 0.870 | 0.919 |
| **oracle (true rates, fit's units)** | **5.7 / 16.6** | 0.955 | 0.947 |
| cap at 10·s | 1.0 / 18.6 | 0.976 | 0.944 |
| one common unit-free rate | 2.2 / 19.2 | 0.979 | 0.949 |

In `wide_rate` and `low_count`, a single rate shared by all genes ranks DE
genes at least as well as the truth does. In `trend`, where the true rates
stay below about 50 and follow expression, the truth ranks best (AUC 0.975
against 0.955 for the common rate) and has 1.5 false discoveries. The
failure is confined to genes whose true rate is above about 30 times the
library size. In eSVD-DE the nuisance acts as the weight
of the fit against the data in the posterior, and the test is calibrated when
that weight is comparable across genes. This is `test-calibration-and-dispersion`'s
point ("biased but stable" beats "unbiased but noisy", `svensson-2025`) in a
stronger form.

---

## 3. The dry-runs

Scripts are in `additional_context/overdispersion_brainstorm/`;
`run_all_claude.sh` reruns everything, and knits the report, in **about 4
minutes** on 8 cores (the knit needs pandoc). It
needs `version_comparison/run_all_claude.sh` to have been run first (devel
library, the six original regimes, master's output).

| Script | What it does |
|---|---|
| `00_simulate_extra_claude.R` | three added regimes (below) |
| `helpers_simulate_claude.R` | the generator with the rate tied to expression (`trend`) |
| `01_cache_fits_claude.R` | fits all 100 data sets once with devel |
| `01b_run_master_claude.R` | runs `master` itself on the four added regimes |
| `helpers_candidates_claude.R` | the likelihood in R, and every candidate as a function returning a rate vector |
| `02_run_candidates_claude.R` | overwrites the rate, reruns `compute_posterior` → `compute_test_statistic` → `compute_pvalue` |
| `03_summarize_claude.R` | `output/summary_discovery.csv`, `summary_boundary.csv`, `summary_null_scale.csv` |
| `04_downstream_claude.R` | Ideas 8 and 9 → `output/summary_downstream.csv` |
| `05_cap_sweep_claude.R` | the cap multiplier and the prior variance → `output/summary_cap_sweep.csv` |
| `06_more_cells_claude.R` | three regimes at 150 cells per individual → `output/summary_more_cells.csv` |
| `07_deseq2_variants_claude.R` | the rates scored as estimates, and the shrinkage with DESeq2's trend and sampling variance → `output/summary_accuracy.csv`, `summary_deseq2_variants.csv` |
| `08_report_tables_claude.R` | the small tables the report reads → `output/report_*.csv` |
| `overdispersion_cap_claude.Rmd` | the report on the chosen cap |

**Regimes.** The six of the version comparison, plus three added because the
original six cannot separate candidates on power (every candidate finds all 15
DE genes): `weak_de` (30 DE genes, effect ±0.3), `low_count` (gene intercepts
−2.5 to −0.5, effect ±0.6), `wide_rate` (true rates 0.5 to 200, effect ±0.3).
A fourth, `trend`, was added for the report (section 3.6). All have 20
individuals × 30 cells, 300 genes, 10 replicates.

**Sanity checks passed.** The `mle` candidate reproduces the comparison
report's devel numbers (2.7 false discoveries under the null, 63.3 and type-I
error 0.28 in `near_poisson`), and the `master` candidate reproduces
`devel_swap`.

### 3.1 False discoveries per data set at FDR < 0.05 (about 285 null genes)

| Rates | null | baseline | weak_de | low_count | wide_rate | near_poisson | generate_null |
|---|---|---|---|---|---|---|---|
| devel MLE (current) | 2.7 | 3.8 | 3.5 | 8.5 | 37.8 | 63.3 | 3.8 |
| master / legacy cap `max_i s` | 0.2 | 1.8 | 1.1 | 0.6 | 1.2 | 2.6 | 0.2 |
| cap at 10 · median s | 0.2 | 2.1 | 1.2 | 0.6 | 1.0 | 1.8 | 0.1 |
| cap at 50 · median s | 1.2 | 2.7 | 1.8 | 3.7 | 2.1 | 2.1 | 1.6 |
| cap at 3 · interior median rate | 0.3 | 1.7 | 1.5 | 2.4 | 1.0 | 2.0 | 0.0 |
| EB shrinkage (MAP) | 0.2 | 1.8 | 1.2 | 0.3 | 1.7 | 1.9 | 0.3 |
| profile lower bound, 90% | 0.3 | 1.9 | 1.2 | 1.2 | 1.3 | 2.6 | 0.1 |
| profile lower bound, 50% | 0.5 | 2.0 | 2.0 | 3.3 | 4.3 | 4.0 | 0.7 |
| common unit-free rate | 0.2 | 1.9 | 1.2 | 0.2 | 2.2 | 1.8 | 1.5 |
| common raw rate | 0.1 | 1.5 | 1.1 | 2.3 | 1.3 | 1.9 | 3.3 |
| boundary genes → median | 1.2 | 3.0 | 2.5 | 3.0 | 4.5 | 4.9 | 1.6 |
| winsorize at the 95th percentile | 0.4 | 1.8 | 1.6 | 4.8 | 37.8 | 55.3 | 1.0 |
| MLE, median in place of geometric mean | 2.9 | 4.0 | 3.5 | 8.6 | 41.2 | 31.2 | 3.8 |
| oracle | 0.2 | 2.0 | 1.2 | 0.2 | 5.7 | 2.9 | — |

### 3.2 True discoveries per data set (DE genes: 15, 30, 30, 30, 15, 10)

| Rates | baseline | weak_de | low_count | wide_rate | near_poisson | generate_null |
|---|---|---|---|---|---|---|
| devel MLE (current) | 15 | 17.5 | 14.5 | 20.1 | 12.7 | 8.2 |
| legacy cap `max_i s` | 15 | 18.2 | 12.3 | 18.9 | 15.0 | 8.8 |
| cap at 10 · median s | 15 | 17.9 | 11.9 | 18.6 | 15.0 | 8.1 |
| cap at 50 · median s | 15 | 17.6 | 13.6 | 17.6 | 15.0 | 8.1 |
| EB shrinkage (MAP) | 15 | 17.4 | 11.2 | 17.9 | 15.0 | 8.2 |
| profile lower bound, 90% | 15 | 17.7 | 12.3 | 17.6 | 15.0 | 8.4 |
| common unit-free rate | 15 | 18.1 | 11.9 | 19.2 | 15.0 | 8.0 |
| oracle | 15 | 18.4 | 12.2 | 16.6 | 15.0 | — |

The MLE's extra discoveries in `low_count` and `wide_rate` come with a worse
ranking (AUC 0.919 and 0.870 against 0.944 and 0.976 for the cap at 10·s), so
they are a looser threshold and not more power.

### 3.3 How tight must the cap be? (`β_j ≤ c · median_i s_ji`)

| c | null | baseline | low_count | wide_rate | near_poisson | generate_null | Genes capped, baseline |
|---|---|---|---|---|---|---|---|
| 1 | 0.2 | 1.9 | 0.4 | 1.7 | 1.8 | 0.2 | 100% |
| 2 | 0.2 | 1.9 | 0.2 | 1.4 | 1.8 | 0.1 | 97% |
| 5 | 0.2 | 2.0 | 0.3 | 1.1 | 1.8 | 0.0 | 41% |
| **10** | 0.2 | 2.1 | 0.6 | 1.0 | 1.8 | 0.1 | **11%** |
| 20 | 0.5 | 2.5 | 1.6 | 0.9 | 1.8 | 0.4 | 3% |
| 50 | 1.2 | 2.7 | 3.7 | 2.1 | 2.1 | 1.6 | 1% |
| 200 | 2.2 | 3.7 | 6.2 | 7.0 | 3.8 | 3.0 | 1% |
| 1000 | 2.9 | 4.0 | 7.2 | 24.2 | 8.3 | 3.5 | 1% |
| none | 2.7 | 3.8 | 8.5 | 37.8 | 63.3 | 3.8 | 0% |

(False discoveries per data set.) The response is monotone and flat up to
about c = 10, then rises. At c = 10 about nine genes in ten keep their own
estimate in ordinary data; at the legacy cap (between c = 2 and c = 3 here)
almost none do.

**A cap justified by biology is too loose.** In NB2 terms the model's size is
`θ_ji = μ_ji β_j` with `μ ≈ 1`, so the technical floor `θ ≈ 100` of
`lause-2021` and `sarkar-stephens-2021`, or nebula's `φ ≤ 1000`
(`pkg-nebula`), would mean a cap near `β = 100` to `1000`, i.e. c = 50 to 500
here. Those leave 1.2 to 2.9 false discoveries under the null. The cap that
calibrates this test is tighter than any statement about technical noise
supports, and the documentation should say that it is a calibration device.

### 3.4 Five times the cells (150 per individual, 3000 in all)

Added after the first draft, to test whether the likelihood-based bounds keep
working as the per-gene estimate sharpens. False / true discoveries per data
set, 10 replicates (SE of the false discoveries in parentheses):

| Rates | null | weak_de (30 DE) | wide_rate (30 DE) | Largest unit-free rate, wide_rate |
|---|---|---|---|---|
| devel MLE | 0.6 / — | 3.3 / 29.4 | 20.1 (1.9) / 26.7 | 1.5e7 |
| legacy cap `max_i s` | 0.4 / — | 2.9 / 29.2 | 3.2 (0.6) / 29.3 | 6 |
| cap at 10 · median s | 0.5 / — | 3.2 / 29.4 | **3.0 (0.6)** / 28.5 | 10 |
| EB shrinkage (MAP) | 0.6 / — | 3.2 / 29.4 | **6.2 (1.1)** / 28.2 | 120 |
| EB shrinkage, then cap at 10 · median s | 0.6 / — | 3.2 / 29.4 | 3.0 (0.6) / 28.5 | 10 |
| profile lower bound, 90% | 0.6 / — | 3.2 / 29.6 | 5.7 (1.0) / 28.6 | 66 |
| oracle | 0.6 / — | 3.2 / 29.6 | 10.1 (1.5) / 27.3 | 197 |

* **With true rates of 2 to 8 and 3000 cells, no gene is at the boundary and
  the MLE needs no fix** (0.6 against 0.4 to 0.6). The problem in ordinary
  regimes is a small-sample one.
* **With widely spread true rates the problem does not go away with more
  cells**, because it is the truth that is spread (oracle: 10.1). 6% of genes
  are still at the boundary.
* **The likelihood-based bounds move toward the oracle as cells are added.**
  The median sampling variance of a gene's log rate fell from 0.19 to 0.056,
  the estimated prior variance was 3.6 (the truth is 3.0), and the prior
  therefore moved the interior genes very little: the largest shrunken rate is
  120 times the library size. The cap does not depend on the number of cells.
* The false discoveries of about 3 in `weak_de` are shared by every candidate
  including the caps (pooled FDP about 0.09), so they are not a nuisance
  effect.

### 3.5 Is the rate a better estimate than `master`'s?

Sections 3.1 to 3.4 score the test. This one scores the rate itself against
the generator's truth, in the fit's units (rate divided by the gene's median
fitted library size), 600 cells.

**Genes within twofold of the true rate:**

| Rates | null | baseline | low_count | wide_rate | near_poisson |
|---|---|---|---|---|---|
| master | 0.60 | 0.59 | 0.60 ᵃ | 0.36 ᵃ | 0.00 |
| devel MLE | 0.89 | 0.89 | 0.72 | 0.60 | 0.07 |
| cap at 10 · median s | 0.96 | 0.95 | 0.84 | 0.60 | 0.00 |
| EB shrinkage | 0.99 | 0.99 | 0.89 | 0.68 | 0.05 |
| profile lower bound, 90% | 0.91 | 0.91 | 0.83 | 0.53 | 0.00 |
| common unit-free rate | 0.93 | 0.94 | 0.79 | 0.25 | 0.03 |

**Spearman correlation with the true rate, within a data set:**

| Rates | null | baseline | low_count | wide_rate | near_poisson |
|---|---|---|---|---|---|
| master | 0.05 | 0.08 | 0.02 ᵃ | 0.46 ᵃ | −0.01 |
| cap at 10 · median s | 0.77 | 0.77 | 0.64 | 0.86 | 0.03 |
| EB shrinkage | 0.77 | 0.77 | 0.65 | 0.88 | 0.08 |

ᵃ When this table was first made `master` had not been run on the added
regimes, and the legacy cap applied to devel's estimate stood in for it.
`master` itself has since been run on them (`01b_run_master_claude.R`) and
gives the same figures.

* **Both the cap at 10 and the shrinkage are better estimates than `master`'s
  wherever the rate is identifiable.** `master`'s rates are a median of 43%
  too low (0.57 of the truth) and carry no ranking; the cap and the shrinkage
  are 9% to 14% high and rank the genes at 0.77.
* **The shrinkage is the most accurate estimate; the cap is a close second**
  (99% against 95% within twofold in ordinary regimes). With the cap applied
  after the shrinkage, the 3000-cell results equal the cap's (section 3.4).
* **In `near_poisson` nothing improves on `master` in any useful sense.**
  Every candidate is off by a factor of 8 or more and none ranks the genes.
  The data carry no information about the rate there; a bounded rate is a
  choice, not an estimate.
* The common rate lands within twofold for most genes in ordinary regimes
  only because the true rates there span a factor of 4.

**DESeq2's two design choices, tried on the shrinkage** (false discoveries
per data set; prior variance in parentheses):

| Variant | null | low_count | wide_rate | near_poisson | generate_null |
|---|---|---|---|---|---|
| prototype (constant centre, observed information) | 0.2 (0.25) | 0.3 (0.25) | 1.7 (2.8) | 1.9 (0.25) | 0.3 (6.2) |
| centre is a fitted trend in the mean count | 0.2 (0.25) | 0.3 (0.25) | 1.7 (2.8) | 1.9 (0.25) | 0.3 (6.8) |
| sampling variance `trigamma((m − p)/2)` | 0.2 (0.28) | 0.5 (0.34) | 1.7 (2.9) | 1.9 (0.67) | 0.3 (6.2) |
| both | 0.2 (0.28) | 0.6 (0.36) | 1.7 (3.0) | 2.0 (0.66) | 0.3 (6.8) |

* **Neither changes the result in these regimes, which have no trend to
  find.** The generator draws the rate independently of the gene's
  expression (fitted slopes −0.10 to 0.01; −0.28 under `generate_null()`).
  **In the `trend` regime of section 3.6 the fitted trend does matter**: 97%
  of rates within twofold against 85%, and 24.3 true discoveries against
  22.7. Whether real data have a trend is part of Idea 11. The prototype's centre already scales with the gene's library
  size, which holds the gene intercept, so it is a trend with slope fixed at
  1 in the raw rate, not a flat prior.
* **`trigamma((m − p)/2)` is 0.003 at 600 cells, against a measured sampling
  variance of 0.10 to 2.1.** The formula is the variance of the log of a
  scaled chi-square with `m − p` degrees of freedom, which describes a bulk
  experiment with moderate counts. At single-cell counts most cells carry
  little information about the rate
  (`identifiability-of-dispersion-at-low-counts`), so it understates the
  noise 30- to 700-fold and the estimated prior variance comes out larger.
  The effect on false discoveries is within noise at 600 cells. It acts in
  the direction of weaker shrinkage, which section 3.4 shows is the unsafe
  direction.

### 3.6 What the report on the cap added

`overdispersion_brainstorm/overdispersion_cap_claude.Rmd` (knitted `.html`
beside it) was written after the decision. It added three things.

**`master` itself on every regime.** Its false and true discoveries on the
added regimes equal those of the legacy cap applied to devel's estimate (the
largest difference is 0.1 false discoveries, in `low_count`).

**A regime with a real trend between overdispersion and expression.**
`trend` uses DESeq2's parametric form `α(m) = α₀ + α₁/m` (`pkg-deseq2`,
`trended-dispersion`), which in the generator's units is
`β_j = 1 / (α₀·exp(intercept_j) + α₁)` times a log-normal factor, with
`α₀ = 0.1`, `α₁ = 0.02`, noise SD 0.3, intercepts −2 to 2, 30 DE genes of
effect ±0.5. True unit-free rates run from about 1 (highly expressed) to
about 50 (lowly expressed).

| Rates | False discoveries | True discoveries (of 30) | AUC | Within twofold of the truth |
|---|---|---|---|---|
| master | 2.3 | 20.7 | 0.954 | 0.33 |
| devel MLE | 28.7 | 20.4 | 0.891 | 0.69 |
| cap at 10 · median s | 1.0 | 21.7 | 0.965 | 0.75 |
| EB shrinkage, constant centre | 1.0 | 22.7 | 0.968 | 0.85 |
| EB shrinkage, fitted trend | 1.6 | 24.3 | — | 0.97 |
| common unit-free rate | 3.4 | 20.9 | 0.955 | — |
| oracle | 1.5 | 24.3 | 0.975 | 1.00 |

* **Here the cap is better than `master` on both counts**: paired over
  replicates, −1.3 false discoveries (SE 0.6) and +1.0 true discoveries (SE
  0.45). `master`'s rate is flat in expression where the truth falls
  twentyfold, and its false discoveries come from the highly expressed third
  (2.1 of 2.3).
* **Too tight a cap is worse here**: 3.5 and 2.8 false discoveries at c = 1
  and c = 2, against 1.0 at c = 10 and 2.4 at c = 50.
* **The cap's cost is in the lowly expressed third**, whose true rate is near
  20: 3.0 true discoveries of 10.7 against 5.9 with the true rates (and 2.3
  under `master`).
* **12% of genes are at the boundary under devel**, mostly lowly expressed.

**The condition for a gene to be at the boundary.** Expanding the
log-likelihood around its Poisson limit,
`ℓ(β) = ℓ_Poisson + D/(2β) + O(β⁻²)` with
`D = Σ_i [(A_i − m_i)² − A_i] / μ_i`. The sign of `D` agrees with whether
devel returned a finite rate for 11,997 of 12,000 genes (four regimes, all
replicates). This is the score statistic for overdispersion, and it replaces
the looser statement in section 2.3 that boundary genes have a Pearson
statistic below 1.

---

## Idea 1: Unit-free cap, `β_j ≤ c · median_i(s_ji)` with c = 10

**Motivation:** The smallest change that restores `master`'s calibration while
letting most genes keep their own estimate.

**What it is:** In `.estimate_nuisance_matrix()` (R, not C++), after the
per-gene estimate: `pmin(nuisance_vec, max_rate_multiplier * apply(library_mat, 2, median))`.
A new argument (default 10), the number of capped genes stored in `param`, and
a per-gene logical stored beside `nuisance_vec`. `gamma_rate` stays an honest
MLE, so T-CPP-GAM-05 and its siblings are untouched. The same line goes into
`compute_test_per_gene()`.

**Existing work:** `min_val` and `nuisance_lower_quantile` already bound the
other tail; this is their counterpart. Before the rescaling of section 2.4 it
reads as "the posterior puts at most `c/(1 + c)` of its weight on the fit";
after the rescaling it bounds how far a gene's weight can exceed the typical
gene's.

**What's missing:** A default for c chosen on something other than the
simulations that evaluate it; a check on real data (Idea 11).

**Go/No-Go task:** Done in simulation (section 3.3): **go**. False discoveries
at c = 10 match the legacy cap in all nine regimes (largest difference 0.8, in
`near_poisson` and in favor of c = 10; SE about 0.5). Remaining check: Idea 11.

**Time estimate:** 1–2 days including tests and documentation
**Risk:** Low — the behavior between c = 1 and c = 10 is flat, so the exact
default matters little.
**Difficulty:** Low
**Impact:** ★★★★★

**Recommendation:** Adopt for 1.1.0, with Idea 6.

## Idea 2: Restore the legacy cap, `β_j ≤ max_i(s_ji)`

**Motivation:** Reproduces the published method's numbers.

**What it is:** The same line as Idea 1 with `max` in place of `10 * median`,
applied in R. Not a revert of `gamma_rate.cpp`: the old bracket returned the
bound as though it were the root, which the new tests rightly reject.

**Existing work:** Verified equivalent to `master` here (rates within 0.1%
median, Welch statistics within 0.025).

**What's missing:** An honest description. With this cap the nuisance is the
largest library size for about 90% of genes, which is hard to call an
estimate of overdispersion in the documentation. `max` is also sensitive to a
single extreme cell.

**Go/No-Go task:** Done: **go** on calibration. Open: whether exact
reproduction of the 2024 results is a requirement for the CRAN release.

**Time estimate:** 1 day
**Risk:** Low
**Difficulty:** Low
**Impact:** ★★★★☆

**Recommendation:** Offer as an option of Idea 1's argument (for example
`max_rate = "legacy"`) rather than as the default.

## Idea 3: Empirical-Bayes shrinkage of `log β` (MAP under a Normal prior)

**Motivation:** The standard answer in genomics to a per-gene dispersion that
is sometimes unidentified (`empirical-bayes-shrinkage-dispersion`,
`log-normal-dispersion-prior`). A boundary gene has a likelihood that is flat
to the right, so the prior decides and the estimate is finite with no cap.

**What it is:** Maximize `ℓ_j(ρ) − (ρ − m_j)² / (2τ²)` in `ρ = log β`, with
`m_j = log median_i(s_ji) + centre`. The centre is the median unit-free log
rate of the interior genes; `τ²` is the spread of the lower half of those
genes minus the median sampling variance, floored at DESeq2's 0.25. One more
term in the C++ derivative class, or a one-dimensional `optimize()` in R.

**Existing work:** `helpers_candidates_claude.R::.map_rates()` is a working
prototype. `pkg-glmpca` and `pkg-deseq2` are the references for the design.

**How the prototype behaved at 600 cells.** In the ordinary regimes the
median sampling variance of a gene's log rate is 0.096 against a prior
variance of 0.25, so a typical gene is pulled about a quarter of the way to
the centre (interquartile range of shrunken / MLE: 0.86 to 1.07). The most
overdispersed tenth of genes moves up by 10%. Boundary genes land at about 10
times the library size (2.3 times the centre), and no shrunken rate exceeds
11 to 15 times. **In effect it is a soft cap at the same place as Idea 1's
hard cap**, which is why the two are indistinguishable there.

**What's missing:** **Its strength depends on the number of cells** (section
3.4): at 3000 cells and widely spread true rates it leaves 6.2 false
discoveries against the cap's 3.0. This is the prior doing what it is
designed to do, recovering the true rates, when the true rates are not what
calibrates the test (section 2.5). The hyperparameters are estimated from
interior genes only. In nearly Poisson data these are the genes that happened to look
overdispersed, so the centre is low (about 26·s against a truth near 240·s).
That errs toward caution, but it is a selection effect that would need
stating. The floor bound in most replicates of seven of the nine regimes, so
in practice `τ² = 0.25` was used, not estimated.

**Go/No-Go task:** Done: **go at 600 cells, no-go as the only bound.** At 600
cells, false discoveries 0.2 to 1.9 across regimes, the same as the caps.
Sweeping `τ²` from 0.05 to 4 moves false discoveries under the null from 0.2
to 0.5 and in `low_count` from 0.1 to 1.9, so the result is insensitive below
`τ² = 1`. At 3000 cells in `wide_rate`: 6.2 against 3.0.

**Time estimate:** 1–2 weeks (C++ or R implementation, tests against the
prototype, documentation of the prior)
**Risk:** Medium — real cohorts have more cells than 3000, so the regime
where the prior stops protecting is the realistic one. Also hyperparameter
estimation on real data with thousands of genes and a real mean–dispersion
trend.
**Difficulty:** Medium
**Impact:** ★★★☆☆ (a better *estimate* than a capped MLE, which is worth
having for reporting the overdispersion itself; not a substitute for a cap)

**Recommendation:** Not as the fix for 1.1.0 on its own. Shrinkage followed
by Idea 1's cap gives the most accurate estimate (section 3.5) with the cap's
calibration (section 3.4); the gain over the cap alone is 99% against 95% of
genes within twofold at 600 cells, and nothing at 3000. The alternative is to
fix `τ²` at a small value by design, which makes it a common rate with
per-gene adjustment and should be described that way.

## Idea 4: Conservative rate — the lower end of the profile-likelihood interval

**Motivation:** `confidence-intervals-for-dispersion` calls the unexposed
profile-likelihood interval "the most concrete unexploited result in the
wiki". The lower end for `β` (the upper end for the overdispersion) exists even
when the maximizer does not, and using it propagates the uncertainty of the
estimate in the cautious direction, which addresses the plug-in problem of
`law-2014`.

**What it is:** `β_j^L = min{β : 2[ℓ_j(β̂) − ℓ_j(β)] ≤ χ²_{1, level}}`, by
`uniroot()` or one more bracketing pass in C++.

**Existing work:** Prototype in `.gene_summaries()`.

**What's missing:** A principled level. No hyperparameters shared between
genes, which is its advantage over Idea 3.

**Go/No-Go task:** Done: **go at level 0.9 and 600 cells, no-go at 0.5, and
no-go as the only bound.** At 0.9 it matches the caps at 600 cells (0.3 under
the null, 2.6 in `near_poisson`). At 0.5 it leaves 3.3 to 4.3 false
discoveries in `low_count`, `wide_rate` and `near_poisson`. At 3000 cells in
`wide_rate` the 0.9 bound leaves 5.7 against the cap's 3.0 (section 3.4): the
weakening predicted under Risk below, now observed.

**Time estimate:** 1 week
**Risk:** Medium — with many more cells per gene the interval narrows and the
bound approaches the MLE, so the protection against the mechanism of section
2.5 (rates that are well estimated but far apart) weakens as data grow.
**Difficulty:** Medium
**Impact:** ★★★☆☆

## Idea 5: One unit-free rate for all genes (or a trend in expression)

**Motivation:** Section 2.5: a common rate ranks at least as well as the
truth. `lause-2021` recommends a single shared `θ` for normalization, and
scVI's default shares one dispersion per gene across all cells
(`pkg-scvi-tools`).

**What it is:** `β_j = c · median_i(s_ji)` with `c` the median unit-free MLE
over interior genes. A trend in mean expression (`trended-dispersion`) is the
elaboration.

**Existing work:** Candidates `common_unit` and `common_raw`.

**What's missing:** It fails on the one generator not written from the eSVD
model: under `generate_null()`, 1.5 (unit-free) and 3.3 (raw) false
discoveries against 0.1 to 0.3 for the capped per-gene estimates. There the
genes differ in within-individual variance by design, and the per-gene
estimate is needed.

**Go/No-Go task:** Done: **no-go as the default.** Useful as a diagnostic: if
a data set's results change much between Idea 1 and a common rate, the rates
are doing real work.

**Time estimate:** 1 day
**Risk:** Medium
**Difficulty:** Low
**Impact:** ★★☆☆☆

## Idea 6: Record and report the genes at the boundary or at the cap

**Motivation:** Q-LFC-1. Whatever bound is chosen, a user should be able to
see which genes it acted on, and a gene whose rate is a bound has a standard
error that reflects the bound.

**What it is:** `estimate_nuisance()` stores a per-gene status (`"interior"`,
`"capped"`, `"boundary"`) and the count in `param`; `report_results()` gets
the column. The boundary test is one comparison of two log-likelihoods.

**Existing work:** `param$nuisance_num_failed` is the precedent.

**What's missing:** —

**Go/No-Go task:** Done, as a *sole* fix: **no-go.** Replacing only the
boundary genes by the median leaves 1.2 false discoveries under the null and
4.9 in `near_poisson`, because the tail is continuous: 0.5% to 13% of genes
are interior with a unit-free rate above 50 (`near_poisson` 13%, `wide_rate`
7%). As a complement to Idea 1 or 3 it costs almost nothing.

**Time estimate:** 1 day
**Risk:** Low
**Difficulty:** Low
**Impact:** ★★★☆☆

**Recommendation:** Do it together with whichever of Ideas 1 to 4 is chosen.

## Idea 7: Rescale by the median in `bool_stabilize_underdispersion`

**Motivation:** Diverged genes inflate the geometric mean that every rate is
divided by (1.8e4 in `near_poisson`).

**What it is:** `median(log10(nuisance_vec))` in place of the mean.

**Go/No-Go task:** Done: **no-go alone** (2.9 under the null, 31.2 in
`near_poisson`; the diverged genes are still far from the rest). **Unnecessary
after a cap** (cap at 50·s with either rescaling: false discoveries within 0.2
of each other in every regime).

**Time estimate:** 1 hour
**Risk:** Low
**Difficulty:** Low
**Impact:** ★☆☆☆☆

A separate question this raised, not answered here: since the rescaling
discards the absolute scale of the rates in most regimes tried, the absolute
level of the posterior weight is set by convention (geometric mean 1) and not
by the data. Whether that is intended deserves a sentence in the
documentation.

## Idea 8: Welch variance between individuals only

**Motivation:** Remove the dependence of the statistic's scale on the rate at
its source (section 2.4; open question 4 of `CLAUDE_kevin.md`).

**What it is:** `var_i(m_i)` with a Bessel correction in the denominator, in
place of the mixture variance.

**Go/No-Go task:** Done: **no-go as tried.** The null SD of the Gaussianized
statistic is 1.45 to 1.5 (3.2 at low counts), so the theoretical null gives
type-I error 0.19 to 0.64: the individuals' means share the fitted
case-control coefficient and the fit's noise, which the between-individual
variance does not see. With the empirical null the test is calibrated but
finds 3 of 30 DE genes in `weak_de` (against 18) and none in `low_count`.

**Time estimate:** —
**Risk:** High
**Difficulty:** Medium
**Impact:** ★☆☆☆☆ in this form

## Idea 9: Rescale the statistic within strata of the rate

**Motivation:** Keep the MLE and let the empirical null see statistics of one
scale.

**What it is:** Divide the Gaussianized statistic by a robust scale computed
within tenths of the rate, then `multtest()`.

**Go/No-Go task:** Done: **no-go.** With 30 genes per stratum the scale is too
noisy: 3.7 false discoveries under the null with the MLE (against 2.7 without
the rescaling) and 0.5 to 1.9 with bounded rates (against 0.2 to 1.2). It may
behave differently with thousands of genes; not worth pursuing before Idea 1.

**Time estimate:** —
**Risk:** High
**Difficulty:** Medium
**Impact:** ★☆☆☆☆

## Idea 10: A variance for the statistic that accounts for the fit

**Motivation:** The root cause is that the Welch denominator is a posterior
variance conditional on the fit, while the numerator's sampling variability
includes the fit's. Every idea above works around this.

**What it is:** Candidates, none tried: a sandwich or delta-method variance
for the fitted case-control coefficient (`test-calibration-and-dispersion`
lists the sandwich route via `sun-2026`); count splitting so that the rate and
the posterior are computed on counts the fit did not see; a permutation null
over individuals' labels, which needs a refit per permutation.

**Existing work:** T-LFC-17's finding that the SE is "about right" against
refitted cohorts (calibration ratio 0.77) on genes with an ordinary rate
suggests the present variance is a workable stand-in when rates are bounded.

**What's missing:** Everything.

**Go/No-Go task:** On the cached `null` fits, compare the SD across the 10
replicates of each gene's numerator (`case_mean − control_mean`) with the
package's denominator, by fifth of the rate. If their ratio is flat in the
rate once rates are capped, the present variance is adequate and this idea is
not needed. Half a day.

**Time estimate:** 2–4 months
**Risk:** High
**Difficulty:** High
**Impact:** ★★★☆☆ for the package now; ★★★★★ as a methods question

## Idea 11: Measure the boundary and the cap on one real data set

**Motivation:** Every number here is from simulation, mostly from the eSVD
model itself, at mean counts of 0.1 to 4.5 and 300 genes. The fraction of
genes at the legacy cap decides how different Idea 1 is from Idea 2 in
practice.

**What it is:** Run devel on one of the paper's data sets (`PAPER_DATA`, no
path recorded in `CLAUDE_kevin.md`; or the ASD vignette's input) and record
the fraction of genes at the boundary, above `max_i s`, and above
`10 · median_i s`; then compare the DE lists under Idea 1, Idea 2 and
`master`.

**Go/No-Go task:** If under 5% of the FDR < 0.05 genes differ between Idea 1
and Idea 2, choose Idea 1 without concern for reproducing the paper. If more
differ, the choice of default is a decision about continuity with the 2024
results. One day once the data are at hand.

**Time estimate:** 1–2 days
**Risk:** Low
**Difficulty:** Low
**Impact:** ★★★★☆

---

## Summary table

| # | Idea | Impact | Time | Risk | Difficulty | Priority | Go/No-Go |
|---|---|:---:|:---:|:---:|:---:|:---:|---|
| 1 | Unit-free cap, c = 10 | ★★★★★ | 1–2 days | Low | Low | **Critical** | **Go**: matches master in 9 of 9 regimes |
| 6 | Record capped / boundary genes | ★★★☆☆ | 1 day | Low | Low | **High** | Go as complement; no-go alone |
| 11 | One real data set | ★★★★☆ | 1–2 days | Low | Low | **High** | Pending: needs `PAPER_DATA` |
| 2 | Legacy cap `max_i s` | ★★★★☆ | 1 day | Low | Low | Medium | Go; 90% of genes at the cap |
| 3 | EB shrinkage (MAP) | ★★★☆☆ | 1–2 wks | Medium | Medium | Low | Same as 1 at 600 cells; 2× the cap's false discoveries at 3000 |
| 4 | Profile lower bound | ★★★☆☆ | 1 wk | Medium | Medium | Low | As 3; no-go at the 50% level |
| 5 | Common rate | ★★☆☆☆ | 1 day | Medium | Low | Low | No-go: fails on `generate_null()` |
| 7 | Median rescaling | ★☆☆☆☆ | 1 hr | Low | Low | Low | No-go alone; unnecessary after a cap |
| 8 | Between-individual variance | ★☆☆☆☆ | — | High | Medium | Drop | No-go: null SD 1.5, power lost |
| 9 | Rate-stratified rescaling | ★☆☆☆☆ | — | High | Medium | Drop | No-go: worse than no rescaling |
| 10 | Variance that accounts for the fit | ★★★☆☆ | 2–4 mo | High | High | Later | Half-day check described |

## Prioritization Recommendation

**Highest priority (do first):**
- Idea 1 with Idea 6: the cheapest change that restores calibration, it keeps
  `gamma_rate` and its tests as they are, and the result is insensitive to the
  default between c = 1 and c = 10.
- Idea 11: one day, and it decides whether Idea 1's default or Idea 2's should
  ship.

**Second priority:**
- Idea 2 as a named option for reproducing the published results.
- Idea 3 only in combination with a cap. It is the better estimator of the
  rate, and the worse guarantee for the test.

**Lowest priority / most speculative:**
- Idea 10, a methods project and not a CRAN item.
- Ideas 5, 7, 8, 9: tried, and no-go for the reasons given.

---

## What this does not settle

* **Simulated data only**, eight of nine regimes from the eSVD model itself.
  `generate_null()` is the one outside check, and it is where the common rate
  failed.
* **The default c = 10 was read off the same simulations that evaluate it.**
  The flat response from c = 1 to c = 10 makes this less worrying, but it is
  not a held-out result.
* **Ten replicates.** The SE of a false-discovery count in the tables is about
  0.5 (0.1 to 1.2 by cell), so the bounded candidates cannot be ranked against
  each other. The contrast with the MLE is far outside that noise.
* **Power was only weakly tested.** Six regimes are saturated (all DE genes
  found). In `low_count` the bounded candidates find 11 to 13 of 30 against
  the MLE's 14.5, with a better AUC; whether some of that is real power lost
  to the cap is not resolved.
* **`compute_test_per_gene()`** (the `bool_diet = TRUE` path) was not run. It
  calls the same estimator, so the same bound applies, but it was not checked.
* **300 genes.** `locfdr` fell back to the truncated MLE in a minority of
  replicates, as in the comparison report. With thousands of genes the tail
  behavior may differ in size.
* **`bool_use_log = TRUE`** already bounds the rate at `exp(10) ≈ 22026`
  through `log_gamma_rate`'s default bracket. That is far above every cap
  that worked here, so it is not a fix.

## Corrections this implies elsewhere (not made)

* `CLAUDE_kevin.md` and section 4.1 of the comparison report say devel
  "overestimates the rates two- to threefold". In the fit's units the excess is
  4% to 16% (section 2.2). The report's `nuisance_true` column would need
  multiplying by `median_i s_ji` to be comparable.
* The same report's statement that master's rates "carry no information about
  the truth" is explained by section 2.1 and could say so.
