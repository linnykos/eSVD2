# Autoresearch alignment memo: reporting a log fold change and standard error from an eSVD-DE fit

Kevin Z. Lin (drafted by Claude), 2026-11-11. Round 1. Branch `2026-11-07_std`.

**Status: a proposal awaiting Kevin's review.** Nothing is built. Section 6
lists ten questions, each with a recommendation; the bank, the scorer and the
queue below assume the recommended answers. When the memo is agreed it is
copied to `autoresearch_lfcse/notes/MEMO_claude.md`, the place the build
reads it from.

Inputs read for this round: `2026-11-07_std-exploration/` (README, the
summary, `output/02_*.csv`, `output/03_*.csv`, `helpers_claude.R`), and the
290 cached fits of that exploration, rescored by
`2026-11-07_std-exploration/04_score-demo_claude.R` (one minute, no fitting).

## 0. What the exploration settled, and what it left open

The exploration compared nine fixed candidates on seven simulated scenarios
and four cohort-pair scenarios. Its conclusions, which this memo takes as the
starting point:

| Finding | Evidence |
|---|---|
| The stored SE of 1.2.0 is 3.4 to 5 times the sampling SD of the logFC | replicate cohorts, median SE / oracle SD |
| The between-donor SE of the posterior donor means is 0.65 to 0.73 of the sampling SD, because the posterior shrinks the donors of an arm toward each other by the weight `w` | same; slope of the SE ratio on `w` |
| The SE of the unshrunk donor means `A / library` (`raw`) is calibrated at the median (0.98 to 1.01) | same |
| The logFC itself is unattenuated (slope 1.00 on the truth) but carries a shift common to all genes when the DE is one-directional, which the empirical null's mean removes | `asymmetric` scenario |
| Rescaling by the empirical null's SD overstates Figure 4's `r` when the true effects are weak (0.96 to 0.97 against 0.83) | `weak_effects` pairs |

What it did not do is search. Each candidate was written down in advance,
and five places where the best of them still loses were measured and left:

1. With 10 donors, `raw_null` scores 0.09 nats per gene worse than DESeq2
   (Section 1 defines the score), the largest gap of any scenario, and the
   stored SE rescaled by the null (`null_rescaled`) beats `raw_null` there
   by 0.08. The per-gene between-donor variance from 5 donors per arm is
   noisy, and nothing in `raw` stabilizes it.
2. With 200 cells per donor, `raw` covers the truth 85% of the time and
   needs the null's rescaling to reach 95%.
3. The shift that corrects the logFC is `null_mean * SE`, proportional to
   each gene's SE, while the depth bias it corrects is one constant on the
   log scale.
4. No scenario had covariates beyond depth, unequal cells per donor, a
   donor effect other than SD 0.1, or true effects in most genes.
5. The two recommendations disagree with each other by use: the rescaled
   SE is best for coverage and worst for `r`.

The loop is for these.

## 1. The goal

**The object searched for** is a *reporting rule*: a function that takes one
fitted eSVD-DE object, with the counts and the design it was fitted to, and
returns for every gene a log2 fold change and its standard error, such that
the pair can be used as a Wald pair. A methods section would state it as an
estimator of the case-control log2 ratio of mean expression and an estimator
of that estimator's sampling SD.

**"Wald-usable" means three things**, in the order the downstream uses need
them: the logFC is unbiased for the true log2 fold change at every effect
size (no shrinkage toward zero); the SE is the SD of the logFC over cohort
replications, gene by gene; and `logFC / SE` is standard normal for a gene
with no effect.

**The score** is the Gaussian negative log-likelihood of the truth, in nats
per gene, lower is better. For one fitted cohort with genes
$g = 1, \dots, p$, true log2 fold change $\beta_g$, and reported
$(\widehat\beta_g, s_g)$,

$$\text{NLL} = \frac{1}{p}\sum_{g=1}^{p}\left[\tfrac12 \log\!\left(2\pi s_g^2\right) + \frac{(\widehat\beta_g - \beta_g)^2}{2 s_g^2}\right].$$

One example of the bank is a pair of cohorts (Section 3), and the example's
score is the mean of the two cohorts' NLL.

**Why this number.** In expectation over cohorts the NLL is smallest, gene
by gene, when $s_g^2$ equals the mean squared error of $\widehat\beta_g$,
where its value is $\tfrac12\log(2\pi e\,\text{MSE}_g)$. It therefore pays
for an SE that is too large (the first term), for one that is too small (the
second), for a bias in the logFC (which enters the MSE), and for a noisier
logFC, all in one currency: an SE off by a factor $c$ costs
$\log c + (c^{-2} - 1)/2$ nats, which is 0.16 at $c = 0.7$ and 0.87 at
$c = 3.8$. It needs no replicate cohorts, because the truth is known.

**Why not its neighbours.**

| Neighbour | Why not |
|---|---|
| Coverage of ±2 SE | The stored SE has coverage 1.00; a wide SE wins |
| Median SE / oracle SD | Blind to the logFC's bias, and to an SE that is right at the median and wrong gene by gene; needs 30 refits per example |
| Error in Figure 4's `k` and `r` | Two moments per pair; a rule can match the moments with SEs that are wrong for every gene |
| False discoveries of the Wald test | Sees only the null genes' tail |
| CRPS | Proper too, but nearly flat in the SE: a factor-0.7 SE costs about 5% |
| The NLL with the squared z-score capped | A cap at 5 SE changed the score of the three calibrated candidates by at most 0.0003, and it is what makes thresholding profitable (below) |

**The demo.** `04_score-demo_claude.R` scores every candidate on the 370
cached cohorts (210 of script 02, 160 of script 03).

| Candidate | Score (nats / gene) | Difference from `raw_null` | Paired SE |
|---|---|---|---|
| DESeq2 pseudobulk (reference) | -0.864 | -0.035 | 0.002 |
| A: null-rescaled stored SE | -0.855 | -0.025 | 0.002 |
| unshrunk + empirical null (`raw_null`) | -0.829 | 0 | |
| between / weight (`deshrunk`) | -0.613 | +0.216 | 0.021 |
| unshrunk donor means (`raw`) | -0.581 | +0.248 | 0.023 |
| between + within / cells | -0.456 | +0.373 | 0.030 |
| stored (1.2.0) | -0.013 | +0.816 | 0.005 |
| D: donor jackknife | -0.011 | +0.819 | 0.047 |
| B: between-donor | +0.067 | +0.897 | 0.050 |

The order is the exploration's, obtained from one number. The paired SD of
a difference between two calibrated candidates is 0.035 to 0.044 per cohort,
so a bank of 144 pairs resolves about 0.003 nats.

Two things in the table bear on the design. First, `null_rescaled` beats
`raw_null` by 0.025 (14 paired SE), with most of the difference in the
10-donor cohorts, and it does so in the `weak_effects` pairs as well (-0.884
against -0.855), where the exploration found its `r` to be 0.12 too high. The
NLL does not see that failure of the moments identity, which is question 1.
Second, the demo's tenth candidate, a thresholding rule that reports
$0 \pm s_g/10$ for every gene with $|z| < 2$, scores -2.70 to -3.48 in the
seven scenarios whose non-DE genes are exactly null and +0.69 to +3.62 in
the four whose effects are Gaussian. A proper score rewards shrinkage
wherever the truth is sparse; Sections 2 to 4 close this with a rule, a
gate, and a bank in which most cells have no exactly-null gene.

**Headline deliverable of a run**: the kept `report_rule_claude.R`, the
self-contained rule file, and `RUN_<tag>_claude.md`.

## 2. The file on trial

**Path**: `autoresearch_lfcse/report_rule_claude.R`. It defines

```r
report_lfc_se <- function(view)
```

**What it receives.** `view` is a list holding what a user of eSVD2 has
after `eSVD()` with `bool_diet = FALSE`, and no element derived from the
truth:

| Element | Content |
|---|---|
| `dat` | counts, cells x genes |
| `covariates` | the design matrix `format_covariates()` built (`Intercept`, `Log_UMI`, `CC_1`, others) |
| `individual_vec`, `case_control_vec` | donor label and arm of each cell |
| `x_mat`, `y_mat`, `z_mat` | the final fit |
| `library_mat` | fitted library size of each cell and gene (`bool_covariates_as_library = TRUE`) |
| `nuisance_vec`, `nuisance_status_vec` | the rate the posterior used, and `estimated` / `capped` / `boundary` |
| `posterior_mean_mat`, `posterior_var_mat` | per cell and gene |
| `stored_df` | per gene: 1.2.0's `log2fc`, `log2fc_se`, Welch statistic, df, Gaussianized z, p-value; with the fitted `null_mean` and `null_sd` as attributes |

**What it returns**: a data frame with one row per gene in the order of
`colnames(view$dat)` and columns `gene`, `log2fc`, `se`. The scorer refuses
a wrong row count or order, a non-finite value, and `se <= 0`.

**What it may call**: base R, `stats`, `Matrix`, `locfdr`, and eSVD2's
exported functions (in practice `multtest()`). It runs sequentially, writes
nothing, and has 30 seconds per cohort.

**Generality requirement.** Every cohort is treated by the same formula.
A numeric constant carries a comment giving the reason for its value, and
the reason is not "it scored better". The rule must be one a methods
section can state in display mathematics.

**Code it leaves alone**: `R/`, `src/`, the harness, the bank. A kept rule
is a proposal for `compute_log_fold_change()`; moving it into the package is
Kevin's decision after the run.

## 3. The bank

**One example** is a pair of cohorts (A, B), each fitted end to end by
eSVD2 1.2.0 in the order `eSVD()` runs its stages (cap 10, defaults). The
two populations share every gene's intercept, loadings and rate; their
case-control effects are correlated at $r = 0.8$ with a size ratio
$k = S_A / S_B = 2$, as in script 03, so that every example also yields the
Figure 4 terms. Each cohort is 600 genes, fitted at the true rank 3.

**Generator**: `helpers_claude.R` of the exploration, version
`lfcse-bank-1`, extended by the four things marked *new* below. The counts
follow the eSVD model, with donor random effects on the log scale.

**Factors.** The factorial is 3 x 2 x 2 x 3 x 2 x 2 = 144 cells.

| Factor | Levels | Reason for the levels |
|---|---|---|
| `donors` (total, half per arm) | 10, 20, 40 | Prater has 22 and the ASD cohort 31; 10 is the level at which the per-gene variance is least stable (Section 0, item 1) |
| `cells` per donor | 50 equal; unequal, log-normal with median 150 and range 30 to 600 (*new*) | the ASD cohort averages 430 cells per donor and real donors are far from equal; 200 equal cells already lowered `raw`'s coverage to 0.85 |
| `rates` (range of the true Gamma rate) | 2 to 8; 0.5 to 200 | about 11% and about half of the genes capped; the real cohorts have 73% |
| `effects` | `sparse`: 10% of genes N(0, 0.6²), the others 0; `sparse_up`: the same sizes, all positive; `dense` (*new*): every gene nonzero, 85% N(0, 0.1²) and 15% N(0, 0.3²) (natural log, cohort A) | `sparse` is script 03; `sparse_up` is the one-directional case in which the depth adjustment shifts every gene; `dense` aims at the real cohorts' reliability of 0.4 to 0.65 and has no exactly-null gene |
| `covariates` (*new*) | none; donor-level sex (binary) and age (continuous), each with gene-specific effects of SD 0.2, sex 60/40 between the arms | the real fits adjust for sex, age and batch, which enter the fitted library that `raw` divides by; the exploration had none |
| `donor_sd` (*new*) | 0.1; 0.3 | 0.1 is the exploration's only value; with it, cell-level noise dominates the between-donor variance at 50 cells. The upper level is a placeholder until question 6 anchors it on the ASD fit |

**Replicates and splits.** One example per cell per split: 144 dev and 144
test, 576 fits in all. The seed of an example is derived from its cell index
and its split.

**Truth stored per example**, apart from the public view: `truth_log2fc_vec`
of each cohort (the generator's case-control coefficient divided by
$\ln 2$, which is the log2 ratio of mean expression at equal depth and equal
covariates), the gene class (exact null, weak with $|\beta| < 0.5$, strong),
the true rates, the population `k` and `r`, and DESeq2's pseudobulk triple
for the reference line. The truth is fixed for the whole loop and is not on
trial.

**Checks at the build**, reported in `BANK.md` as realized against target:
the fraction capped per `rates` level, the realized reliability
$\rho_A = \text{mean}(\beta^2) / \text{mean}(\widehat\beta^2)$ per `effects`
level under the baseline (target 0.4 to 0.65 for `dense`), the arm imbalance
of sex, and the cells-per-donor range. The reliability needs the fits, so it
is a build-time check, not a generation-only one.

**Cost.** A 1000-cell x 300-gene fit takes 3.6 seconds and an 8000-cell one
about 30; doubling the genes roughly doubles both. The estimate is 30 to 45
minutes on 8 cores for the bank, measured at the build. Each public view
holds four cells x genes matrices, so the bank is several gigabytes; if it
passes 5 GB the posterior matrices are dropped from the files and recomputed
by the scorer from the stored fit.

**Home**: outside the repository, at `ESVD2_BANK_DIR`, default
`tools::R_user_dir("eSVD2", "cache")`, folder `lfcse-bank-1/`. The
repository sits in Dropbox, which rewrites files underneath a run. The name
`LFCSE_BANK` is added to the master `CLAUDE.md` at the build, and the path
to `CLAUDE_kevin.md`.

## 4. Directions of variation

**D1. Which donor summary is averaged.** The baseline averages
`A / library` over a donor's cells. Alternatives: the ratio of sums
$\sum_c y_c / \sum_c \ell_c$ per donor, which weights cells by depth; donors
weighted by their cell count or by their precision; the posterior means with
per-donor weights undone. *Should help* in the unequal-`cells` and
large-`donor_sd` cells, where equal weighting of a 30-cell and a 600-cell
donor is inefficient. *A discard if* the logFC's RMSE does not fall, or the
null genes' median moves from zero.

**D2. The variance form and the small-cohort correction.** Delta method on
the ratio against the variance of log donor means; a Welch-df t-quantile
inflation of the SE; a within-donor Poisson term added for donors with few
cells. *Should help* at 10 donors, where a Wald interval from 5 donors per
arm is too narrow. *A discard if* it widens the SE at 40 donors.

**D3. Centring the logFC.** A constant shift in log2 units (the mode or
median of the logFC over genes, or the empirical null's mean converted at
the median SE) against the baseline's `null_mean * SE`. *Should help* in
`sparse_up`, because the depth bias is one constant on the log scale.
*A discard if* it loses in the `dense` cells, where every gene has an
effect and the centre of the logFC distribution is not a null gene's value.

**D4. Borrowing across genes for the SE.** The global null SD (baseline),
none, a null SD fitted on a narrower central window, a null SD by stratum
(nuisance status, expression), and moderation of the per-gene between-donor
variance toward a trend in expression with a prior weight estimated from
the genes. *Should help* at 10 donors (Section 0, item 1) and in `dense`,
where `locfdr`'s null SD is inflated by weak true effects. *A discard if*
the gain at 10 donors is paid for at 40, or `r_err` rises.

**D5. Using more of the fit.** The stored within term as a stabilizer of
the between term (the reason `null_rescaled` wins at 10 donors); the
shrinkage weight per donor instead of per gene; the fit's own case-control
coefficient `z_mat[, "CC_1"]` as the logFC. *Should help* wherever the fit
is more precise than the raw ratio, which is the low-count and capped
genes. *A discard if* the SE's calibration then differs between the capped
and the estimated genes, the failure of the stored SE.

**Baseline**: `raw_null` (equation 8 of the summary), the exploration's
recommendation. `stored`, `null_rescaled`, `deshrunk` and `raw` are scored
once at the build as reference rows of `BANK.md` and kept in `candidates/`.

**Pre-registered queue**, in order:

1. [D4] Keep the null-mean shift, drop the null-SD rescaling.
2. [D3] Shift by a constant: the null mean times the median SE.
3. [D2] Inflate the SE by $t_{0.975,\,\text{df}} / 1.96$ with the Welch df of the unshrunk donor means.
4. [D1] Donor means as ratios of sums.
5. [D4] Moderate the between-donor variance toward its trend in mean expression, prior df estimated across genes.
6. [D5] SE² as the unshrunk between term plus the stored within term divided by the donor's cell count, then the null fit.
7. [D2] The SE from the variance of log donor means, with a half-count offset for zero donors.
8. [D4] The null SD fitted separately for capped and estimated genes.
9. [D1] Donors weighted by cell count in the arm means and in the variance.
10. [D5] The fit's case-control coefficient as the logFC, with the baseline SE.

After the queue, ideas come from the best commit's strata lines.

**Out of bounds this round.**

| Not tried | Reason |
|---|---|
| Any shrinkage of the logFC toward zero or toward other genes' effects (thresholding, empirical-Bayes posterior means, `apeglm`-style priors). A shift common to all genes is allowed | Figure 4 and any meta-analysis need an unbiased estimate at every effect size; a proper score rewards shrinkage, so this is a rule, backed by the gate of Section 5 |
| A logFC that ignores the eSVD fit, or a call to DESeq2, edgeR, limma, NEBULA or dreamlet | the object is eSVD-DE's reported effect; DESeq2 stays the reference line |
| Refitting: bootstrap or split-half refits of eSVD | 4 to 30 seconds each, times hundreds per cohort; a later round on a smaller bank could lift this |
| Returning or changing the p-value | see question 3 |
| Any change to `R/`, `src/`, the cap, or the posterior | the fit is not on trial |
| A branch on a factor's level inferred from the data in order to treat cells of the bank differently | generality |

## 5. The protocol

**Keep rule**: keep when the dev score falls by at least one paired SE;
keep a tie (gain in [0, 1 SE)) when the rule got shorter or twice as fast;
discard otherwise. The dependence unit is the example (the pair).

**Validity gate** (the truth is not on trial; the gate guards the
no-shrinkage rule). Over the weak-effect genes of the split
($0 < |\beta_g| < 0.5$, over 50,000 genes, most of them in the `dense`
cells), the
slope through the origin of the reported logFC on the truth must lie in
[0.95, 1.05]. The nine candidates of the demo give 1.01 to 1.02 and the
thresholding rule 0.85. A commit outside the interval is logged
`gate_failed` and never kept.

**Summary block.** `score`, `score_se`; the diagnostics `rmse`,
`coverage_95`, `z_sd_null` (SD of `logFC / SE` over exactly-null genes),
`bias_null` (their median logFC), `slope_weak`, `wald_false` (BH 0.05
discoveries per 100 exactly-null genes), `k_err` and `r_err` (mean absolute
error of the Figure 4 terms over pairs, computed as Writeup50 computes
them), `rho_err`, `score_vs_deseq2`; one `score_<factor>_<level>` line per
level; and `score_status_capped` / `score_status_estimated`. A diagnostic
that moved the wrong way goes in the row's description; the decision stays
on the score.

**Machine**: Kevin's laptop (12 cores, 32 GB). **Cores**: 8, leaving 4.
**Hours**: 4. **Experiment cap**: 40. **Stall rule**: twelve consecutive
non-keeps after fifteen experiments. **Scorer budget**: 600 seconds; the
expected scoring is 1 to 3 minutes, since nothing is refitted.

**Package version**: the scorer loads eSVD2 1.2.0 from
`additional_context/version_comparison/lib/devel`, as the exploration does;
`R/` has not changed since that library was built.

**Loop folder**: `autoresearch_lfcse/` at the repository root, tracked in
git and added to `.Rbuildignore`; `autoresearch_lfcse/output/` and
`autoresearch_lfcse/cache*/` added to `.gitignore`. **Branches**:
`autoresearch_lfcse/<mmmdd>`, cut from `2026-11-07_std` once that branch's
exploration is committed.

## 6. What I am uncertain about, stated plainly

All ten are open. Each states the default the sections above assume.

1. **Should the Figure 4 terms enter the score?** The NLL prefers
   `null_rescaled` to `raw_null` even in the `weak_effects` pairs, where
   its `r` is 0.12 too high. Options: (a) NLL alone, `k_err` and `r_err`
   as diagnostics; (b) NLL plus a penalty on the pair's moment error
   $|\text{mean}(\widehat\beta^2) - \text{mean}(s^2) - \text{mean}(\beta^2)|$
   in some exchange rate; (c) the moment error as a second gate.
   *Recommendation: (a) for this round*, with the `dense` cells in the
   bank, because the `r` overshoot is a symptom of an SE inflated among
   weak effects, which the NLL does charge for when most genes have one;
   if the first run keeps a commit whose `r_err` rose, move to (c).

2. **Is "no shrinkage of the logFC" the right rule, and is [0.95, 1.05]
   the right gate?** The alternative is to let the loop find a shrunken
   estimate and report both. *Recommendation: the rule as written.*

3. **Is the p-value on trial?** The stored empirical-null p-value makes 2
   to 3 false discoveries per 270 null genes and is not broken.
   *Recommendation: no.* The rule returns the logFC and SE only; the
   scorer reports the Wald test's false discoveries as a diagnostic, and
   whether `report_results()` should show a Wald p-value beside the stored
   one is decided after the run.

4. **Must the logFC be built from the fit?** A pseudobulk ratio with its
   own normalization would probably score near DESeq2's -0.864.
   *Recommendation: yes, it must use the fitted library or the posterior*;
   otherwise the loop converges on reimplementing pseudobulk.

5. **Which baseline?** `raw_null` (the exploration's recommendation) or
   `stored` (what 1.2.0 ships). *Recommendation: `raw_null`*, so that the
   run's experiments are spent beyond the known candidates; the gain over
   `stored` is read from the reference rows.

6. **What is a realistic donor SD?** The exploration used 0.1 throughout,
   and the answer changes how much of the sampling variance is between
   donors. *Recommendation: before the build, a ten-minute check on the
   fitted ASD object in `~/Downloads/esvd_asd_tmp/`* (the SD across donors
   of the log donor mean of `A / library`, per gene, less its Poisson
   part), with the two levels set at its lower and upper quartiles.

7. **600 genes at the true rank.** `locfdr` on 300 genes often falls back
   to its cruder estimators, which handicaps every null-based direction
   relative to the real 6,599 genes; 600 doubles the build time. The
   fitted rank is fixed at the truth because `overfit_k` moved no
   candidate. *Recommendation: 600 genes, rank fixed.*

8. **The truth under one-directional DE.** In `sparse_up` the truth is the
   generator's coefficient, an absolute change at equal depth, which no
   method can identify without assuming that most genes are unchanged.
   With 10% of genes DE that assumption holds. *Recommendation: keep it,
   and do not add a cell in which most genes move one way.*

9. **Covariates as a two-level factor** doubles the bank for a question
   the exploration only flagged. *Recommendation: keep it*; the fitted
   library that `raw` divides by is where covariates enter, and it is the
   least-tested part of the baseline.

10. **Cores, hours, cap**: 8 cores, 4 hours, 40 experiments.
    *Recommendation: as stated*, confirmed again at the run's setup.

## Vocabulary

| Term | Meaning |
|---|---|
| reporting rule | the function on trial, `report_lfc_se(view)` |
| score, NLL | Gaussian negative log-likelihood of the true log2FC, nats per gene, lower is better |
| example | one pair of fitted cohorts (A, B) |
| `donors`, `cells`, `rates`, `effects`, `covariates`, `donor_sd` | the six factors of the bank |
| `sparse`, `sparse_up`, `dense` | the levels of `effects` |
| gate, `slope_weak` | the slope of the logFC on the truth over weak-effect genes, required in [0.95, 1.05] |
| `rmse`, `coverage_95`, `z_sd_null`, `bias_null`, `wald_false`, `k_err`, `r_err`, `rho_err`, `score_vs_deseq2` | the diagnostics of the summary block |
| `stored`, `null_rescaled`, `between`, `between_plus`, `jackknife`, `deshrunk`, `raw`, `raw_null`, `deseq2` | the exploration's candidates; `raw_null` is the baseline, `deseq2` the reference |
| D1 to D5 | donor summary; variance form; centring; borrowing across genes; using more of the fit |
| `lfcse-bank-1` | the bank and its generator version |
| `autoresearch_lfcse/<mmmdd>` | a run's branch |
