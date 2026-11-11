# Exploration: a Wald-usable log fold change and standard error for eSVD-DE

Kevin Z. Lin (drafted by Claude), 2026-11-07. Branch `2026-11-07_std`.

## Why

Writeup50's Figure 4 (Was2CoDE, `Writeup50_figure4-from-summary-statistics_claude.Rmd`)
needs, per gene and cohort, a triple (logFC, SE, p) in which `SE^2` is the
sampling variance of logFC across cohort replications, so that
`mean(logFC^2) = mean(beta^2) + mean(SE^2)` across genes. eSVD-DE 1.2.0 stores
an SE whose dominant term is the average per-cell posterior variance (the
"within" term), which is 2 to 3 times the Welch denominator's calibrated size;
the empirical null (`locfdr`) rescales the test statistic but neither the
logFC nor the SE. `Writeup50_esvd-wald-se_claude.md` lists options A to E.

Kevin's additional question (2026-11-07): is eSVD-DE's lower sign agreement
in Figure `Writeup50_sign-agreement.png` (the pairs with Prater, about 60%)
caused by not correcting the logFC for the empirical null, as the p-value is?

## What is here

| File | What it does |
|---|---|
| `helpers_claude.R` | Generator (population truth shared across replicate cohorts; donor random effects; Gamma-Poisson counts), the eSVD 1.2.0 pipeline in the order `eSVD()` runs it, every candidate (logFC, SE) computed from one fitted object, a DESeq2 pseudobulk reference, and the Figure 4 moments decomposition copied from the Writeup50 Rmd. |
| `01_real-data-checks_claude.R` | Tati's four cohorts from the Writeup50 CSVs and her diet `eSVD` objects: sign agreement under logFC recentering, same-gene cross-method sign comparison, and the empirical null fitted within nuisance-status strata. No refit. |
| `02_replicate-cohorts_claude.R` | The oracle study: `R` replicate cohorts drawn from one population truth, each refitted, so the SD of logFC across replicates is the sampling SE Figure 4 assumes. Every candidate SE is compared with it, by nuisance status. |
| `03_figure4-pairs_claude.R` | Pairs of cohorts (A, B) with known `k` and `r`; the moments decomposition under each candidate against the truth, with DESeq2 as the positive control. |
| `04_score-demo_claude.R` | Added 2026-11-11 for the autoresearch memo (`../AUTORESEARCH_LFC-SE_MEMO_claude.md`): the Gaussian negative log-likelihood of the truth for every candidate on the cached fits of 02 and 03, with a thresholding rule as the exploit to guard against. No fitting; not run by `run_all_claude.sh`. |
| `summary_claude.Rmd` | The write-up. Reads only `output/*.csv`. |
| `run_all_claude.sh` | Runs 01 to 03 and knits the summary. |

`cache/` (per-replicate fits) and `output/*.csv` are `.gitignore`d; the
knitted `summary_claude.html` is tracked.

## Candidate triples (all log2 scale)

| Name | logFC | SE |
|---|---|---|
| `stored` | `log2(m1/m0)`, ratio of arm means of per-donor posterior means | 1.2.0's `log2fc_se_vec` (mixture variance: within + between) |
| `null_rescaled` (option A) | `stored - null_mean * se_stored` | `null_sd * se_stored` |
| `between` (option B) | `stored` | delta method with the between-donor sample variance only (`n - 1` divisor) |
| `between_plus` | `stored` | `between` with the within term divided by the donor's cell count added |
| `jackknife` (option D, fit fixed) | `stored` | leave-one-donor-out SD of the logFC |
| `raw` | `log2` ratio of arm means of per-donor means of `A / library` (unshrunk) | between-donor, delta method |
| `deseq2` | DESeq2 on pseudobulk, `~ CC` | `lfcSE` |

The oracle is the SD of each candidate's logFC over the replicate cohorts.

## Findings (2026-11-07, first pass)

The write-up is `summary_claude.html`. In short: recentering the logFC does
not change the sign agreement on Tati's cohorts (the non-replicating eSVD-DE
DEGs are its smallest effects); the stored SE is 3.4 to 5 times the
replicate-cohort SD and the between-donor SE of option B is 0.65 to 0.73 of
it, shrunk by the posterior's weight on the fit; the SE of the unshrunk
per-donor means (`raw`), or the between-donor SE divided by the shrinkage
weight, is calibrated and recovers Figure 4's k and r; the logFC should be
shifted by the empirical null's mean, which carries the depth adjustment's
bias (large only when the DE is one-directional).

## Running

```sh
cd additional_context/2026-11-07_std-exploration
bash run_all_claude.sh          # 01, 02, 03, knit; several hours on a laptop
```

Each script loads eSVD2 1.2.0 from `additional_context/version_comparison/lib/devel`
(built 2026-09-30 from the working tree; `R/` has not changed since).
