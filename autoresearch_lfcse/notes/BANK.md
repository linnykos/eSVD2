# The bank `lfcse-bank-1`, as built

Kevin Z. Lin (drafted by Claude), 2026-11-11. The design is section 3 of
`MEMO_claude.md`; this file records what the build produced.

## Where and how

| | |
|---|---|
| Home | `LFCSE_BANK`: `<ESVD2_BANK_DIR>/lfcse-bank-1/`, default root `tools::R_user_dir("eSVD2", "cache")`. Outside the repository |
| Generator | `harness/generator_claude.R`, version `lfcse-bank-1` |
| Builder | `VECLIB_MAXIMUM_THREADS=1 Rscript autoresearch_lfcse/harness/prepare_claude.R` |
| Built | 2026-11-11, eSVD2 1.2.0 (`additional_context/version_comparison/lib/devel`), R 4.5.1 |
| Cost | 47 minutes on 5 cores (3.6 core-hours of fitting over 576 cohorts); 699 MB on disk |
| Examples | 144 dev and 144 test, one per cell of the 3 x 2 x 2 x 3 x 2 x 2 factorial per split; every fit succeeded |
| A scoring | 23 seconds on 5 cores; nothing is refitted |

One example is two files: `<split>/<id>_public.rds` (the two fitted cohorts:
counts, design, donor labels, the fit, the nuisance rates and their status,
and what 1.2.0 reports per gene) and `<split>/<id>_truth.rds` (the true
log2 fold changes, gene classes, the population, DESeq2's pseudobulk
triple, the factor labels). The public file holds no cells x genes matrix
of doubles; the scorer rebuilds the library and the posterior with
`.build_view()`, and the builder checked for all 576 cohorts that the
rebuilt posterior equals the fit's (largest difference 0).

A rebuild reproduces the files: the seed of an example is
`1000000 * split + 100 * cell`, and the reviewer regenerated one cohort
from its seed and obtained identical counts.

## Realized against target

| Quantity | Target (memo) | Realized |
|---|---|---|
| Donor SD levels | lower and upper quartiles on the ASD fit | 0.2 and 0.5 (quartiles 0.21 and 0.52; `diagnostics/donor_sd_asd_claude.R`) |
| Cells per donor, `unequal` | median 150, range 30 to 600 | smallest donor 30, largest 600; cohorts of 1,025 to 9,213 cells |
| Sex among cases / controls | 60% / 40% | exactly 0.6 / 0.4 at every donor count |
| Population `k` | 2 | median 1.99 to 2.01 per `effects` level (range 1.62 to 2.48) |
| Population `r` | 0.8 | median 0.80 for `sparse` and `dense`; 0.85 for `sparse_up`, where taking absolute values raises the uncentred correlation. The scorer uses each example's realized `k` and `r` |
| Genes capped, `rates` narrow | about 11% | mean 0.8% (3% at 10 donors x 50 cells, under 1% elsewhere) |
| Genes capped, `rates` wide | about half | **mean 13%** (median 13%, range 0 to 39%; 15 to 18% at 10 donors, 11% at 40), plus 1.3% at the boundary |
| Reliability of cohort A under the baseline, `dense` | 0.4 to 0.65 | median 0.46 on dev, 0.44 on test (0.32 at 10 donors, 0.64 at 40) |
| Reliability, `sparse` / `sparse_up` | not set | median 0.58 / 0.60 on dev |
| Fit warnings | none expected | 3 cohorts of 576, each `locfdr` falling back to the truncated MLE in `compute_pvalue()`; all three at donor SD 0.5 |

**The capped share is the one target the bank misses.** The memo's "about
half" came from cohorts of 600 cells and the real cohorts sit at 73%; here
a cohort has 500 to 9,000 cells and the cap binds less as cells are added.
The `score_status_capped` line therefore rests on few genes (the dev and
test baselines give -0.01 and +0.27), and a direction whose case is the
capped genes (D5, queue item 8) is weakly tested by this bank. Raising the
share needs a different rate range or a tighter cap, which is a new bank.
Kevin accepted the 13% for this bank (2026-11-11).

**The Figure 4 diagnostics are defined for 80% of pairs.** In the others
the estimated signal `mean(lfc^2) - mean(se^2)` of one cohort is not
positive, nearly all at donor SD 0.5, where cohort B's effects (half of
A's) are below the noise. `k_err` and `r_err` average over the valid pairs
and are dominated by a few near-zero denominators (median `k_err` 0.10 and
`r_err` 0.05 at donor SD 0.2); read them by stratum in the per-example
table, not from the block.

## Baseline (`raw_null`)

Scored on 5 cores at commit `5f8bc61` plus the uncommitted loop folder,
trial file md5 `ef2dcab387e893514794e167eb5e0d29`. The per-example tables
are `output/baseline_dev_scores.csv` and `output/baseline_test_scores.csv`.

```
                            dev        test
score:                      0.0097     0.0178
score_se:                   0.0433     0.0440
rmse:                       0.2937     0.2972
coverage_95:                0.9484     0.9484
z_sd_null:                  1.0269     1.0295
bias_null:                  -0.0085    -0.0079
slope_weak:                 0.9843     0.9741
wald_false:                 0.567      0.592
k_err:                      0.5677     0.4999
r_err:                      1.4330     0.7122
figure4_valid:              0.799      0.806
rho_err:                    0.1427     0.1310
score_vs_deseq2:            0.0340     0.0350
score_donors_10:            0.3973     0.4181
score_donors_20:            -0.0116    -0.0044
score_donors_40:            -0.3566    -0.3603
score_cells_equal:          0.0273     0.0390
score_cells_unequal:        -0.0079    -0.0034
score_rates_narrow:         0.0083     0.0243
score_rates_wide:           0.0111     0.0113
score_effects_dense:        0.0279     0.0351
score_effects_sparse:       -0.0042    0.0058
score_effects_sparse_up:    0.0054     0.0124
score_covariates_none:      -0.0243    -0.0181
score_covariates_sexage:    0.0437     0.0537
score_donor_sd_0.2:         -0.4008    -0.3932
score_donor_sd_0.5:         0.4202     0.4287
score_status_capped:        -0.0107    0.2747
score_status_estimated:     0.0097     0.0179
gate:                       pass       pass
seconds:                    22.7       22.4
n_examples:                 144        144
status:                     ok         ok
```

The donor count and the donor SD move the score by 0.75 to 0.8 nats each;
every other factor moves it by under 0.07. The gap to DESeq2 is 0.034
overall and is largest at 10 donors with donor SD 0.2 (0.084); at 40
donors with donor SD 0.5 the baseline matches DESeq2 (-0.001).

## Reference rows (dev, scored once at the build)

Lower is better. Each is a file of `candidates/`, scored with
`ESVD2_TRIAL=<file>`.

| Candidate | Score | SE | RMSE | Coverage | `z_sd_null` | `bias_null` | `slope_weak` | `wald_false` | vs DESeq2 | 10 donors | 40 donors | `dense` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `null_rescaled` | -0.006 | 0.042 | 0.296 | 0.950 | 1.02 | -0.009 | 0.983 | 0.44 | +0.018 | 0.339 | -0.344 | 0.008 |
| **`raw_null`** (baseline) | 0.010 | 0.043 | 0.294 | 0.948 | 1.03 | -0.009 | 0.984 | 0.57 | +0.034 | 0.397 | -0.357 | 0.028 |
| `deshrunk` | 0.083 | 0.045 | 0.299 | 0.906 | 1.15 | -0.043 | 0.962 | 2.06 | +0.107 | 0.514 | -0.310 | 0.039 |
| `raw` | 0.099 | 0.047 | 0.297 | 0.894 | 1.20 | -0.042 | 0.961 | 2.59 | +0.124 | 0.568 | -0.317 | 0.055 |
| `stored` (1.2.0) | 0.223 | 0.029 | 0.299 | 0.986 | 0.65 | -0.043 | 0.962 | 0.04 | +0.247 | 0.563 | -0.116 | 0.211 |

The order of the exploration holds on this bank, with smaller gaps: the
stored SE costs 0.21 nats against the baseline (0.82 on the exploration's
cohorts, whose donor SD was 0.1), and `null_rescaled` is again ahead of the
baseline at 10 donors (0.058) and behind it at 40 (0.012). All five pass
the gate.

## Review

A fresh-context review of the harness (2026-11-11) confirmed the
generator's math, the seed scheme, the score against an independent
implementation, and that the public files carry nothing of the truth. It
found that a token scan alone let a trial file reach the scorer's call
stack (`dynGet`) or build a forbidden name at run time
(`lapply(x, paste0("read", "RDS"))`). The scorer now loads the trial file
into an environment whose only ancestor is a curated copy of base R with
its own `::`, seeds the child process per example and cohort, and the
builder writes each file under a temporary name. The bypasses the review
listed were rerun and are refused. It remains a tripwire: a run report
lists every constant a kept commit introduced.
