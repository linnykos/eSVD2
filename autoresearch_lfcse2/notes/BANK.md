# The bank `lfcse-bank-2`, as built

Kevin Z. Lin (drafted by Claude), 2026-11-11. The design is section 3 of
`MEMO_claude.md`; this file records what the build produced.

## Where and how

| | |
|---|---|
| Home | `LFCSE_BANK`: `<ESVD2_BANK_DIR>/lfcse-bank-2/`, default root `tools::R_user_dir("eSVD2", "cache")`. Outside the repository |
| Generator | `harness/generator_claude.R`, version `lfcse-bank-2` |
| Builder | `VECLIB_MAXIMUM_THREADS=1 Rscript autoresearch_lfcse2/harness/prepare_claude.R` |
| Built | 2026-11-11, eSVD2 1.2.0 (`additional_context/version_comparison/lib/devel`), R 4.5.1 |
| Cost | 36 minutes on 8 cores (4.5 core-hours of fitting over 576 cohorts); 628 MB on disk |
| Examples | 288 dev and 288 test, one cohort each, every fit succeeded; the rebuilt posterior equals the fit's for all 576 (largest difference 0) |
| A scoring | 17 to 20 seconds on 8 cores; nothing is refitted |

One example is two files: `<split>/<id>_public.rds` (the fitted cohort:
counts, design, donor labels, the fit, the nuisance rates and their status,
and what 1.2.0 reports per gene) and `<split>/<id>_truth.rds` (the true
log2 fold changes, gene classes, the population, each donor's state and
whole effect on every gene, DESeq2's pseudobulk triple, `donor_wls`'s pair,
the factor labels).

## Realized against target

| Quantity | Target (memo) | Realized |
|---|---|---|
| Cells per cohort | 500 to about 9,000 | 500 to 9,683 |
| Genes capped or at the boundary, `expression = high` | not set (first bank: 13% at its wide rates) | 31% (39% at 10 donors, 25% at 40) |
| Genes capped or at the boundary, `expression = low` | about half; the real cohorts have 73% | **49%** (52% at 10 donors, 46% at 40) |
| Mean count per cell and gene | | 3.15 at `high`, 0.37 at `low` |
| Share of the donor variance carried by the states | 75% | about 72% (the idiosyncratic part is scaled by the gene's own donor SD); the total donor variance at `donor` and `donor_arm` is about 12% below that at `none` |
| Exactly null genes | 540 of 600 at `sparse` and `sparse_up` | 0 at `dense` and at `state = donor_arm`, where every gene has an effect |
| Fit warnings | none expected | 4 cohorts of 576 |

The capped share is higher than the first bank's at both levels (rank 5
and the donor-level states both raise it) and is still below the real
cohorts' 73%.

**Aliasing in the design.** Within a split every level and every pair of
levels of two factors is balanced, but `donors` x `effects` x `state` is a
third fraction: at a given `state`, each `donors` level occurs with exactly
one `effects` level. A breakdown by `effects` within one `state` is also a
breakdown by `donors`.

## The second gate

`se_strong` on dev, for the reference rules and for the rule the first run
withdrew:

| Rule | `se_strong` |
|---|---|
| `deshrunk` | 0.975 |
| `donor_wls` | 1.020 |
| `null_rescaled`, `stored` | 1.043 |
| `gap_trend` (withdrawn in the first run) | **1.203** |

The interval is [0.92, 1.08]: a margin of about 0.04 beyond the valid
rules on each side, and 0.12 short of the invalid one. A first version of
the gate classed genes by the marginal truth, which at `state = donor_arm`
selects genes with large loadings and so with a truly larger donor
variance; the fresh-context review found it, and the gate now classes
genes by the direct case-control effect.

## What the embedding could give: oracle rows (dev)

`donor_wls` with a donor x gene offset subtracted from the log donor means
(`diagnostics/bank_checks_claude.R`). `oracle_state` subtracts the donor's
true state effect and block shift, which is what an ideal use of the
embedding could remove; `oracle_all` subtracts the donor's whole true
effect and leaves sampling noise only.

| `state` | `donor_wls` | `oracle_state` | `oracle_all` |
|---|---|---|---|
| `none` | 0.184 | 0.148 | -0.748 |
| `donor` | 0.155 | -0.210 | -0.751 |
| `donor_arm` | 0.114 | -0.222 | -0.733 |
| all | 0.151 | -0.094 | -0.744 |

The states hold about 0.35 nats per gene at the two levels that have them,
and 0.04 at `none` through the block shift. The baseline is 0.16 behind
`donor_wls`, so a rule that used the embedding ideally would gain about
0.4 nats on the baseline.

## Baseline (`null_rescaled`)

Scored on 8 cores before the loop folder's first commit; the per-example
tables are `output/baseline_dev_scores.csv` and
`output/baseline_test_scores.csv`.

```
                            dev                 test
score:                      0.3066              0.2969
score_se:                   0.0293              0.0291
rmse:                       0.3780              0.3802
coverage_95:                0.9339              0.9371
z_sd_null:                  1.1421              1.1548
bias_null:                  -0.0048             -0.0098
slope_weak:                 0.9665              0.9468
wald_false:                 2.247               2.370
se_strong:                  1.0431              1.0416
score_vs_deseq2:            0.0993              0.1028
score_vs_donor_wls:         0.1558              0.1650
score_donors_10:            0.6374              0.6126
score_donors_20:            0.2767              0.3086
score_donors_40:            0.0057              -0.0304
score_cells_equal:          0.3026              0.3333
score_cells_unequal:        0.3106              0.2605
score_expression_high:      0.1041              0.1173
score_expression_low:       0.5091              0.4765
score_effects_dense:        0.3473              0.2659
score_effects_sparse:       0.2616              0.3252
score_effects_sparse_up:    0.3109              0.2997
score_covariates_none:      0.2684              0.2751
score_covariates_sexage:    0.3448              0.3188
score_donor_sd_0.2:         0.0293              0.0159
score_donor_sd_0.5:         0.5839              0.5779
score_state_donor:          0.2330              0.2278
score_state_donor_arm:      0.3117              0.2966
score_state_none:           0.3751              0.3664
score_arm_sd_equal:         0.3121              0.3231
score_arm_sd_unequal:       0.3011              0.2708
score_depth_error:          0.3526              0.3686
score_depth_none:           0.2606              0.2252
score_status_capped:        0.2360              0.2195
score_status_estimated:     0.3387              0.3366
seconds:                    17.1                17.3
n_examples:                 288                 288
status:                     ok                  ok
gate:                       pass                pass
```

The donor count, the donor SD and the expression level move the score
most (0.63, 0.55 and 0.41 nats between levels). The baseline is 0.10
behind DESeq2 and 0.16 behind `donor_wls` on both splits, and its Wald
statistic is too wide on exactly null genes (SD 1.14).

## Reference rows (dev, scored once at the build)

Lower is better. Each is a file of `candidates/`, scored with
`ESVD2_TRIAL=<file>`. The last three use nothing of the embedding and are
not candidates of this loop.

| Rule | Uses the embedding | Score | RMSE | Coverage | `z_sd_null` | `slope_weak` | `se_strong` | 10 donors | 40 donors | `state = donor` |
|---|---|---|---|---|---|---|---|---|---|---|
| **`null_rescaled`** (baseline) | through the posterior | 0.307 | 0.378 | 0.934 | 1.14 | 0.967 | 1.043 | 0.637 | 0.006 | 0.233 |
| `stored` (1.2.0) | through the posterior | 0.468 | 0.386 | 0.984 | 0.65 | 0.963 | 1.043 | 0.816 | 0.133 | 0.457 |
| `deshrunk` | through the posterior | 0.577 | 0.386 | 0.957 | 0.80 | 0.963 | 0.975 | 0.999 | 0.201 | 0.620 |
| `donor_wls` (first run's kept rule) | no | 0.151 | 0.372 | 0.936 | 1.05 | 0.989 | 1.020 | 0.518 | -0.197 | 0.155 |
| `deseq2` | no | 0.207 | | | | | | | | |
| `gap_trend` | no | 0.138, gate failed | 0.373 | 0.941 | 1.05 | 0.996 | 1.203 | 0.479 | -0.188 | 0.142 |

`raw_null`, the first loop's baseline, returns a non-finite value for one
gene in four of the sparse cohorts and was not scored.

## Review

A fresh-context review of the harness (2026-11-11) confirmed the design's
balance, the seed scheme, the truth by a Monte Carlo through the generator
(no gene off beyond Monte Carlo error at any `state`, `arm_sd` or `depth`
level), the score against an independent implementation, and that the
public files carry nothing of the truth. It found two things that were
fixed in the scorer before the baseline: the trial function's random seed
coincided with the stream that drew the true effects (it is now an offset
no generator stream uses; the first loop's scorer has the same collision
and was not changed), and the second gate's class of strong genes (above).
