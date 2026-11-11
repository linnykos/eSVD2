# autoresearch_lfcse

An autoresearch loop that searches for a **reporting rule**: a function that
takes one fitted eSVD-DE cohort and returns, for every gene, a log2 fold
change and a standard error usable as a Wald pair. One file is on trial, one
scorer judges every version of it on a fixed bank of fitted cohorts, and
one script logs each version and keeps or discards it. The design and every
decision behind it are in `notes/MEMO_claude.md`; the bank as built is in
`notes/BANK.md`.

Nothing here ships to CRAN (`.Rbuildignore`). All commands run from the
repository root.

## Layout

```
autoresearch_lfcse/
  README.md                 this page
  report_rule_claude.R      the ONE file an experiment edits
  harness/                  fixed machinery; read-only during a run
    generator_claude.R      bank location, design, generator, the eSVD-DE fit, the public view
    prepare_claude.R        builds the bank outside the repository
    evaluate_claude.R       the scorer
    record_result.sh        the only writer of output/results.tsv
    loop.conf               the four keys record_result.sh reads
  candidates/               drop-in bodies for the trial function (the reference rows of BANK.md)
  diagnostics/              truth-reading scripts; their code never enters the trial file
  notes/                    MEMO_claude.md, BANK.md, and each run's report, rule file and figures
  output/                   gitignored: logs, per-example tables, results.tsv
```

## The trial function

`report_rule_claude.R` defines `report_lfc_se(view)`.

`view` is what a user has after `eSVD()` with `bool_diet = FALSE`: `dat`
(counts, cells x genes), `covariates`, `individual_vec`,
`case_control_vec`, the fit `x_mat`, `y_mat`, `z_mat`, `library_mat`,
`nuisance_vec` and `nuisance_status_vec`, `posterior_mean_mat`,
`posterior_var_mat`, and `stored_df` (what 1.2.0 reports per gene, with the
fitted `null_mean` and `null_sd` as attributes). Nothing in it comes from
the truth.

It returns a data frame with one row per gene in the order of
`colnames(view$dat)` and columns `gene`, `log2fc`, `se`, every value finite
and every `se` positive.

It may call base R, `stats`, `Matrix`, `locfdr` and eSVD2's exported
functions, each with `::`. The file's top level holds function definitions
and named constants only. It reads and writes no file, and has 30 seconds
per cohort. Every cohort is treated by the same formula, every numeric
constant carries the reason for its value, and the logFC is not shrunk
toward zero (a shift common to all genes is allowed).

## The scorer

```
Rscript autoresearch_lfcse/harness/evaluate_claude.R [dev|test] > autoresearch_lfcse/output/run.log 2>&1
```

It checks the bank against its manifest, runs the tripwire on the trial
file, and for each example rebuilds the two cohorts' views, calls the trial
function in a child process, checks the returned value, and only then
reads the truth. The score of an example is the Gaussian negative
log-likelihood of the true log2 fold change under the reported pair, in
nats per gene, averaged over its two cohorts; **lower is better**. The log
ends with the summary block: `score`, `score_se`, the diagnostics, one
`score_<factor>_<level>` line per stratum, `delta_vs_best` and `paired_se`
against the best commit, the `gate` (the slope of the reported logFC on the
truth over weak-effect genes must lie in [0.9, 1.1]), and `status`.

| Variable | Default | Meaning |
|---|---|---|
| `ESVD2_BANK_DIR` | `tools::R_user_dir("eSVD2", "cache")` | root of the bank |
| `ESVD2_CORES` | 8 | examples scored in parallel |
| `ESVD2_BUDGET` | 600 | wall-clock seconds of one scoring |
| `ESVD2_TRIAL` | unset | score another file (a candidate) without touching a run's files |

## The bank

```
VECLIB_MAXIMUM_THREADS=1 Rscript autoresearch_lfcse/harness/prepare_claude.R [--resume]
```

builds `lfcse-bank-1` under the bank root: 144 dev and 144 test examples,
each a pair of cohorts generated and fitted end to end by eSVD2 1.2.0, so a
scoring refits nothing. It uses `ESVD2_BUILD_CORES` cores (default 5).

## Logging an experiment

```
bash autoresearch_lfcse/harness/record_result.sh <keep|discard|gate_failed|crash|timeout> "<description>" [--tie] [--killed]
```
