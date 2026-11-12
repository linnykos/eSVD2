# autoresearch_lfcse2

The second autoresearch loop on the **reporting rule**: a function that
takes one fitted eSVD-DE cohort and returns, for every gene, a log2 fold
change and a standard error usable as a Wald pair. It differs from
`autoresearch_lfcse/` in two ways that Kevin set on 2026-11-11: the rule
must use the fitted embedding (explicitly, or through the posterior), and
the bank contains the settings the first run's report listed as untested,
above all donor-level cell states that the embedding can recover. The
design and every decision behind it are in `notes/MEMO_claude.md`; the bank
as built is in `notes/BANK.md`.

Nothing here ships to CRAN (`.Rbuildignore`). All commands run from the
repository root.

## Layout

```
autoresearch_lfcse2/
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

**The embedding rule.** The logFC or the SE must be computed from the
fitted embedding: from `x_mat` and `y_mat`, or from the posterior mean or
variance, which contain them. A rule that reads only the counts, the
design and the fitted library is out of bounds, and so is a rule in which
the embedding enters only as a token term that could be deleted without
changing the method.

## The scorer

```
Rscript autoresearch_lfcse2/harness/evaluate_claude.R [dev|test] > autoresearch_lfcse2/output/run.log 2>&1
```

It checks the bank against its manifest, runs the tripwire on the trial
file, and for each example rebuilds the cohort's view, calls the trial
function in a child process, checks the returned value, and only then
reads the truth. The score of an example is the Gaussian negative
log-likelihood of the true log2 fold change under the reported pair, in
nats per gene; **lower is better**. The log ends with the summary block:
`score`, `score_se`, the diagnostics, `score_vs_deseq2` and
`score_vs_donor_wls` (the two references stored with the truth), one
`score_<factor>_<level>` line per stratum, `delta_vs_best` and `paired_se`
against the best commit, the `gate`, and `status`.

The gate has two parts. `slope_weak`, the slope of the reported logFC on
the truth over weak-effect genes, must lie in [0.9, 1.1]. `se_strong`, the
SE of strong-effect genes over the SE of genes of the same expression level
without a strong effect, must lie in the interval `BANK.md` records, so
that an SE which follows the estimated logFC is refused.

| Variable | Default | Meaning |
|---|---|---|
| `ESVD2_BANK_DIR` | `tools::R_user_dir("eSVD2", "cache")` | root of the bank |
| `ESVD2_CORES` | 8 | examples scored in parallel |
| `ESVD2_BUDGET` | 600 | wall-clock seconds of one scoring |
| `ESVD2_TRIAL` | unset | score another file (a candidate) without touching a run's files |

## The bank

```
VECLIB_MAXIMUM_THREADS=1 Rscript autoresearch_lfcse2/harness/prepare_claude.R [--resume]
```

builds `lfcse-bank-2` under the bank root: 288 dev and 288 test examples,
each one cohort generated and fitted end to end by eSVD2 1.2.0, so a
scoring refits nothing. It uses `ESVD2_BUILD_CORES` cores (default 8).

## Logging an experiment

```
bash autoresearch_lfcse2/harness/record_result.sh <keep|discard|gate_failed|crash|timeout> "<description>" [--tie] [--killed]
```
