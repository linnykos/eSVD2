# Does `eSVD2` need `Rmpfr`? — decision experiment

**Verdict: no. Drop it from `DESCRIPTION` entirely.**

Run 2026-08-29 on R 4.5.1 / macOS (Darwin 23.2.0), `Rmpfr` 1.1.2, GMP 64-bit
limbs. Reproduce with:

```
Rscript additional_context/rmpfr_experiment_claude.R
```

This is the one-off experiment called for by `UNIT_TEST_PLAN.md` §2.15 and
question **Q-MPFR-1** ("let's do this outside of the package"). Tests T-MPFR-01
through T-MPFR-04 and T-MPFR-07 live only in that script, because they require
the very dependency they exist to remove. Their shipped counterparts —
T-MPFR-05, T-MPFR-06, T-MPFR-08 in `tests/testthat/test_tail_precision.R` — need
no `Rmpfr` and pin the same conclusion permanently.

## Where `Rmpfr` is used today

Eight calls, all `Rmpfr::pnorm`, all on plain doubles:

| File | Lines | `log.p` | Feeds |
|---|---|---|---|
| `R/multtest.R` | 26, 28 | `TRUE` | `logpvalue_vec` |
| `R/multtest.R` | 36, 38 | **`FALSE`** | `pvalue_vec` → `stats::p.adjust()` |
| `R/compute_pvalue.R` | 87, 89 | `TRUE` | `log10pvalue` |
| `R/compute_test_per_gene.R` | 329, 334 | `TRUE` | `log10pvalue` |

## The four findings

### 1. The calls as written are literal no-ops (T-MPFR-01)

`Rmpfr::pnorm` is an S4 generic. Given a plain `numeric` it dispatches straight
back to `stats::pnorm`.

```
class(Rmpfr::pnorm(<double>)) = numeric
identical(stats::pnorm(z, log.p=TRUE), Rmpfr::pnorm(z, log.p=TRUE)) = TRUE
```

`TRUE` at every point of the grid `z ∈ {−5, −10, −20, −37, −38.5, −40, −100,
−1000, −1e4}`. No argument anywhere in `eSVD2` is ever an `mpfr` object, so the
extra precision is never engaged.

**This much was already in `CRAN_READINESS.md` §2.2.** It shows the current calls
buy nothing; it does not show that arbitrary precision would be useless if used
*correctly*. The next three findings do.

### 2. Doubles are already correct to the last bit, in log space (T-MPFR-02)

Against a genuine 200-bit MPFR oracle (`Rmpfr::pnorm(Rmpfr::mpfr(z, 200))`):

| `z` | `stats::pnorm(z, log.p=TRUE)` | relative error vs. 200-bit |
|---|---|---|
| −5 | −1.506500e+01 | 2.82e−17 |
| −10 | −5.323129e+01 | 1.49e−17 |
| −20 | −2.039172e+02 | 9.97e−18 |
| −37 | −6.890306e+02 | 7.32e−17 |
| −38.5 | −7.456953e+02 | 6.74e−17 |
| −40 | −8.046084e+02 | 1.69e−17 |
| −100 | −5.005524e+03 | 4.21e−17 |
| −1000 | −5.000078e+05 | 4.86e−17 |
| −10000 | −5.000001e+07 | 4.08e−17 |

**Maximum relative error 7.3e−17** — below one unit in the last place of a
double (2.2e−16). The accuracy does not degrade as `z` grows: the error at
`z = −10⁴` is no worse than at `z = −5`. There is no regime `eSVD2` can reach
where a double's log-tail is wrong.

### 3. `Rmpfr` does not prevent the one real underflow (T-MPFR-03)

`multtest()` computes `pvalue_vec` with `log.p = FALSE`. That branch underflows,
and it underflows *identically* with and without `Rmpfr`:

```
2*stats::pnorm(-38.5) = 0     2*stats::pnorm(-40) = 0
2*Rmpfr::pnorm(-38.5) = 0     2*Rmpfr::pnorm(-40) = 0
```

The zero p-values that reach `stats::p.adjust()` are caused by the
**parameterization**, not by the precision. `Rmpfr` is positioned as if it were
guarding this and is not.

### 4. MPFR precision cannot reach the user anyway (T-MPFR-04)

This is the finding that closes the question. Suppose the calls were rewritten
to use MPFR properly — carrying `mpfr` objects rather than doubles. The final
step of the pipeline is `stats::p.adjust()`, which is double-only and coerces
silently:

```
input  (200-bit mpfr):  7.3117870818300594e-350
                        9.8134278542963741e-198
                        0.0455002638963584144

p.adjust(., "BH"):      0.000000e+00
                        1.472014e-197
                        4.550026e-02
```

The `1e−350` entry returns as **exactly 0**. Delivering MPFR precision to the
user would require reimplementing Benjamini–Hochberg in `mpfr` — for a
difference that changes no rejection at any FDR threshold anyone uses.

## What *is* the fix

**`log.p = TRUE`, which is `CRAN_READINESS.md` §1.1.** Confirmed here (T-MPFR-07):

```
qnorm(pt(40, 18))                          = Inf
qnorm(pt(-40, 18, log.p=TRUE), log.p=TRUE) = -8.9152934
```

and the resulting `log10pvalue` computed in double agrees with the same
computation carried out entirely in 200-bit MPFR to **relative difference 0**.
Two lines of implementation replace a GMP/MPFR system dependency.

## What *is* genuinely lost, and it is not precision (T-MPFR-08)

`report_results()` returns `pvalue = 10^(-log10pvalue)`, which underflows to `0`
above `log10pvalue = 308`:

| `log10pvalue` | reported `pvalue` |
|---|---|
| 400 | 0 |
| 500 | 0 |
| 800 | 0 |

Three genes differing by 400 orders of magnitude report the same `p = 0` and
**cannot be ranked** — while `pvalue_list$log10pvalue` has held the distinction
at full double precision the entire time. More bits would not help: the column is
a double either way.

**Fix: add a `log10pvalue` column to `report_results()`.** One line. Recommended
alongside the §1.1 change, since both concern the same extreme-tail regime.

## The permanent, dependency-free replacement (T-MPFR-05)

`stats::pnorm(z, log.p = TRUE)` can be checked against the Mills-ratio
asymptotic expansion of the normal log-tail, which needs no package at all:

```
log Φ(−z) = −z²/2 − log z − log(2π)/2 + log(1 − z⁻² + 3z⁻⁴ − 15z⁻⁶ + 105z⁻⁸)
```

| `z` | relative error vs. `stats::pnorm` |
|---|---|
| −10 | 1.62e−09 |
| −20 | 4.42e−13 |
| −40 | 1.41e−16 |
| −100 | 0 |
| −1000 | 1.16e−16 |
| −10000 | 0 |

The series is asymptotic, so it is loose at `z = −10` and tightens from there;
the shipped test asserts agreement for `|z| ≥ 20` at relative 1e−12. That pins
double-precision adequacy for good, with nothing in `Suggests`.

## Recommended actions

1. Replace all eight `Rmpfr::pnorm` calls with `stats::pnorm`.
2. Remove `Rmpfr` from `Imports:` in `DESCRIPTION`. Do **not** move it to
   `Suggests:` — nothing that ships needs it.
3. Apply the `log.p = TRUE` fix of §1.1 in the same commit; it touches the same
   call sites.
4. Add a `log10pvalue` column to `report_results()`.
5. Drop the `Rmpfr` install caveat from `README.md` (it currently warns that
   `Rmpfr` "is sometimes tricky to install due to its required C++ libraries" —
   which is true, and is now moot).
6. No `SystemRequirements:` entry is needed for GMP/MPFR, since the dependency
   is gone.

## Raw output

The full console transcript is reproducible with the command at the top of this
file. The tables above are transcribed from it verbatim.
