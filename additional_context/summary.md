# additional_context/ — Paper Summaries

This folder contains reference PDFs for **eSVD2**, the R package implementing
eSVD-DE. Read this file before opening any PDF — it exists so that future
sessions do not re-read the same papers.

Right now the folder contains exactly one paper (in two files: main text and
supplement): the eSVD-DE methods paper itself. It is the *specification* for the
package, so the "key mathematical ideas" section below is written to be
sufficient for reasoning about the code without reopening the PDF — including
the mapping from each equation to the function that implements it.

The papers fall into 1 category:
- **The method itself**: the published description of what `eSVD2` is supposed to compute.

Papers are grouped by topic, not by date. Each entry is keyed by its
`generate-paper-id` citation key, the same key used in `brainstorming_[name].md`.
Superseded papers are marked `[SUPERSEDED by X]` rather than deleted.

This folder also holds working documents that are not paper summaries:

| File | What it is |
|---|---|
| `CRAN_READINESS.md` | The audit driving the CRAN submission. Its §0 is the current status (1.2.0, `R CMD check` result, open questions) |
| `UNIT_TEST_PLAN.md` | The test suite, with an ID and an oracle for each test. §2.18 to §2.21 cover the features of 1.1.0 and 1.2.0 |
| `OVERDISPERSION_BRAINSTORM.md` | Options for bounding the nuisance rate, and the decision taken (the cap of 1.2.0) |
| `TEST_RUN_REPORT.md`, `RMPFR_REPORT.md` | Dated records from 2026-08-29 to 2026-09-01: the suite's first run, and the `Rmpfr` experiment |
| `version_comparison/` | master (`3d5f7bf`) against the current `devel` on simulated cohorts, as a knitted report |
| `overdispersion_brainstorm/` | The dry-runs behind the brainstorm, and the report on the cap (`overdispersion_cap_claude.html`) |
| `lfc-se-comparison_2026-09-28_claude.R`, `rmpfr_experiment_claude.R` | One-off scripts |

Last updated: 2026-09-29 (the paper index itself is unchanged since 2026-08-27)

---

## The method itself

### lin2024esvd — Lin, Qiu & Roeder (2024)

**Files:** `s12859-024-05724-7.pdf` (main text, 30 pp.);
`12859_2024_5724_MOESM1_ESM.pdf` (Additional file 1 — supplement).

**Full citation:** Lin, K. Z., Qiu, Y., & Roeder, K. (2024). eSVD-DE: cohort-wide
differential expression in single-cell RNA-seq data using exponential-family
embeddings. *BMC Bioinformatics*, 25, 113. doi:10.1186/s12859-024-05724-7

**What it does:** Tests for differential expression **among individuals** (not
among cells) in cohort scRNA-seq data, by (i) fitting a Poisson matrix
factorization that pools information across genes to remove individual-level
confounders, and then (ii) testing on the *Gamma–Poisson posterior* rather than
on the fitted low-rank mean — which is what prevents the Type-1 error inflation
that ordinarily follows DE testing after dimension reduction.

**Key mathematical ideas** (with the implementing function in the package):

- **Hierarchical model** (Eq. 1–3). For gene `j`, cell `i`:
  `A_ji | λ_ji ~ Poisson(ℓ_ji · λ_ji)`, `λ_ji ~ Gamma(α = μ_ji/γ_j, β = 1/γ_j)`,
  so `E[λ_ji] = μ_ji` and `V[λ_ji] = γ_j·μ_ji` (constant Fano factor, as in SAVER).
  Low-rank mean `μ_ji = exp(Y_j·ᵀ X_i· + Z_j,(cc)·C_i,(cc))`;
  covariate-adjusted depth `ℓ_ji = exp(Z_j,-(cc)ᵀ C_i,-(cc))`.
  Marginally `A_ji ~ NB(r = μ_ji/γ_j, p = ℓ_ji/(ℓ_ji + 1/γ_j))`, giving
  `E[A_ji] = ℓ_ji·μ_ji` and `V[A_ji] = (ℓ_ji γ_j + 1)·ℓ_ji·μ_ji`.
  **⚠️ Parameterization gotcha:** the paper's overdispersion `γ_j` is the Gamma
  *scale*; the code's `nuisance_vec` is the Gamma **rate** `β = 1/γ_j` (hence the
  C++ function is named `gamma_rate`). Larger `nuisance_vec` = *less*
  overdispersion. Every comparison between paper equations and code must apply
  this inversion.
- **Covariate matrix `C`** (§ "Statistical model and method"). Required columns:
  intercept `C·,(int)`; log sequencing depth `C·,(lib) = log(Σ_j A_ji)`;
  case–control indicator `C·,(cc) ∈ {0,1}`; then one-hot-encoded categoricals
  (drop one level) and numerical covariates scaled to sd 1 but **deliberately not
  centered**, so that `Z` stays interpretable. Individual one-hot vectors are
  optional and often harmful (they make `C` collinear).
  → `format_covariates()`.
- **Initialization.** Per-gene Poisson ridge regression
  `Ẑ_j· = argmin_z −Σ_i [A_ji z ᵀC_i· − exp(zᵀC_i·)] + τ‖z_{−1}‖²₂` with τ ≈ 0.01
  (intercept unpenalized), then SVD of `R = log(A+1) − ẐᵀC = UDVᵀ`, setting
  `X̂ = V√D`, `Ŷ = U√D`. → `initialize_esvd()` (uses `glmnet::glmnet`).
- **Optimization** (Eq. 9–12). Minimize
  `−(1/np)Σ log P(A_ji | X,Y,C,Z) + τ(‖X‖²_F + ‖Y‖²_F + ‖Z‖²_F)` by alternating
  between `X` (convex given `Y,Z`) and `(Y,Z)` (convex given `X`). Each
  alternation decomposes into `n` k-dimensional and `p` (k+r)-dimensional
  problems, solved by **Newton's method** — chosen because it is deterministic
  (no stochasticity → reproducible fits across users) and because the gradients
  and Hessians of constrained exponential families naturally keep iterates
  feasible. Run in **two phases**: phase one holds `Ẑ·,(cc)` fixed at its
  initialized value, phase two frees it. Motivated by warm-starting results for
  non-convex problems. → `opt_esvd()`, `src/optimization.cpp`,
  `src/constrained_newton.cpp`.
- **Reparameterization** (identifiability, two steps).
  *Step 1:* regress each latent dimension `X̂·,d` on `C` (coefficients `β`,
  residual `ε`), update `Ẑ_j· ← Ẑ_j· + β·Ŷ_{j,d}` and `X̂·,d ← ε`. After all `d`,
  `X̂ᵀC = 0` and the fitted `ŶᵀX̂ + ẐᵀC` is unchanged.
  *Step 2:* SVD `R = ŶᵀX̂ = UDVᵀ`, set `X̂ ← (n/p)^{1/4}V√D`,
  `Ŷ ← (p/n)^{1/4}U√D`, making `X̂ᵀX̂/n` and `ŶᵀŶ/p` equal diagonal matrices.
  → `reparameterization_esvd_covariates()` (step 1) and `.reparameterize()` (step 2).
- **Overdispersion MLE.** `γ̂_j = argmax_γ Σ_i [A_ji log ℓ̂_i + μ̂_ji γ log γ
  − log Γ(γμ̂_ji) + log Γ(A_ji + γμ̂_ji) − (A_ji + γμ̂_ji)log(ℓ̂_i + γ)]`,
  solved by Newton's method. → `estimate_nuisance()`, `src/gamma_rate.cpp`
  (`gamma_rate` on the natural scale, `log_gamma_rate` on `ρ = log θ = −log β`
  for stability; the R wrapper tries the former and falls back to the latter).
- **Posterior** (Eq. 4, 14) — *the step that fixes the Type-1 inflation*.
  `λ_ji | A_ji ~ Gamma(α = A_ji + μ̂_ji/γ̂_j, β = ℓ̂_ji + 1/γ̂_j)`, so
  `μ̂^(post)_ji = (μ̂_ji/γ̂_j + A_ji)/(1/γ̂_j + ℓ̂_ji)` and
  `v̂^(post)_ji = (μ̂_ji/γ̂_j + A_ji)/(1/γ̂_j + ℓ̂_ji)²`.
  Testing on `μ̂_ji` directly instead would inflate Type-1 error, because the
  low-rank embedding makes genes highly correlated in `μ̂`.
  → `compute_posterior()`.
- **Test statistic** (Eq. 15). Average posteriors within individual `s`:
  `μ̂^(s)_j = mean_{i∈I(s)} μ̂^(post)_ji`, `v̂^(s)_j = mean_{i∈I(s)} v̂^(post)_ji`.
  Summarize the case group as one Gaussian:
  `μ̂^(case)_j = mean_{s∈A} μ̂^(s)_j` and
  `v̂^(case)_j = mean_{s∈A} v̂^(s)_j + mean_{s∈A}(μ̂^(s)_j)² − (mean_{s∈A} μ̂^(s)_j)²`
  (i.e. mixture-of-Gaussians variance = within + between). Then Welch:
  `T̂_j = (μ̂^(case)_j − μ̂^(control)_j)/√(v̂^(case)_j/|A| + v̂^(control)_j/|B|)`.
  → `compute_test_statistic()`, `.compute_mixture_gaussian_variance()`.
- **Multiple testing.** Welch–Satterthwaite df
  `df_j = (v^c/|A| + v^k/|B|)² / [(v^c/|A|)²/(|A|−1) + (v^k/|B|)²/(|B|−1)]`;
  transform `Ẑ_j = Φ⁻¹(F_{df_j}(T̂_j))`; estimate an **empirical null** with
  `locfdr::locfdr` (`fp0["mlest", c("delta","sigma")]`); recompute two-sided
  p-values against `N(delta, sigma²)`; BH-adjust; call FDR < 0.05.
  → `.compute_df()`, `compute_pvalue()`, `multtest()`.

**The hypothesis actually being tested** (supplement + § "High-level description").
With `μ^(s)_j ~ N(μ^(case)_j, v^(case)_j)` for each case individual `s ∈ A`, and
each individual contributing cells `λ_ji ~ F^(s)_j` with mean `μ^(s)_j`, eSVD-DE
tests `H_{0,j}: μ^(case)_j = μ^(control)_j`. This is *deliberately different from*
`H′_{0,j}: mean over all case cells = mean over all control cells` (which
pseudobulk targets but which ignores within-individual variance) and from
`H″_{0,j}` stated on `μ_ij` ignoring individual membership (which SCTransform-style
cell-level tests target — a handful of extreme case individuals can reject it).

**Temporal & methodological context:** 2024; human tissue (lung IPF — Adams
GSE136831 and Habermann GSE135893; colon ulcerative colitis — Smillie; brain
autism — Velmeshev); 10x Chromium droplet scRNA-seq; ~4,000–7,100 genes and
2,600–15,200 cells per analysis, 10–34 individuals, 4.4–39% sparsity. Baselines
compared: DESeq2 (pseudobulk), MAST (mixed effects), SCTransform (cell-level),
GLM-PCA (dimension reduction then Wilcoxon). Findings are tied to droplet
UMI-count data; the Gamma–Poisson assumption is the load-bearing one.

**Runtime context (supplement Section D):** eSVD-DE was slower than DESeq2/
SCTransform (minutes) but far faster than MAST (8.6 h and 84.8 h on the two
reported datasets); the authors explicitly flag speeding up the optimization
sub-routine as future work. Relevant when judging whether a proposed
CRAN-readiness change costs meaningful runtime.

**Project relevance (CRAN conversion):**
- **Specification of correctness**: this paper is the ground truth for
  `additional_context/CRAN_READINESS.md`'s correctness tests. Any proposed unit
  test that pins down "what the code should compute" should cite an equation here.
- **Three assumptions the package rests on** (stated in § Results): (1) counts are
  Gamma–Poisson, (2) covariate effects are removable in a GLM framework, (3) DE
  genes differ in *mean* between case and control individuals. Documentation
  should say these out loud; a CRAN reviewer reading `?opt_esvd` currently cannot
  find them.
- **Reproducibility claim to protect**: "there is no need to consider stochastic
  optimization schemes… different practitioners using our method would necessarily
  obtain the same resulting fit." Any change touching `src/` must preserve
  bitwise-deterministic output, and a regression test should enforce it.

**Key quotes or claims:**
- "a naive application of DE testing on the dimension-reduced scRNA-seq data has
  been observed to inflate the Type-1 error… because the dimension reduction
  introduces correlations among genes that contaminate the signal."
- "we do not perform a differential expression test based on μ_ji's because
  different genes are highly correlated based on their values in μ_ji due to the
  low dimensional embedding, which will artificially inflate the Type-1 error."
- "Since (10) and (11) decompose into n and p smaller optimization problems… there
  is no need to consider stochastic optimization schemes. This is appealing as
  this means different practitioners using our method would necessarily obtain the
  same resulting fit."
