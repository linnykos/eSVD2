# Memo: a second autoresearch loop, on rules that use the eSVD-DE embedding

Kevin Z. Lin (drafted by Claude), 2026-11-11. This memo fixes the design of
`autoresearch_lfcse2/`. It follows the first loop's memo
(`autoresearch_lfcse/notes/MEMO_claude.md`), and only what differs from it
is argued here.

**How this memo was agreed.** Kevin set the goal, the constraint and the
bank's purpose in one message on 2026-11-11 and asked for the loop to be
built and run overnight without further rounds. The choices below that his
message does not settle were made by Claude and are listed in section 6 for
him to confirm or reverse.

## 0. What the first loop left open

The first run (`autoresearch_lfcse/notes/RUN_nov11_claude.md`) kept
`donor_wls`, a donor-level weighted regression that reads the counts, the
design and the fitted library. It uses nothing of the embedding
$\widehat x_c^\top \widehat y_g$, the fitted case-control coefficient, the
nuisance rate or the posterior. Two facts limit that result. The first bank
drew every cell's state independently of its donor and its arm, so the
embedding had no donor-level structure to recover. Its donor effects were
independent across genes, which is the one structure a low-rank term cannot
represent. Kevin's instruction for this loop:

> I don't want to consider a method to compute the LFC and SE that relies
> only on $\widehat z_{gk}$ for the covariates $k$ that isn't the
> case-control arm, since this defeats the purpose of fitting eSVD. [...]
> My intent is to explore LFC and SE rules that better use the actual eSVD
> fit.

## 1. The goal and the score

**The object searched for** is unchanged: a rule that reports, from one
fitted cohort, a log2 fold change that is unbiased for the true one at
every effect size, and a standard error that is the sampling SD of that log
fold change.

**The constraint is new.** The log fold change or the SE must be computed
from the fitted embedding, explicitly (from `x_mat` and `y_mat`) or
implicitly (from the posterior mean or variance, which contain
$\exp(\widehat x_c^\top \widehat y_g + \widehat z_{g,\text{cc}} C_{c,\text{cc}})$).
A rule that reads only the counts, the design and the fitted library is
out of bounds. So is a rule in which the embedding enters as a token term
whose removal leaves the method unchanged; the run report states, for the
kept rule, how its output changes when the embedding is scrambled across
cells.

**The score** is unchanged: the Gaussian negative log-likelihood of the
truth under the reported pair, in nats per gene, lower is better,

$$\text{NLL} = \frac{1}{p}\sum_{g=1}^{p}\left[\tfrac12 \log\!\left(2\pi s_g^2\right) + \frac{(\widehat\beta_g - \beta_g)^2}{2 s_g^2}\right].$$

One example of this bank is one cohort (section 3), and its score is that
cohort's NLL.

**Headline deliverable of a run**: the kept `report_rule_claude.R`, the
self-contained rule file, `RUN_<tag>_claude.md`, and an `.Rmd` writeup.

## 2. The file on trial

`autoresearch_lfcse2/report_rule_claude.R` defines `report_lfc_se(view)`
under the first loop's contract: the same `view`, the same returned data
frame, the same allowed calls, the same sandbox, 30 seconds per cohort,
every cohort treated by one formula, every numeric constant with its
reason, and no shrinkage of the log fold change toward zero. The embedding
constraint of section 1 is added.

**The baseline is `null_rescaled`**: the pair that 1.2.0 stores, computed
from donor averages of the posterior mean and variance, shifted and
rescaled by the empirical null that 1.2.0 fitted. It is the best rule of
the first bank's reference rows that uses the posterior. The first loop's
baseline, `raw_null`, reads only the counts and the fitted library and is
out of bounds here.

**Two references are stored with the truth** and printed in every block:
`deseq2` (DESeq2 on donor pseudobulk counts) and `donor_wls` (the first
loop's kept rule). Neither uses the embedding. They show what the
embedding has to beat, and `donor_wls` is not a candidate.

## 3. The bank `lfcse-bank-2`

**One example** is one cohort of 600 genes, generated from the eSVD model
and fitted end to end by eSVD2 1.2.0. The first bank's examples were pairs
of cohorts, kept for two moments of Was2CoDE's Figure 4 that were never in
the score and were dominated by a few near-zero denominators; they are
dropped, which halves the fits per example.

**The generator.** For cell $c$ of donor $d$ in arm $a \in \{0, 1\}$ and
gene $g$,

$$\log \mu_{cg} = \alpha_g + x_c^\top y_g + a\,\beta_g + \delta_{dg} + (\text{sex and age terms}),
\qquad x_c = u_d + a\,s + e_c,$$

with $\lambda_{cg} \sim \text{Gamma}(\text{mean } \mu_{cg}, \text{rate } \gamma_g)$
and $A_{cg} \sim \text{Poisson}(\ell_c \lambda_{cg})$. Here $e_c \sim N(0, 0.4^2 I_3)$
is the cell's own state as in the first bank, $u_d$ is the state of its
donor, $s$ is a shift of the state of every case cell, and $\delta_{dg}$
is what remains of the donor's effect on gene $g$. The loadings $y_g$ are
centred over the genes. Each donor term is centred so that its exponential
has mean 1 within an arm, and the truth

$$\beta^{\text{true}}_g = \frac{\beta_g + s^\top y_g}{\ln 2}$$

is the log2 ratio of mean expression between the populations of case and
control donors at equal depth and equal covariates. This is the first
bank's definition, and it is the quantity that eSVD-DE's own
reparameterization assigns to the case-control coefficient, since the
fitted embedding is made orthogonal to the covariates.

**Factors.** Nine factors; the first six are the first bank's, with one
replaced.

| Factor | Levels | What it answers |
|---|---|---|
| `donors` | 10, 20, 40 | as before |
| `cells` per donor | 50 equal; unequal, 30 to 600 | as before |
| `effects` | `sparse`, `sparse_up`, `dense` | as before |
| `covariates` | none; donor sex and age | as before |
| `donor_sd` | 0.2; 0.5, times a log-normal factor per gene of log SD 0.3 | as before, and the donor SD now varies between genes |
| `expression` (*replaces `rates`*) | `high`: $\alpha_g$ from 0 to 1.5; `low`: from -3.5 to 0 | **capped genes.** The first bank capped 13% of genes against 73% in the real cohorts. The cap is 10 times a gene's median fitted library, so it binds for sparse genes; at `low` about half of the genes are capped or at the boundary. The rate is log-uniform on 0.5 to 200 at both levels |
| `state` (*new*) | `none`: $u_d = 0$, $s = 0$; `donor`: $u_d \sim N(0, \sigma_u^2 I_3)$ carrying about 72% of the donor variance on average, $\delta_{dg}$ the rest (the total donor variance is about 12% lower than at `none`); `donor_arm`: the same with $s \sim N(0, 0.25^2 I_3)$ | **the benefit of the low-rank term.** Donor-to-donor variation that is shared across genes is what an embedding can estimate, and removing it can lower the error of the log fold change. At `donor_arm` the arms differ in cell state and every gene has a nonzero effect |
| `arm_sd` (*new*) | `equal`; `unequal`: the donor variance of the cases is 1.5 times, and of the controls 0.5 times, the `equal` value | **the estimand.** A contrast of mean log expression differs from the log ratio of means by half the difference of the arms' log-scale variances, gene by gene |
| `depth` (*new*) | `none`; `error`: each donor shifts a block of 15% of the genes (one unit higher in $\alpha_g$) by $N(0, 0.5^2)$ | **donor-level depth error.** The shift moves the donor's total count, so the depth adjustment moves every other gene of that donor by a common offset |

**The design.** The full factorial has 1728 cells. Each split holds a
balanced 1/6 fraction of 288: a third of the combinations of the three
3-level factors crossed with a half of the combinations of the six 2-level
factors, chosen so that every level, and every pair of levels of two
factors, appears equally often. Dev and test use two disjoint fractions
and share no cell, so the test split checks a rule on combinations of
levels it was never scored on.

**The rank.** Every cohort is fitted at rank 5. The generator has three
state factors and, at `depth = error`, a fourth donor-level direction; a
user's rank is above the unknown truth. The first bank fitted at the true
rank 3.

**Truth stored per example**: the true log2 fold changes, the gene class,
the population, each donor's state and whole effect on each gene (for the
oracle rows of `BANK.md`), and the two references.

**A pilot of 16 fitted cohorts** (before the build) showed that the fitted
embedding tracks the donor-level states: the donor mean of
$\exp(\widehat x_c^\top \widehat y_g)$, centred within arm, has correlation
0.72 to 0.91 with the true low-rank donor effect, with slope 0.83 to 1.02.
An oracle that subtracts the true low-rank effect before `donor_wls` scores
-0.21 where `donor_wls` scores 0.52 (`state = donor`), so the bank holds a
gain of about 0.7 nats for a rule that uses the embedding well. The naive
rule that subtracts the fitted embedding offset scores far worse than
`donor_wls`, so the gain is not free.

## 4. Directions of variation

**E1. Donor summaries of the posterior.** The posterior mean of a cell is
$w\,A/\ell + (1 - w)\,\widehat\mu$ with $w = \ell / (\ell + \gamma)$. Its
donor averages pull donors toward the fit, so their spread understates the
noise by about $w$; a between-donor SE divided by the shrinkage weight is
the exploration's calibrated form. Fails if the fitted mean carries its
own error that the between-donor spread does not show.

**E2. The embedding as an offset.** Subtract from each donor's log
expression the log of its mean $\exp(\widehat x_c^\top \widehat y_g)$, in
full or by a coefficient estimated per gene or per cohort. Fails if the
embedding has absorbed part of the case-control effect or of the noise,
which the pilot's naive version suggests.

**E3. The embedding as donor-level covariates.** The donor means of
$\widehat x_c$ enter the donor-level design, as surrogate variables do.
Costs degrees of freedom: 5 columns against 10 donors.

**E4. The variance after the adjustment.** The residual variance, its
degrees of freedom, the uncertainty of the estimated embedding, and the
moderation across genes, once E2 or E3 has removed part of the
donor-to-donor variation.

**E5. The arm and the centring.** How the arm's share of the embedding is
treated (the reparameterization puts it in the case-control coefficient),
whether the fitted case-control coefficient helps once donors are weighted
properly, and the empirical null's shift.

**E6. Cell-level use.** The posterior mean as a denoised cell value, the
posterior variance as a within-donor term, the nuisance rate in the
weights.

**The queue of first experiments**, in order:

1. [E1] The between-donor SE of posterior donor means divided by the
   shrinkage weight, with the empirical null's shift (`deshrunk`).
2. [E1] A donor-level regression on the log of posterior donor means, with
   moderated variance.
3. [E2] The full embedding offset in the donor-level regression.
4. [E2] The embedding offset with one coefficient per gene, estimated in
   the regression.
5. [E3] Donor means of the embedding scores as columns of the design.
6. [E4] The moderated variance with degrees of freedom reduced for the
   estimated embedding.
7. [E5] The offset centred within arm, so that the arm contrast is
   unchanged.
8. [E6] The posterior variance as a within-donor term.

**Out of bounds this round**: a rule without the embedding (section 1);
shrinking the log fold change toward zero; returning or changing the
p-value; refitting the embedding inside the rule; any branch on the bank's
factor levels.

## 5. The protocol

The first loop's, with one addition.

| | |
|---|---|
| Keep rule | paired gain over the best commit positive and at least one paired SE over the 288 dev examples; a tie is kept when the file got shorter or at least twice as fast |
| Gate 1 | `slope_weak` in [0.9, 1.1] |
| Gate 2 (*new*) | `se_strong` in the interval `BANK.md` records. It is the pooled ratio of the SE of strong-effect genes (direct case-control effect $|\beta_g| / \ln 2 \ge 0.5$, which is independent of the loadings) to the SE that a smoothing spline, fitted to the other genes against the log mean count, gives at the same expression level. The first run withdrew by hand two commits whose SE followed the estimated log fold change; this gate lets the keep rule refuse them |
| Cores | 8 for the build and for a scoring |
| Budget | 40 experiments or 4 hours; stall after 12 consecutive non-keeps past 15 experiments |
| Scorer budget | 600 seconds |
| Bank home | `ESVD2_BANK_DIR`, default `tools::R_user_dir("eSVD2", "cache")`, folder `lfcse-bank-2/` |

## 6. Decisions made without a round, for Kevin to confirm

| # | Decision | Reason | What reversing it costs |
|---|---|---|---|
| 1 | The truth is the marginal effect $\beta_g + s^\top y_g$, cell-state differences between the arms included | it is the first bank's definition, it is identifiable (an effect at fixed cell state is not separable from a shift of the states), and it is where eSVD-DE's reparameterization puts the arm's share of the embedding | a truth at fixed cell state is a new bank |
| 2 | The benefit of the embedding is built in as donor-level states (variance reduction), with an arm shift at one level | an arm shift alone changes the truth and gives the embedding nothing to remove under decision 1 | new `state` levels are a new bank |
| 3 | Rank 5 for every fit | above | a rebuild |
| 4 | One cohort per example; the Figure 4 moments are dropped | above | pairs need twice the fits |
| 5 | `expression` replaces `rates` to raise the capped share | the cap binds through sparsity, as in the real cohorts | a rebuild |
| 6 | The second gate, with its interval set from the reference rules at the build | the first run's recommendation | editing one constant of the scorer and re-scoring |
| 7 | The baseline is `null_rescaled` | the best embedding-using reference of the first bank | re-baselining |
| 8 | The embedding constraint is enforced by the contract and by a scrambling check in the run report, not by the scorer | a rule that reads `stored_df` uses the posterior without touching `x_mat`, so a mechanical check inside the scorer would refuse the baseline | a scorer change |
| 9 | Dev and test use two disjoint fractions of the design | a stronger held-out check | none; the manifests record the cells |
| 10 | (Added during the session, after the first experiments of the first run.) The scorer prints the score against two truths and `ESVD2_TRUTH` selects which one is the score: `marginal` (decision 1, the default) or `direct`, the effect at fixed cell state $\beta_g / \ln 2$. A second run, `nov11-direct`, was scored against the direct effect | under the marginal truth the embedding cannot change the log fold change: the fitted scores are orthogonal to the arm, and the part of the arm contrast along the loadings mixes a real shift with donor sampling. The direct effect is the target for which removing that part is unbiased | none; both runs and both rule files are kept, and the choice of target is Kevin's |
