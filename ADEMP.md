# ADEMP: Storey (2002) reproduction and signal-strength extension

This document records the implemented primary study. The baseline reproduces
Section 5, Table 2 and Figure 1 of Storey (2002); the extension varies signal
strength. The additional discrepancy checks are identified separately below.
Results and interpretation are reported in [ANALYSIS.md](ANALYSIS.md).

## Aims

1. Assess how closely the implemented experiment reproduces the paper's
   fixed-threshold FDR, both methods' power, mean estimated FDR, mean estimated
   null proportion, and mean BH rejection threshold.
2. Evaluate the accuracy and variability of fixed-rule FDR estimation across
   null proportions, reporting achieved FDR alongside power.
3. Examine whether the baseline conclusions persist with weaker and stronger
   alternatives while retaining the same comparison protocol.

The study evaluates a specified numerical comparison. It does not attempt
to prove general finite-sample FDR control for an adaptive procedure.

## Data-generating mechanisms

Each replicate is a family of $m=1000$ independent tests. Define $H_i=0$ for
a true null and $H_i=1$ for an alternative. The primary study fixes
$m_0=m\pi_0$ nulls and $m_1=m-m_0$ alternatives; all selected counts are
integers. Independently generate

$$
Z_i\mid H_i=0\sim N(0,1),\qquad
Z_i\mid H_i=1\sim N(\mu,1),\qquad
p_i=1-\Phi(Z_i).
$$

The implementation uses `scipy.stats.norm.sf` for numerical stability.
Section 5's proportion-based description leaves the original label-generation
implementation incompletely specified. The theoretical model and errata use
Bernoulli labels. Fixed counts are an explicit numerical-reproduction choice;
analytic benchmarks match that choice. Random-label diagnostics do not
replace the primary design.

### Factors and replication

| Component | Specification |
|---|---|
| Tests per family | $m=1000$ |
| Null proportions | $\pi_0=0.1,0.2,\ldots,0.9$ |
| Rejection thresholds | $\gamma\in\{0.001,0.01\}$ |
| Estimation parameter | $\lambda=0.5$ |
| Baseline alternative mean | $\mu=2$ |
| Extension means | $\mu\in\{1,3\}$ |
| Replicates per condition | $K=1000$ |
| Main seed | `20261006` |
| Label mode and caps | `fixed`, `cap_pi0=false`, `cap_fdr=true` |

The complete primary matrix appears in the appendix. It has 18 baseline and
36 extension conditions. For each $(\pi_0,\mu,b)$ combination, one generated
family is shared across both thresholds and both methods. Thus, 27,000
families produce 54,000 condition-replicate evaluations and 108,000
method-level rows. Shared-family pairing is retained in comparisons.

`configs/full.json` defines the formal study. `pilot.json` uses the same
conditions with 100 replicates. `smoke.json` uses $\pi_0=0.1,0.9$, $\mu=2$,
both thresholds and 20 replicates; it checks execution, not statistical
conclusions.

### Random numbers and saved output

`family_rng` creates a NumPy `SeedSequence` from the master seed, $m$,
replicate index, and the binary representations of $\pi_0$ and $\mu$.
Neither the method nor $\gamma$ changes the family stream. The implemented
ordering and parameter representation are therefore part of reproducibility.

Raw output records family and scenario identifiers, replicate index,
parameters, truth counts, decisions' performance, and FDR estimates. Metadata
records configuration, versions, available Git revision, hashes and elapsed
time. Raw files are saved before analysis and regenerated after a fresh clone;
committed summaries and figures retain the reported results.

## Estimands and targets

For a method's rejection set, let $R$ be total rejections, $V$ false
rejections and $S$ true rejections. Primary performance targets are

$$
\mathrm{FDR}=E\!\left[\frac{V}{\max(R,1)}\right],\qquad
\mathrm{Power}=E\!\left[\frac{S}{m_1}\right].
$$

All primary conditions have $m_1>0$. Generic functions return missing power
when there are no alternatives. Secondary targets are $\pi_0$, the fixed-rule
FDR, rejection probability, and
$\mathrm{pFDR}=E[V/R\mid R>0]$. The latter is estimated only among eligible
replicates and is undefined if none have rejections.

For independent fixed-count tests, write
$G_\mu(\gamma)=\Phi(\mu-\Phi^{-1}(1-\gamma))$. Fixed-rule power is
$G_\mu(\gamma)$, and its FDR is calculated from independent
$V\sim\mathrm{Binomial}(m_0,\gamma)$ and
$S\sim\mathrm{Binomial}(m_1,G_\mu(\gamma))$. The implementation evaluates
the corresponding one-dimensional integral, with independent small-family
enumeration as a check. Bernoulli-label runs use their own mixture-model
benchmark. These targets apply to the fixed rule, not to BH.

## Methods

### Fixed rejection threshold and FDR estimation

Reject all $p_i\leq\gamma$. Let $R_\gamma=\#\{p_i\leq\gamma\}$ and
$W_\lambda=\#\{p_i>\lambda\}$. Calculate

$$
\widehat\pi_0=\frac{W_\lambda}{m(1-\lambda)},\qquad
\widehat{\mathrm{FDR}}_{\rm raw}(\gamma)
=\frac{\widehat\pi_0 m\gamma}{\max(R_\gamma,1)},\qquad
\widehat{\mathrm{FDR}}(\gamma)
=\min\{1,\widehat{\mathrm{FDR}}_{\rm raw}(\gamma)\}.
$$

The primary estimate of $\pi_0$ is untruncated and may exceed 1. The FDR
estimate is capped at 1 for operational use; raw and capped values are saved.
This cap and the choice not to cap $\widehat\pi_0$ are explicit implementation
decisions. A theorem for an untruncated estimate is not automatically a
guarantee for its truncated version or for each observed FDP.

### BH at the replicate-specific estimated level

For the same p-values, set
$\alpha_b=\widehat{\mathrm{FDR}}_b(\gamma)$. Sort them and choose the largest
$k$ such that $p_{(k)}\leq k\alpha_b/m$, rejecting the first $k$ hypotheses.
If no such index exists, reject none. Save the largest rejected p-value as
the BH threshold, or 0 when there are no rejections.

The level is estimated from the same data. This is the paper's numerical
comparison protocol, not a comparison at a common prespecified level such as
0.05. Ordinary fixed-alpha BH guarantees cannot automatically be applied.
Both methods receive p-values and method parameters only; truth labels,
true $\pi_0$ and $\mu$ are reserved for generation and evaluation.

## Performance measures and analysis

For every replicate-method pair, save $m_0,m_1,R,V,S$, FDP, power, $W_\lambda$,
raw and used null-proportion estimates, raw and capped FDR estimates, the BH
input level and its rejection threshold. Retain zero-rejection families.

For each condition and method, report achieved FDR as mean FDP, mean power,
rejection probability, pFDR and their Monte Carlo standard errors (MCSEs).
For the null-proportion estimate, report mean, bias and RMSE. For both raw
and capped fixed-rule FDR estimates, report mean, bias and RMSE against the
matching analytic FDR. No such estimation-error claim is made for achieved
BH FDR: the fixed-rule estimate is merely BH's input level.

For a replicate statistic $q_b$,
$\mathrm{MCSE}=\mathrm{SD}(q_b)/\sqrt{K}$, using sample SD. Paired power
differences are formed within replicate before calculating their mean and
MCSE. Approximate pointwise 95% Monte Carlo intervals are mean $\pm1.96$
MCSE, not simultaneous bands or individual-family FDP guarantees.
FDR is not estimated by the pooled ratio $\sum V/\sum R$.

`references/table2.csv` retains rounded published values. The comparison
output reports published value, simulation value, MCSE and difference for
six fields in each baseline condition. Original seeds and MCSEs are
unavailable, so differences are not treated as exact-equality failures or
formal two-study significance tests based only on the new MCSE.

## Visualizations and validation

The pipeline produces four PNG/PDF pairs: a preselected-replicate DGP
diagnostic, baseline power with published BH points, signal-strength power
panels, and FDR calibration. Result plots distinguish methods, order factors
consistently and display Monte Carlo uncertainty. FDR calibration includes
achieved BH FDR, which was absent from the original Table 2 comparison.

The 26 automated tests cover hand-computed and independently implemented BH,
edge cases, DGP properties, deterministic streams, analytic targets versus
enumeration, paired MCSE, cache integrity and a small end-to-end pipeline.
Raw validation additionally checks row counts, scenario-replicate keys,
parameter consistency, paired families, count identities and metric formulas.

## Additional discrepancy diagnostics

After identifying the unmatched BH power entry at
$\pi_0=0.8,\gamma=0.001,\mu=2$, two separate 5,000-replicate checks use seed
`20261016`. `check-row-fixed.json` retains fixed counts;
`check-row-bernoulli.json` independently assigns each test a null label with
probability 0.8. Other method settings remain unchanged. These checks are
post hoc, do not replace primary results, and are reported separately.

## Scope and limitations

The primary study assumes independent normal test statistics, a common
positive signal mean per condition and fixed $\lambda$. It does not assess
correlated tests, heterogeneous effects, bootstrap tuning or a general
q-value implementation. The original comparison's data-dependent level,
fixed-count interpretation and estimate truncation limit the conclusions.

## References

Storey, J. D. (2002). A direct approach to false discovery rates.
*JRSS B*, 64(3), 479–498.
[Paper](https://genomics.princeton.edu/storeylab/papers/directfdr.pdf);
[author's errata](https://genomics.princeton.edu/storeylab/papers/directfdr_corrections.pdf).

## Appendix: complete primary condition matrix

All rows use $m=1000$, $\lambda=0.5$, fixed null counts, and 1,000 replicates.
Each row evaluates both methods. This matrix is generated from the same
factor values as `configs/full.json`.

| ID | Analysis | mu | pi0 | gamma |
|---:|---|---:|---:|---:|
| 1 | Extension | 1 | 0.1 | 0.001 |
| 2 | Extension | 1 | 0.1 | 0.01 |
| 3 | Extension | 1 | 0.2 | 0.001 |
| 4 | Extension | 1 | 0.2 | 0.01 |
| 5 | Extension | 1 | 0.3 | 0.001 |
| 6 | Extension | 1 | 0.3 | 0.01 |
| 7 | Extension | 1 | 0.4 | 0.001 |
| 8 | Extension | 1 | 0.4 | 0.01 |
| 9 | Extension | 1 | 0.5 | 0.001 |
| 10 | Extension | 1 | 0.5 | 0.01 |
| 11 | Extension | 1 | 0.6 | 0.001 |
| 12 | Extension | 1 | 0.6 | 0.01 |
| 13 | Extension | 1 | 0.7 | 0.001 |
| 14 | Extension | 1 | 0.7 | 0.01 |
| 15 | Extension | 1 | 0.8 | 0.001 |
| 16 | Extension | 1 | 0.8 | 0.01 |
| 17 | Extension | 1 | 0.9 | 0.001 |
| 18 | Extension | 1 | 0.9 | 0.01 |
| 19 | Baseline | 2 | 0.1 | 0.001 |
| 20 | Baseline | 2 | 0.1 | 0.01 |
| 21 | Baseline | 2 | 0.2 | 0.001 |
| 22 | Baseline | 2 | 0.2 | 0.01 |
| 23 | Baseline | 2 | 0.3 | 0.001 |
| 24 | Baseline | 2 | 0.3 | 0.01 |
| 25 | Baseline | 2 | 0.4 | 0.001 |
| 26 | Baseline | 2 | 0.4 | 0.01 |
| 27 | Baseline | 2 | 0.5 | 0.001 |
| 28 | Baseline | 2 | 0.5 | 0.01 |
| 29 | Baseline | 2 | 0.6 | 0.001 |
| 30 | Baseline | 2 | 0.6 | 0.01 |
| 31 | Baseline | 2 | 0.7 | 0.001 |
| 32 | Baseline | 2 | 0.7 | 0.01 |
| 33 | Baseline | 2 | 0.8 | 0.001 |
| 34 | Baseline | 2 | 0.8 | 0.01 |
| 35 | Baseline | 2 | 0.9 | 0.001 |
| 36 | Baseline | 2 | 0.9 | 0.01 |
| 37 | Extension | 3 | 0.1 | 0.001 |
| 38 | Extension | 3 | 0.1 | 0.01 |
| 39 | Extension | 3 | 0.2 | 0.001 |
| 40 | Extension | 3 | 0.2 | 0.01 |
| 41 | Extension | 3 | 0.3 | 0.001 |
| 42 | Extension | 3 | 0.3 | 0.01 |
| 43 | Extension | 3 | 0.4 | 0.001 |
| 44 | Extension | 3 | 0.4 | 0.01 |
| 45 | Extension | 3 | 0.5 | 0.001 |
| 46 | Extension | 3 | 0.5 | 0.01 |
| 47 | Extension | 3 | 0.6 | 0.001 |
| 48 | Extension | 3 | 0.6 | 0.01 |
| 49 | Extension | 3 | 0.7 | 0.001 |
| 50 | Extension | 3 | 0.7 | 0.01 |
| 51 | Extension | 3 | 0.8 | 0.001 |
| 52 | Extension | 3 | 0.8 | 0.01 |
| 53 | Extension | 3 | 0.9 | 0.001 |
| 54 | Extension | 3 | 0.9 | 0.01 |
