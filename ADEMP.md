# ADEMP: Reproducing Storey (2002) with a signal-strength extension

**Status:** The primary design was specified before the full run. Results are reported separately in ANALYSIS.md.

## Aims

The primary aim is to reproduce the numerical experiment in Section 5 of Storey (2002), using Table 2 and Figure 1 as reference outputs. We will assess agreement in power, false discovery rate (FDR), mean estimated FDR, mean estimated null proportion, and the BH rejection threshold.

The secondary aim is to examine how signal strength affects these quantities while retaining the original comparison protocol. We will report achieved FDR alongside power rather than interpreting higher power alone as superior performance. The study does not aim to establish a general finite-sample FDR guarantee for an adaptive procedure.

## Data-generating mechanisms

For each replicate, generate a family of $m=1000$ independent test statistics. Let $H_i=0$ denote a true null and $H_i=1$ an alternative. The primary implementation fixes $m_0=m*\pi_0$ and $m_1=m-m_0$. All selected values yield integer counts. Null statistics follow $N(0,1)$; alternative statistics follow $N(\mu,1)$. Compute one-sided p-values with the standard normal survival function.

Section 5's description does not completely resolve whether null counts were fixed in the original code. The theoretical model and errata specify independent Bernoulli labels. We adopt the fixed-count interpretation for the numerical reproduction, record this choice, and will report any additional Bernoulli-label sensitivity analysis separately. Analytic benchmarks will match the chosen label mechanism.

| Experiment | $m$ | $\pi_0$ | $\gamma$ | $\lambda$ | $\mu$ | Replicates per condition | Conditions |
|---|---:|---|---|---:|---|---:|---:|
| Reproduction | 1000 | 0.1,0.2,...,0.9 | 0.001,0.01 | 0.5 | 2 | 1000 | 18 |
| Extension | 1000 | 0.1,0.2,...,0.9 | 0.001,0.01 | 0.5 | 1,3 | 1000 | 36 |

The full factorial matrix is exported as ``results/full/raw/design_matrix.csv``. Both gamma values and both methods use the same generated family within each ($\pi_0$, $\mu$, replicate) combination. The 54 conditions therefore require 27,000 generated families and produce 108,000 method-level rows.

The master seed is 20261006. Stable integer identifiers and NumPy SeedSequence will define independent family-level streams. Neither the method nor gamma changes the DGP stream. Configuration, dependency versions, seed construction, elapsed time and available Git revision will be recorded.

## Estimands and targets

For a rejection rule, let $R$ be the total number of rejections, $V$ the number of rejected true nulls, and $S$ the number of rejected alternatives. FDR is $E[V/\max(R,1)]$. Power is $E[S/m_1]$ for the primary conditions, which all have $m_1>0$.

We also evaluate the estimation of the null fraction $\pi_0$ and the FDR of the fixed-threshold rule. The fixed-count analytic FDR target is computed from independent $V\sim\text{Binomial}(m_0, \gamma)$ and $S\sim\text{Binomial}(m_1,G_{\mu}(\gamma))$, where $G_{\mu}(\gamma)=\Phi(\mu-\Phi^{-1}(1-\gamma))$. It can be evaluated by binomial enumeration for small families and by one-dimensional quadrature for the full design.

Optional pFDR summaries condition on $R>0$ and will report the number of eligible replicates. If none are eligible, pFDR is undefined. This is distinct from the FDR estimate used to set the BH comparison level.

## Methods

### Fixed-threshold rule with Storey's FDR estimate

Reject when $p_i\leq\gamma$. Let $W_{\lambda}$ count $p_i>\lambda$ and $R_{\gamma}$ count $p_i\leq\gamma$. Compute:

```text
pi0_hat_raw = W_lambda / [m*(1-lambda)]
fdr_hat_raw = pi0_hat_raw*m*gamma / max(R_gamma,1)
fdr_hat = min(1,fdr_hat_raw)
```

The main implementation preserves the untruncated null-proportion formula from Section 3 and caps the FDR estimate at 1 for its operational use. Both raw and capped estimates will be retained. A variant that caps ``pi0_hat`` at 1, as discussed in the errata's comparison, would be identified separately rather than substituted silently.

### BH with a replicate-specific comparison level

Apply BH to the same p-values with $\alpha_{rep}=\widehat{FDR}$. For ordered p-values, select the largest k satisfying $p_{(k)}\leq k*\alpha_{rep}/m$; select none if no such k exists. Record the rejection threshold $p_{(k)}$, or 0 if no hypotheses are rejected.

This reproduces the comparison protocol rather than comparing methods at a prespecified $\alpha=0.05$. Fixed-alpha BH guarantees cannot automatically be applied to this data-dependent level. Methods receive p-values and method parameters, but not truth labels, true pi0 or the alternative mean.

## Performance measures

At replicate level, save $R$,$V$,$S$,$m_0$,$m_1$,FDP,power, the raw null-proportion estimate, the raw and capped FDR estimates, the BH input level and its rejection threshold.

For each condition and method report:

- Mean FDP (achieved FDR), mean power and rejection probability.
- Monte Carlo standard errors (MCSEs) and approximate 95% Monte Carlo intervals where appropriate.
- Mean, bias and RMSE of ``pi0_hat_raw``.
- Mean, bias and RMSE of the fixed-threshold FDR estimate relative to its matching analytic target; identify raw and capped versions.
- The paired power difference between methods, with its paired MCSE.
- Comparison to the relevant rounded entries of Table 2, acknowledging unavailable original seeds and original Monte Carlo uncertainty.

For a replicate statistic q, MCSE is ``sd(q)/sqrt(K)``. For paired differences, compute the SD of the replicate differences. Original paper discrepancies will not be treated as exact equality failures or as formal significance tests using only our MCSE. The FDR is the mean of replicate FDPs, not the pooled ratio sum(V)/sum(R).

Generic metric functions return missing power when $m_1=0$. All zero-rejection replicates remain in FDR summaries.

## Visualizations

1. A DGP diagnostic displaying null and alternative p-value distributions, with a uniform null reference.
2. A power comparison across pi0, faceted by gamma and mu, with Monte Carlo uncertainty and an explicit explanation of the BH level.
3. An FDR calibration plot comparing achieved fixed-threshold FDR, the corresponding analytic target, mean estimated FDR and achieved BH FDR.

Baseline and extension panels will be distinguished. Figures will include m,K,lambda and label-mode information, use a consistent facet order, and be exported as PDF and PNG.

## Validation and limitations

Tests will cover hand-computed BH decisions, index alignment, zero-rejection metrics, zero-alternative handling, FDR formula edge cases, deterministic seeds and complete small-run outputs. Independent analytic checks will validate fixed-threshold power and FDR calculations.

The study assumes independent tests, a single positive alternative mean within each condition, normal test statistics and a fixed lambda. It does not examine arbitrary dependence, heterogeneous effects, modern q-value implementations or bootstrap tuning. The fixed-count label interpretation and estimate truncation choices are explicitly recorded. Deviations from this prospective design will be logged with reasons and separated from planned analyses.

## Post-run diagnostic addition (2026-10-16)

After the primary run identified a BH power discrepancy at pi0=0.8, gamma=0.001, mu=2, two additional 5,000-replicate configurations were added: check-row-fixed.json and check-row-bernoulli.json. Both use the new seed 20261016 and retain the original estimation and comparison rules. These post hoc checks examine whether the label mechanism explains the mismatch; they do not replace or enlarge the prespecified primary 54-condition analysis. The discrepancy and both checks are reported in ANALYSIS.md.

## References

- Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B, 64(3), 479–498. https://genomics.princeton.edu/storeylab/papers/directfdr.pdf
- Storey, J. D. Errata dated April 16, 2003. https://genomics.princeton.edu/storeylab/papers/directfdr_corrections.pdf
