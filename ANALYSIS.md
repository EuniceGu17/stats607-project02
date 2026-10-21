# Critical analysis

The primary study broadly reproduced Storey's (2002) Section 5 experiment.
Across the baseline conditions, fixed-threshold power ranged from 0.1364
to 0.1384 at gamma=0.001 and from 0.3707 to 0.3728 at gamma=0.01, close to
the published values near 0.138 and 0.372. Seventeen of eighteen BH power
entries differed by less than 0.003. The exception was pi0=0.8,
gamma=0.001: BH power was 0.112985 (MCSE 0.000831), compared with 0.129 in
the paper. Two additional 5,000-replicate checks gave 0.114279 with fixed
null counts and 0.113573 with Bernoulli labels. The label mechanism therefore
does not appear to explain this discrepancy. The published value is retained;
without original code and seeds, the difference remains unresolved rather
than establishing a paper error. Numerical comparisons are available in
[the Table 2 comparison](results/full/summaries/comparison_table2.csv).

The comparison's interpretation requires care. BH uses the fixed rule's
estimated FDR as a replicate-specific input level, so the methods are not
compared at a common prespecified nominal level. Higher power alone cannot
establish better error control. Including BH provides a relevant predecessor,
but the experiment is not a broad comparison of adaptive FDR procedures.
Independent normal tests, one alternative mean and fixed lambda make the
design interpretable while limiting its realism. Examining both favorable
and unfavorable signal strengths is more informative than treating the
baseline as representative of all multiple-testing problems.

The extension showed why signal strength matters. At pi0=0.1, the mean null
proportion estimate was 0.3842 for mu=1, 0.1407 for mu=2 and 0.1023 for mu=3.
Weak alternatives contribute large p-values that inflate this estimate and
change the BH comparison level. A surprising reversal occurred at mu=1,
pi0=0.9, gamma=0.001: BH power was 0.2414 versus 0.0179 for the fixed rule,
but achieved FDR was also higher, 0.4914 versus 0.3310. This reverses the
baseline power ordering without demonstrating better error control.
[Committed summaries](results/full/summaries/summary.csv) support these values.

The figures recreate the baseline comparison while adding Monte Carlo
uncertainty, analytic fixed-rule references, signal-strength panels and
achieved BH FDR. The calibration plot makes the cost of the power reversal
visible, which a power-only display would obscure. The baseline plot also
retains the unmatched published point rather than hiding it. Approximate
pointwise intervals quantify simulation precision; they are not simultaneous
bands or guarantees about one family's false discovery proportion.

The main implementation challenge was distinguishing the rejection rule,
estimated FDR, achieved FDR and BH input level. Confidence comes from shared
families, truth labels restricted to evaluation, retention of zero-rejection
families, 26 passing tests, independent BH and binomial checks, and agreement
with analytic fixed-rule targets. Fixed versus random labels and truncation
remain explicit design choices. A stronger follow-up would vary dependence
and effect heterogeneity and separately compare procedures at prespecified
error levels, keeping those questions distinct from this reproduction.
