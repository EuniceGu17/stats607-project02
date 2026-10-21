# Reproducing Storey (2002): FDR estimation and power

My repository URL is [here](https://github.com/EuniceGu17/stats607-project02).

This project reproduces the simulation comparison in Section 5 of Storey
(2002), using Table 2 and Figure 1 as references, and extends it to weaker and
stronger signals. The Python pipeline generates independent tests, applies
both methods to shared data, saves replicate-level results, and produces
performance summaries and figures.

**Key findings:** The baseline reproduces fixed-threshold power near 0.138
and 0.372, and 17 of 18 BH power entries differ from the published values by
less than 0.003. One BH entry remains unmatched; the signal-strength
extension also shows that higher power can accompany substantially higher
achieved FDR.

## Setup and complete run

Use Python 3.11+ and `make`. Run commands from the repository root. All data
are simulated; no external dataset is required.

```bash
git clone https://github.com/EuniceGu17/stats607-project02.git
cd stats607-project02
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make test PYTHON=python
make all PYTHON=python
```

The default configuration is `configs/full.json`. For a quick execution
check before the complete study:

```bash
make all CONFIG=configs/smoke.json PYTHON=python
```

The tests and full pipeline were independently verified using the repository
code: **26 tests passed**, the rebuilt raw file matched the SHA-256 recorded
in committed metadata, and full summaries matched within numerical precision.
On macOS 15.7.4 ARM64 with Python
3.11.9, a verification run took approximately **10.1 seconds for simulation,
0.5 seconds for analysis, and 1.4 seconds for figures**, excluding dependency
installation and process startup. Allow roughly a minute for a first run on
similar hardware; runtime varies. The verification environment used NumPy
2.4.6, SciPy 1.17.1, pandas 3.0.6, Matplotlib 3.11.2, and pytest 9.1.1.
The repository's `requirements.txt` specifies dependency ranges rather than
a frozen environment.

## Study design

The baseline uses 1,000 tests per family, alternative mean `mu=2`,
`pi0=0.1,...,0.9`, p-value thresholds `gamma=0.001,0.01`, and `lambda=0.5`.
The extension adds `mu=1,3`. Each of the 54 conditions uses 1,000 replicates.
There are 27,000 generated families because both thresholds and both methods
share data within each `(pi0, mu, replicate)` combination, producing 108,000
method-level raw rows and 108 method-condition summaries.

The primary DGP fixes the null count at `m*pi0`. Null statistics follow
`N(0,1)` and alternative statistics follow `N(mu,1)`. P-values are one-sided
normal survival probabilities. The master seed is `20261006`; family-level
streams depend on the seed, DGP parameters and replicate index, independently
of the method or threshold. Truth labels are used only for evaluation.

The fixed rule rejects `p<=gamma` and estimates that rule's FDR. BH is applied
to the same p-values using the estimated FDR as its replicate-specific input
level. This follows the paper's comparison protocol; it is not a comparison
at a prespecified `alpha=0.05`. Fixed-alpha BH guarantees do not automatically
apply to this data-dependent level. The null-proportion estimate is not
capped; the operational FDR estimate is capped at 1, and both raw and capped
values are retained.

See [ADEMP.md](ADEMP.md) for the complete design matrix, estimands and
implementation choices, and [ANALYSIS.md](ANALYSIS.md) for interpretation.

## Commands and configurations

| Command | Action |
|---|---|
| `make all` | Run simulation, analysis and figure generation in order |
| `make simulate` | Generate raw results or verify and reuse their cache |
| `make analyze` | Ensure raw results, validate rows and calculate summaries |
| `make figures` | Ensure upstream results and generate all figures |
| `make test` | Run unit tests and a small end-to-end test |
| `make clean` | Remove generated raw, summary and figure directories for the selected output profile |

Add `PYTHON=python` when using the activated virtual environment.

| Configuration | Conditions | Replicates per condition | Method-level rows |
|---|---:|---:|---:|
| `configs/smoke.json` | 4 | 20 | 160 |
| `configs/pilot.json` | 54 | 100 | 10,800 |
| `configs/full.json` | 54 | 1,000 | 108,000 |

Outputs default to `results/<configuration-name>`. For a separate run:

```bash
make all CONFIG=configs/full.json OUT=results/new-run PYTHON=python
```

Valid caches require matching configuration, simulation-source hashes,
dependency versions and raw-file hashes. To force a new full simulation:

```bash
python -m src.simulation --config configs/full.json --force
make figures CONFIG=configs/full.json PYTHON=python
```

## Repository structure and outputs

```text
configs/          # smoke, pilot, full, and two discrepancy-check designs
src/
  dgps.py         # data generation and reproducible random streams
  methods.py      # fixed rule, FDR estimation, and BH
  metrics.py      # replicate metrics and analytic checks
  simulation.py   # simulation driver and raw output
  analyze.py      # validation, summaries, and paper comparison
  figures.py      # diagnostic and result plots
  common.py       # configuration, file operations, and cache checks
  clean.py        # restricted cleanup of generated outputs
tests/            # method, DGP, metric, and pipeline tests
examples/         # one-family and one-condition examples
references/       # rounded Table 2 values and provenance
data/simulated/   # optional simulated-data cache
results/full/
  summaries/
    summary.csv
    paired_power.csv
    comparison_table2.csv
    validation.json
    metadata.json
  figures/        # four figures, each as PNG and PDF
```

`make simulate` also creates `results/full/raw/replicates.csv`,
`design_matrix.csv`, `diagnostic.csv` and `metadata.json`. Raw directories are
Git-ignored and are regenerated after a fresh clone. Committed summaries and
figures can be inspected without running the simulation.

In raw output, `R,V,S` are rejection counts, `FDP` and `power` are
replicate-level performance, `fdr_hat` is the estimated fixed-rule FDR, and
`alpha_rep` is the BH input. In summaries, `fdr` is achieved FDR, calculated
as the mean FDP including zero-rejection families. It is not the pooled
ratio of false to total discoveries. `true_fdr` and `true_power` always refer
to the fixed-rule analytic targets, including when repeated on a BH row;
they are not BH analytic values. Likewise, mean `fdr_hat` on a BH row is its
mean input level, not an estimate of achieved BH FDR.

## Figures

| Figure | Purpose |
|---|---|
| [DGP diagnostic](results/full/figures/01_dgp_diagnostic.png) | Null and alternative p-value distributions from preselected families |
| [Baseline power](results/full/figures/02_baseline_power.png) | Reproduction with published BH points and analytic fixed-rule power |
| [Signal-strength extension](results/full/figures/03_signal_extension.png) | Power for `mu=1,2,3` |
| [FDR calibration](results/full/figures/04_fdr_calibration.png) | Achieved FDR, mean estimated fixed-rule FDR, and analytic references |

Result plots use approximate 95% Monte Carlo intervals. These quantify the
precision of simulation averages; they are neither simultaneous bands nor
guarantees for the FDP of an individual family.

## Unmatched paper entry and additional checks

At `pi0=0.8, gamma=0.001, mu=2`, the primary BH power is **0.112985**
(MCSE 0.000831), versus the published **0.129**. Two additional
5,000-replicate checks use seed `20261016`:

```bash
make analyze CONFIG=configs/check-row-fixed.json PYTHON=python
make analyze CONFIG=configs/check-row-bernoulli.json PYTHON=python
```

Their committed summaries report BH power **0.114279** (MCSE 0.000387) with
fixed counts and **0.113573** (MCSE 0.000392) with Bernoulli labels. These are
post hoc diagnostics, separate from the primary design. The published value
is retained; this unresolved difference does not establish a paper error.

## References

Storey, J. D. (2002). A direct approach to false discovery rates.
*Journal of the Royal Statistical Society: Series B*, 64(3), 479–498.
[Paper](https://genomics.princeton.edu/storeylab/papers/directfdr.pdf).
[Author's errata](https://genomics.princeton.edu/storeylab/papers/directfdr_corrections.pdf).
