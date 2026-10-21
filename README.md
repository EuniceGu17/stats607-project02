# Reproducing Storey (2002): FDR estimation and power

This Python project reproduces the Section 5 simulation comparison in Storey
(2002) and extends it to alternative means mu=1 and mu=3. The main study uses
54 conditions, 1,000 replicates per condition, and two methods evaluated on
shared families of independent tests. It generates 108,000 method-level rows.

**Key findings:** Fixed-threshold power is approximately 0.138 at gamma=0.001
and 0.372 at gamma=0.01, closely reproducing the published trend. Most BH power
entries are close, but the published value 0.129 at pi0=0.8, gamma=0.001 is not
reproduced; separate fixed-count and Bernoulli-label checks give about 0.114.

For a step-by-step explanation in Chinese, read [IMPLEMENTATION_GUIDE_ZH.md](IMPLEMENTATION_GUIDE_ZH.md).
Start with `python -m examples.one_replicate` to follow one complete replicate.

## Setup

Python 3.11+ and `make` are recommended. No external data download is required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make test PYTHON=python
make all CONFIG=configs/smoke.json PYTHON=python
make all CONFIG=configs/full.json PYTHON=python
```

`requirements.txt` specifies supported dependency ranges. `requirements-lock.txt`
records the package versions used in the validated Python 3.11.9 environment;
install it instead to request those versions. The interpreter and dependency
versions used for every simulation are recorded in raw metadata.

With dependencies already installed for `python3`, simply run `make test` and
`make all`. A failed shell command stops the Make pipeline. Work from the
project root so Python can import the `src` package.

## Commands and workflow

| Command | Behavior |
|---|---|
| `make all` | Simulate, analyze, generate figures; full configuration by default |
| `make simulate` | Generate raw replicate-level results or verify/reuse the cache |
| `make analyze` | Ensure simulation results, validate and summarize them |
| `make figures` | Ensure upstream results and create plots |
| `make test` | Run unit tests and a small end-to-end pipeline test |
| `make clean` | Delete generated stages under results/full only |

Use `CONFIG=configs/pilot.json` for 100 replicates or
`CONFIG=configs/smoke.json` for a small execution check. The output folder
defaults to `results/<configuration-file-stem>`. Set `OUT=results/my-run` to
keep a separate run. Never direct concurrent writers to the same OUT.

To force a simulation rerun:

```bash
python -m src.simulation --config configs/full.json --output results/full --force
make figures
```

Raw reuse requires matching configuration, simulation-source hashes,
dependency versions and raw-file hashes. Increasing the number of replicates
rebuilds the file, but leaves the previously indexed random streams unchanged.
Analysis and plotting are rerun cheaply without rerunning valid simulations.

## Design and interpretation

Read [ADEMP.md](ADEMP.md) for the prospective design and implementation choices.
The baseline uses m=1000, mu=2, lambda=0.5, pi0=0.1,...,0.9 and
gamma in {0.001,0.01}. The extension adds mu=1 and mu=3. The main seed is
20261006. Each family has a fixed null count m*pi0; truth labels are available
only for evaluation, not to the methods.

The fixed rule rejects p<=gamma and estimates its FDR. BH is applied to the
same p-values at that replicate's estimated FDR level. This is **not** a
comparison at a prespecified alpha=0.05, and ordinary fixed-alpha BH guarantees
cannot automatically be applied. Raw null-proportion estimates are retained;
the FDR estimate used as a level is capped at 1.

FDR is the mean of replicate FDPs, including zero-rejection families. It is
not the ratio of pooled false discoveries to pooled discoveries. Intervals
are approximate Monte Carlo intervals, not simultaneous bands or guarantees
for the FDP in an individual family.

## Outputs

```text
results/full/
  raw/replicates.csv             # 108,000 rows before any summarization
  raw/design_matrix.csv         # all 54 condition combinations
  raw/diagnostic.csv             # preselected diagnostic families
  raw/metadata.json              # config, versions, seed, runtime, hashes
  summaries/summary.csv          # 108 method-condition summaries and MCSEs
  summaries/paired_power.csv     # paired differences and MC intervals
  summaries/comparison_table2.csv# 108 comparisons with rounded paper values
  summaries/validation.json      # integrity checks
  figures/01_dgp_diagnostic.png  # also PDF
  figures/02_baseline_power.png  # also PDF; includes published BH power points
  figures/03_signal_extension.png
  figures/04_fdr_calibration.png
```

In `replicates.csv`, `R,V,S,m0,m1` are counts; `FDP` and `power` are
replicate metrics. `pi0_hat_raw` is not capped. `pi0_hat_used` records any
optional cap; default cap_pi0=false. `fdr_hat_raw` and `fdr_hat` are the raw
and capped fixed-rule estimates. `alpha_rep` is the BH input, while
`bh_threshold` is its largest rejected p-value. `family_id` and replicate ID
identify shared data. The seed is reconstructed from metadata and DGP parameters.

In `summary.csv`, `fdr` is each method's achieved FDR. `true_fdr`,
`true_power`, `true_pfdr` and other analytic targets always refer to the
**fixed-threshold rule**, even where repeated on a BH row. The mean `fdr_hat`
on a BH row is its input level, not an estimate of achieved BH FDR. Empty
cells represent undefined quantities, such as power with no alternatives.

## Local validation and runtime

The implementation was tested on macOS 15.7.4 ARM64 with Python 3.11.9.
The full simulation took **59.6 seconds**, analysis 2.9 seconds and plotting
about 8.6 seconds, excluding environment setup and process startup. Allow
roughly 1-2 minutes on a similar machine; other machines can differ.

The automated tests cover hand-computed BH decisions, agreement with SciPy's
independent BH implementation, empty-rejection metrics, deterministic streams,
DGP distribution checks, analytic FDR versus independent enumeration,
paired MCSEs, cache invalidation and a complete small pipeline. The full run
passed row-level validation. See [VALIDATION.md](VALIDATION.md) for execution
evidence and [ANALYSIS.md](ANALYSIS.md) for results and limitations.

## Additional investigation of one discrepant table entry

After observing the mismatch at pi0=0.8, gamma=0.001, we ran two separate
5,000-replicate checks using a new seed, 20261016:

```bash
make analyze CONFIG=configs/check-row-fixed.json
make analyze CONFIG=configs/check-row-bernoulli.json
```

These are post hoc diagnostics, not part of the 54-condition primary study.
They gave BH powers 0.114279 (MCSE 0.000387) and 0.113573 (MCSE 0.000392),
respectively. The printed entry was checked visually in the source PDF.
We retain the published 0.129 and report an unresolved discrepancy; this does
not establish that the paper contains an error. Full raw output is retained
separately for both checks.

## Submission

The local project includes code, configuration, tests, figures and documentation.
No public repository has been created or uploaded by this implementation.
For submission, create your own public repository, verify that it is accessible,
review ANALYSIS.md in your own words, and upload the analysis and ADEMP document
as required. Generated raw files are Git-ignored but can be rebuilt with make all;
summary tables and figures can be committed.

## References

Storey, J. D. (2002). A direct approach to false discovery rates. JRSS B,
64(3), 479-498. [Paper](https://genomics.princeton.edu/storeylab/papers/directfdr.pdf).
[Author's errata](https://genomics.princeton.edu/storeylab/papers/directfdr_corrections.pdf).
Rounded Table 2 values and their provenance are in references/.
