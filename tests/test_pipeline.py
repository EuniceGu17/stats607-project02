import json
import numpy as np
import pandas as pd
import pytest

from src.analyze import analyze, summarize, validate_raw
from src.clean import clean
from src.common import ROOT, load_config
from src.figures import make_figures
from src.simulation import run_simulation


def tiny_config():
    config = load_config(ROOT / "configs/smoke.json")
    return dict(config, m=40, pi0_values=[.5], gamma_values=[.01, .1], replicates=4)


def test_small_end_to_end_and_cache(tmp_path):
    config = tiny_config()
    raw_file = run_simulation(config, tmp_path)
    before = raw_file.read_bytes()
    before_mtime = raw_file.stat().st_mtime_ns
    run_simulation(config, tmp_path)
    assert raw_file.stat().st_mtime_ns == before_mtime
    raw = pd.read_csv(raw_file)
    assert len(raw) == 16
    assert validate_raw(raw, config)["status"] == "passed"
    analyze(config, tmp_path)
    make_figures(config, tmp_path)
    assert (tmp_path / "figures/01_dgp_diagnostic.png").stat().st_size > 1000
    assert (tmp_path / "summaries/comparison_table2.csv").exists()
    run_simulation(config, tmp_path, force=True)
    assert raw_file.read_bytes() == before


def test_validation_detects_corrupt_output_and_raw_hash(tmp_path):
    config = tiny_config()
    path = run_simulation(config, tmp_path)
    raw = pd.read_csv(path)
    bad = raw.copy()
    bad.loc[0, "FDP"] = .12345
    with pytest.raises(ValueError, match="FDP"):
        validate_raw(bad, config)
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="changed"):
        analyze(config, tmp_path)


def test_paired_mcse_uses_differences(tmp_path):
    config = tiny_config()
    raw = pd.read_csv(run_simulation(config, tmp_path))
    _, paired = summarize(raw, config)
    condition = raw[raw.scenario_id == paired.iloc[0].scenario_id]
    wide = condition.pivot(index="replicate_id", columns="method", values="power")
    differences = wide.fixed - wide.bh
    assert paired.iloc[0].mcse == pytest.approx(differences.std(ddof=1) / 2)


def test_config_change_invalidates_results(tmp_path):
    config = tiny_config()
    run_simulation(config, tmp_path)
    with pytest.raises(ValueError, match="stale"):
        analyze(dict(config, master_seed=config["master_seed"] + 1), tmp_path)


def test_clean_rejects_external_directory(tmp_path):
    protected = tmp_path / "important.txt"
    protected.write_text("keep")
    with pytest.raises(ValueError):
        clean(tmp_path)
    assert protected.read_text() == "keep"


def test_full_design_size():
    from src.common import design_rows
    config = load_config(ROOT / "configs/full.json")
    design = design_rows(config)
    assert len(design) == 54
    assert sum(row["experiment"] == "baseline" for row in design) == 18
