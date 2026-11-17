"""Full-row source integrity checks and streaming append verification."""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .pipeline import SOURCE_COLUMNS

BINS = np.asarray([0, 1.5, 2, 5, 10, 20, 30], dtype=float)


def compare_append(original, high, combined):
    """Check exact original prefix and exact high-z suffix (excluding its header)."""
    with Path(original).open("rb") as old, Path(high).open("rb") as extra, Path(combined).open("rb") as joined:
        for block in iter(lambda: old.read(1024 * 1024), b""):
            if joined.read(len(block)) != block:
                return False
        extra.readline()
        for block in iter(lambda: extra.read(1024 * 1024), b""):
            if joined.read(len(block)) != block:
                return False
        return joined.read(1) == b""


def inspect_sources(path, years, z_min, z_max, open_min=False):
    """Check every numeric source field, not only a sample of the first rows."""
    if tuple(pd.read_csv(path, nrows=0).columns) != SOURCE_COLUMNS:
        raise ValueError(f"unexpected schema: {path}")
    count, valid = 0, True
    bins = np.zeros(len(BINS) - 1, dtype=np.int64)
    for frame in pd.read_csv(path, chunksize=100_000):
        if frame.empty:
            continue
        count += len(frame)
        valid &= bool(np.isfinite(frame.to_numpy(dtype=float)).all())
        valid &= bool((frame.z.gt(z_min) if open_min else frame.z.ge(z_min)).all())
        valid &= bool(frame.z.le(z_max).all())
        valid &= bool((frame.m1_Msun >= frame.m2_Msun).all() and frame.m2_Msun.gt(0).all())
        valid &= bool(frame.q.gt(0).all() and frame.q.le(1).all())
        valid &= bool(np.allclose(frame.m2_Msun, frame.m1_Msun * frame.q, rtol=1e-12, atol=1e-12))
        for column in ("a_1", "a_2"):
            valid &= bool(frame[column].between(0, 1).all())
        for column in ("cos_tilt_1", "cos_tilt_2", "sin_dec", "cos_iota"):
            valid &= bool(frame[column].between(-1, 1).all())
        for column in ("phi_1_rad", "phi_2_rad", "ra_rad"):
            valid &= bool(frame[column].between(0, 2 * np.pi).all())
        valid &= bool(frame.psi_rad.between(0, np.pi).all())
        valid &= bool(frame.t_c_detector_yr.between(0, years).all())
        bins += np.histogram(frame.z, BINS)[0]
    return {"rows": count, "valid": valid, "bins": bins.tolist()}


def _inspect_task(task):
    name, low, high, combined, years, low_max, high_max, counts = task
    low_result = inspect_sources(low, years, 0, low_max)
    result = {"catalog": name, "low": low_result,
              "passed": low_result["valid"] and low_result["rows"] == counts[0]}
    if high is not None:
        high_result = inspect_sources(high, years, low_max, high_max, open_min=True)
        exact = compare_append(low, high, combined)
        result.update(high=high_result, append_byte_exact=exact)
        result["passed"] &= high_result["valid"] and high_result["rows"] == counts[1] and exact
        result["passed"] &= counts[2] == counts[0] + counts[1]
    return result


def validate_run(run_dir, workers=1):
    """Validate all selected catalogs against their run manifest and count tables."""
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "run.json").read_text())
    settings = manifest["settings"]
    low_dir = Path(manifest["original_dir"]) if manifest["original_dir"] else run_dir / "lowz"
    low_table = pd.read_csv(low_dir / "catalog_parameters_summary.csv").set_index("catalog_id")
    extended = settings["model"] == "paper-append"
    extra_table = (pd.read_csv(run_dir / "combined/extension_parameters_summary.csv").set_index("catalog_id")
                   if extended else None)
    ids = manifest["catalog_ids"]
    expected_names = {f"catalog_{i:03d}.csv" for i in ids}
    checked_dirs = [run_dir / "highz", run_dir / "combined"] if extended else []
    if not manifest["original_dir"]:
        checked_dirs.append(low_dir)
    file_sets_match = all({p.name for p in directory.glob("catalog_[0-9]*.csv")} == expected_names
                         for directory in checked_dirs)
    tasks = []
    for catalog_id in ids:
        name = f"catalog_{catalog_id:03d}.csv"
        counts = [int(low_table.loc[catalog_id, "N_sources"])]
        if extended:
            counts += [int(extra_table.loc[catalog_id, key]) for key in ("N_high", "N_total")]
            if int(extra_table.loc[catalog_id, "N_old"]) != counts[0]:
                raise ValueError("low and extension summaries disagree")
        tasks.append((name, low_dir / name, run_dir / "highz" / name if extended else None,
                      run_dir / "combined" / name if extended else None, settings["years"],
                      settings["low_z_max"], settings["extension"]["z_max"], counts))
    if workers < 1:
        raise ValueError("workers must be positive")
    if workers == 1:
        rows = [_inspect_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as executor:
            rows = list(executor.map(_inspect_task, tasks))
    totals = np.sum([row["low"]["bins"] for row in rows], axis=0)
    if extended:
        totals += np.sum([row["high"]["bins"] for row in rows], axis=0)
    low_count = sum(row["low"]["rows"] for row in rows)
    high_count = sum(row.get("high", {}).get("rows", 0) for row in rows)
    result = {
        "passed": manifest["status"] == "complete" and file_sets_match and all(row["passed"] for row in rows),
        "scope": "all source rows and all 15 numeric columns; exact append prefix AND suffix",
        "catalog_count": len(rows), "file_sets_match": file_sets_match,
        "N_low": low_count, "N_high": high_count, "N_total": low_count + high_count,
        "mean_N_total": (low_count + high_count) / len(rows),
        "redshift_bins": {f"{left:g}..{right:g}": int(n) for left, right, n in zip(BINS[:-1], BINS[1:], totals)},
        "bin_convention": "first interval includes 0; boundary 1.5 has no generated high-z rows; NumPy histogram otherwise",
        "catalogs": rows,
    }
    destination = run_dir / "validation"
    destination.mkdir(exist_ok=True)
    (destination / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
