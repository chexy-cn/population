"""Portable, deterministic generation of the two supported catalog ensembles."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .extension import ExtensionParameters, high_z_catalog, rate_integrals
from .sampling import generate_catalog, load_official_posterior

SOURCE_COLUMNS = (
    "z", "m1_Msun", "m2_Msun", "q", "a_1", "a_2", "cos_tilt_1", "cos_tilt_2",
    "phi_1_rad", "phi_2_rad", "t_c_detector_yr", "ra_rad", "sin_dec", "cos_iota", "psi_rad",
)


def sha256_file(path):
    """Hash a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensemble_draws(posterior, count, seed):
    """Choose distinct joint posterior rows; retain historical per-catalog RNGs."""
    if not 1 <= count <= len(posterior):
        raise ValueError(f"catalogs must lie in 1..{len(posterior)}")
    rng = np.random.default_rng(seed)
    draws = rng.choice(len(posterior), size=count, replace=False)
    seeds = rng.integers(0, np.iinfo(np.uint64).max, size=count, dtype=np.uint64)
    return [(posterior[draw], int(seed)) for draw, seed in zip(draws, seeds)]


def _low_worker(task):
    catalog_id, params, years, z_max, seed, directory = task
    frame = generate_catalog(params, years, z_max, np.random.default_rng(seed))
    path = Path(directory) / f"catalog_{catalog_id:03d}.csv"
    temporary = path.with_suffix(".csv.tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)
    return {
        "catalog_id": catalog_id, **asdict(params), "source_seed": seed,
        "expected_source_count": frame.attrs["expected_source_count"],
        "N_sources": len(frame), "observation_time_yr": years, "z_max": z_max,
    }


def _append_worker(task):
    catalog_id, params, low_path, low_count, years, config, ratio, seed, high_dir, combined_dir = task
    low_path = Path(low_path)
    if tuple(pd.read_csv(low_path, nrows=0).columns) != SOURCE_COLUMNS:
        raise ValueError(f"unexpected original schema: {low_path}")
    count = 0
    for chunk in pd.read_csv(low_path, usecols=["z"], chunksize=100_000):
        if not chunk.z.between(0, config.z_join).all():
            raise ValueError(f"original redshift outside [0, z_join]: {low_path}")
        count += len(chunk)
    if count != low_count:
        raise ValueError(f"original count disagrees with summary: {low_path}")
    rng = np.random.default_rng(seed)
    expected = count * ratio
    n_high = int(rng.poisson(expected))
    high = high_z_catalog(params, n_high, config, years, rng)
    name = low_path.name
    high_path = Path(high_dir) / name
    high_tmp = high_path.with_suffix(".csv.tmp")
    high.to_csv(high_tmp, index=False)
    high_tmp.replace(high_path)
    combined_path = Path(combined_dir) / name
    combined_tmp = combined_path.with_suffix(".csv.tmp")
    original_hash = hashlib.sha256()
    with low_path.open("rb") as source, combined_tmp.open("xb") as target:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            original_hash.update(block)
            target.write(block)
        with high_path.open("rb") as extra:
            extra.readline()  # one shared header; original bytes and ordering are retained
            shutil.copyfileobj(extra, target)
    combined_tmp.replace(combined_path)
    return {
        "catalog_id": catalog_id, "N_old": count, "N_high_expected": expected,
        "N_high": n_high, "N_total": count + n_high, **asdict(params),
        "extension_seed": seed, "original_sha256": original_hash.hexdigest(),
        "high_sha256": sha256_file(high_path), "combined_sha256": sha256_file(combined_path),
    }


def _run_tasks(worker, tasks, workers):
    rows = []
    if workers == 1:
        for index, task in enumerate(tasks, 1):
            rows.append(worker(task))
            print(f"completed {index}/{len(tasks)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as executor:
            futures = [executor.submit(worker, task) for task in tasks]
            for index, future in enumerate(as_completed(futures), 1):
                rows.append(future.result())
                print(f"completed {index}/{len(tasks)}", flush=True)
    return pd.DataFrame(rows).sort_values("catalog_id")


def generate_ensemble(settings, posterior_path, output_dir, workers=1, ids=None, original_dir=None):
    """Generate GWTC-4 low-z catalogs or the manuscript's append ensemble.

    ``catalogs`` specifies the whole posterior ensemble. ``ids`` selects work
    from that fixed ensemble, so subset work reproduces the full-run RNG mapping.
    """
    model = settings["model"]
    if model not in ("gwtc4", "paper-append"):
        raise ValueError("model must be gwtc4 or paper-append")
    if not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be positive")
    years = settings["years"]
    z_max = settings["low_z_max"]
    if not np.isfinite([years, z_max]).all() or years <= 0 or z_max <= 0:
        raise ValueError("years and low_z_max must be finite and positive")
    extension = ExtensionParameters(**settings["extension"])
    if model == "paper-append" and z_max != extension.z_join:
        raise ValueError("low_z_max must equal extension.z_join")
    if original_dir is not None and model != "paper-append":
        raise ValueError("--original-dir is only supported by paper-append")
    posterior_path, output_dir = Path(posterior_path), Path(output_dir)
    posterior = load_official_posterior(posterior_path)
    all_draws = ensemble_draws(posterior, settings["catalogs"], settings["seed"])
    selected = sorted(set(ids if ids is not None else range(1, len(all_draws) + 1)))
    if not selected or selected[0] < 1 or selected[-1] > len(all_draws):
        raise ValueError("catalog IDs must lie in 1..catalogs")

    original_table = None
    if original_dir is not None:
        original_dir = Path(original_dir).resolve()
        original_table = pd.read_csv(original_dir / "catalog_parameters_summary.csv").set_index("catalog_id")
        if not original_table.index.is_unique:
            raise ValueError("duplicate catalog IDs in original summary")
        for catalog_id in selected:
            row = original_table.loc[catalog_id]
            expected_params = asdict(all_draws[catalog_id - 1][0])
            for key, value in expected_params.items():
                if not np.isclose(row[key], value, rtol=1e-12, atol=1e-12):
                    raise ValueError(f"original {catalog_id}: posterior field {key} does not match seeded ensemble")
            if not np.isclose(row["observation_time_yr"], years) or row["z_max"] != z_max:
                raise ValueError("original duration/redshift boundary does not match configuration")
            if not (original_dir / f"catalog_{catalog_id:03d}.csv").is_file():
                raise FileNotFoundError(f"missing original catalog {catalog_id}")

    # Exclusive creation prevents both overwrites and competing writers.
    output_dir.mkdir(parents=True, exist_ok=False)
    metadata = {
        "status": "running", "settings": settings, "catalog_ids": selected,
        "posterior_sha256": sha256_file(posterior_path), "posterior_file": posterior_path.name,
        "valid_posterior_rows": len(posterior), "python": platform.python_version(),
        "versions": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "astropy", "h5py")},
        "cosmology": "astropy.cosmology.Planck18 (full sky)",
        "original_dir": str(original_dir) if original_dir else None,
    }
    manifest = output_dir / "run.json"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n")
    low_dir = original_dir or output_dir / "lowz"
    if original_table is None:
        low_dir.mkdir()
        tasks = [(i, all_draws[i - 1][0], years, z_max, all_draws[i - 1][1], low_dir) for i in selected]
        low_table = _run_tasks(_low_worker, tasks, workers)
        low_table.to_csv(low_dir / "catalog_parameters_summary.csv", index=False)
        original_table = low_table.set_index("catalog_id")
    if model == "paper-append":
        high_dir, combined_dir = output_dir / "highz", output_dir / "combined"
        high_dir.mkdir()
        combined_dir.mkdir()
        i_low, i_high = rate_integrals(extension)
        rng = np.random.default_rng(settings["extension_seed"])
        extra_seeds = rng.integers(0, np.iinfo(np.uint64).max, size=len(all_draws), dtype=np.uint64)
        tasks = [(
            i, all_draws[i - 1][0], low_dir / f"catalog_{i:03d}.csv",
            int(original_table.loc[i, "N_sources"]), years, extension,
            i_high / i_low, int(extra_seeds[i - 1]), high_dir, combined_dir,
        ) for i in selected]
        table = _run_tasks(_append_worker, tasks, workers)
        table.to_csv(combined_dir / "extension_parameters_summary.csv", index=False)
        metadata.update(I_low_Gpc3=i_low, I_high_Gpc3=i_high, high_low_ratio=i_high / i_low)
    metadata["status"] = "complete"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
