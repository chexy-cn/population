"""Command-line entry point; paths and resource choices are caller supplied."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .pipeline import generate_ensemble


def parse_ids(value):
    """Parse comma-separated IDs and inclusive ranges, e.g. 1,3-5."""
    ids = []
    for part in value.split(","):
        bounds = [int(item) for item in part.split("-")]
        if len(bounds) == 1:
            ids.extend(bounds)
        elif len(bounds) == 2 and bounds[0] <= bounds[1]:
            ids.extend(range(bounds[0], bounds[1] + 1))
        else:
            raise ValueError("invalid catalog ID range")
    return ids


def main():
    parser = argparse.ArgumentParser(description="Generate intrinsic GWTC-4 BP2P or manuscript catalogs.")
    commands = parser.add_subparsers(dest="command", required=True)
    generate = commands.add_parser("generate")
    generate.add_argument("--config", required=True, type=Path)
    generate.add_argument("--posterior", required=True, type=Path)
    generate.add_argument("--output-dir", required=True, type=Path)
    generate.add_argument("--original-dir", type=Path, help="Reuse original low-z CSVs and their parameter summary")
    generate.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    generate.add_argument("--catalogs", type=int, help="Size of the posterior ensemble (changes RNG mapping)")
    generate.add_argument("--seed", type=int)
    generate.add_argument("--extension-seed", type=int)
    generate.add_argument("--catalog-ids", help="Select a subset of the fixed ensemble, e.g. 1-3,7")
    validate = commands.add_parser("validate")
    validate.add_argument("--run-dir", required=True, type=Path)
    validate.add_argument("--workers", type=int, default=1)
    plot = commands.add_parser("plot", help="Reference curves and actual source histograms; requires [plots]")
    plot.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "plot":
        from .plots import make_plots
        make_plots(args.run_dir)
        return
    if args.command == "validate":
        from .validation import validate_run
        result = validate_run(args.run_dir, args.workers)
        print(json.dumps(result, indent=2))
        if not result["passed"]:
            raise SystemExit(1)
        return
    settings = json.loads(args.config.read_text())
    if set(settings) != {"model", "catalogs", "years", "low_z_max", "seed", "extension_seed", "extension"}:
        parser.error("configuration must match the documented preset schema")
    for key in ("catalogs", "seed", "extension_seed"):
        if getattr(args, key) is not None:
            settings[key] = getattr(args, key)
    generate_ensemble(settings, args.posterior, args.output_dir, args.workers,
                      parse_ids(args.catalog_ids) if args.catalog_ids else None, args.original_dir)
