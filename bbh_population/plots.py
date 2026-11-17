"""Optional manuscript-model diagnostics (requires the 'plots' extra)."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .extension import ExtensionParameters, event_density, integration_grid, reference_rate_shape


def make_plots(run_dir):
    """Plot reference evolution and actual generated distributions separately."""
    root = Path(run_dir)
    manifest = json.loads((root / "run.json").read_text())
    if manifest["settings"]["model"] != "paper-append" or manifest["status"] != "complete":
        raise ValueError("plots require a completed paper-append run")
    config = ExtensionParameters(**manifest["settings"]["extension"])
    low_dir = Path(manifest["original_dir"]) if manifest["original_dir"] else root / "lowz"
    bins = np.linspace(0, config.z_max, 121)
    low_hist = np.zeros(len(bins) - 1, dtype=np.int64)
    high_hist = np.zeros_like(low_hist)
    for catalog_id in manifest["catalog_ids"]:
        name = f"catalog_{catalog_id:03d}.csv"
        for directory, histogram in ((low_dir, low_hist), (root / "highz", high_hist)):
            for frame in pd.read_csv(directory / name, usecols=["z"], chunksize=100_000):
                histogram += np.histogram(frame.z, bins)[0]
    out = root / "plots"
    out.mkdir(exist_ok=True)
    z = integration_grid(0, config.z_max, config.grid_size)
    high_z = integration_grid(config.z_join, config.z_max, config.grid_size)
    high_pdf = event_density(high_z, config)
    high_pdf /= np.trapz(high_pdf, high_z)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for axis in axes:
        axis.plot(z, (1 + z) ** config.gamma, label="reference power law")
        axis.plot(z, reference_rate_shape(z, config), label="piecewise reference")
        axis.axvline(config.z_join, color="gray", linestyle=":")
        axis.set(xlabel="z", ylabel="reference R(z)/R0")
        axis.legend()
    axes[0].set_xlim(0, 4)
    axes[1].set_yscale("log")
    fig.tight_layout()
    fig.savefig(out / "reference_rate.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots()
    ax.plot(z, event_density(z, config))
    ax.axvline(config.z_join, color="gray", linestyle=":")
    ax.set(xlabel="z", ylabel="reference observer-frame rate shape (Gpc^3)")
    fig.tight_layout()
    fig.savefig(out / "reference_observer_rate.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots()
    if high_hist.sum():
        ax.stairs(high_hist / high_hist.sum() / np.diff(bins), bins, label="generated high-z")
    ax.plot(high_z, high_pdf, label="conditional target PDF")
    ax.set(xlabel="z", ylabel="p_high(z)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "highz_sampling.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots()
    for hist, label in ((low_hist, "original low-z"), (high_hist, "added high-z"),
                        (low_hist + high_hist, "combined")):
        ax.stairs(hist / len(manifest["catalog_ids"]) / np.diff(bins), bins, label=label)
    ax.axvline(config.z_join, color="gray", linestyle=":")
    ax.set(xlabel="z", ylabel="mean sources per catalog per unit z")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "combined_redshift.png", dpi=180)
    plt.close(fig)
