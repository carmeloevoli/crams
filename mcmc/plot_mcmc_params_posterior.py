#!/usr/bin/env python3
"""Plot MCMC parameter posterior diagnostics.

Produces posterior plots for the three injection slopes, halo size when sampled,
and the two-dimensional D0/H versus delta posterior.

Usage
-----
    python plot_mcmc_params_posterior.py h_he_chain.npz
    python plot_mcmc_params_posterior.py chain.npz --output figs/params_posterior
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"
POSTERIOR_BINS = 70

COLORS = {
    "H":  "#0072B2",
    "He": "#D55E00",
    "C":  "#8D5524",
}

SLOPE_POSTERIORS = [
    ("hslope",  r"$\gamma_\mathrm{H}$",       COLORS["H"]),
    ("heslope", r"$\gamma_\mathrm{He}$",      COLORS["He"]),
    ("slope",   r"$\gamma_{Z\geq3}$",         COLORS["C"]),
]


def _load_chain(path: Path, discard: int = 0, thin: int = 1):
    d = np.load(path, allow_pickle=True)
    chain = d["chain"]
    if discard:
        chain = chain[discard:]
    if thin > 1:
        chain = chain[::thin]
    param_names = [str(name) for name in d["param_names"]]
    acceptance = float(d["acceptance"])
    return chain, param_names, acceptance


def _param_columns(
    param_names: list[str],
    required: list[str],
    plot_name: str,
) -> dict[str, int] | None:
    missing = [name for name in required if name not in param_names]
    if missing:
        print(f"Skipping {plot_name}: missing chain parameter(s) {missing}")
        return None
    return {name: param_names.index(name) for name in required}


def _finite_column(chain: np.ndarray, column: int) -> np.ndarray:
    values = chain[:, column]
    return values[np.isfinite(values)]


def _plot_slope_posteriors(chain: np.ndarray, param_names: list[str]) -> plt.Figure | None:
    names = [name for name, _, _ in SLOPE_POSTERIORS]
    columns = _param_columns(param_names, names, "slope posterior plot")
    if columns is None:
        return None

    fig, ax = plt.subplots(figsize=(12.5, 10.5))
    for name, label, color in SLOPE_POSTERIORS:
        values = _finite_column(chain, columns[name])
        if len(values) == 0:
            print(f"Skipping {label}: no finite samples")
            continue
        lo, median, hi = np.nanpercentile(values, [16, 50, 84])
        ax.hist(
            values,
            bins=POSTERIOR_BINS,
            density=True,
            histtype="step",
            linewidth=2.0,
            color=color,
            label=rf"{label}: {median:.3f}$^{{+{hi - median:.3f}}}_{{-{median - lo:.3f}}}$",
        )
        ax.axvline(median, color=color, linestyle="--", linewidth=1.2, alpha=0.85)

    ax.set_xlabel("Injection spectral index")
    ax.set_ylabel("Posterior density")
    ax.legend()
    fig.tight_layout()
    return fig


def _plot_halo_size_posterior(chain: np.ndarray, param_names: list[str]) -> plt.Figure | None:
    columns = _param_columns(param_names, ["h"], "Halo Size posterior plot")
    if columns is None:
        return None

    h = _finite_column(chain, columns["h"])
    if len(h) == 0:
        print("Skipping Halo Size posterior plot: no finite samples")
        return None

    lo, median, hi = np.nanpercentile(h, [16, 50, 84])
    color = "tab:green"

    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    ax.hist(
        h,
        bins=POSTERIOR_BINS,
        density=True,
        histtype="stepfilled",
        linewidth=1.6,
        edgecolor=color,
        facecolor=color,
        alpha=0.28,
        label=rf"$H = {median:.2f}^{{+{hi - median:.2f}}}_{{-{median - lo:.2f}}}\ \mathrm{{kpc}}$",
    )
    ax.axvline(median, color=color, linestyle="--", linewidth=1.4, alpha=0.9)
    ax.set_xlabel(r"Halo Size $H$ [kpc]")
    ax.set_ylabel("Posterior density")
    ax.legend()
    fig.tight_layout()
    return fig


def _credible_count_levels(counts: np.ndarray, fractions: tuple[float, ...]) -> list[float]:
    flat = counts.ravel()
    flat = flat[flat > 0]
    if len(flat) < 2:
        return []

    ordered = np.sort(flat)[::-1]
    cdf = np.cumsum(ordered) / np.sum(ordered)
    levels = []
    for fraction in fractions:
        idx = min(np.searchsorted(cdf, fraction), len(ordered) - 1)
        levels.append(float(ordered[idx]))

    max_count = float(np.max(flat))
    return sorted({level for level in levels if 0.0 < level < max_count})


def _plot_d0h_delta_posterior(chain: np.ndarray, param_names: list[str]) -> plt.Figure | None:
    columns = _param_columns(param_names, ["d0_h", "delta"], "D0/H vs delta posterior plot")
    if columns is None:
        return None

    d0_h = chain[:, columns["d0_h"]]
    delta = chain[:, columns["delta"]]
    finite = np.isfinite(d0_h) & np.isfinite(delta)
    d0_h, delta = d0_h[finite], delta[finite]
    if len(d0_h) == 0:
        print("Skipping D0/H vs delta posterior plot: no finite samples")
        return None

    counts, x_edges, y_edges = np.histogram2d(d0_h, delta, bins=POSTERIOR_BINS)
    posterior_mass = counts / np.sum(counts)
    posterior_mass_percent = 100.0 * posterior_mass
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    fig, ax = plt.subplots(figsize=(12.5, 10.5))
    mesh = ax.pcolormesh(x_edges, y_edges, posterior_mass_percent.T, cmap="Blues", shading="auto")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Posterior mass per bin")
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=100.0))

    levels = _credible_count_levels(counts, (0.68, 0.95, 0.997))
    if levels:
        ax.contour(x_centers, y_centers, counts.T, levels=levels,
                   colors="black", linewidths=1.1, alpha=0.75)

    ax.axvline(np.nanmedian(d0_h), color="tab:red", linestyle="--", linewidth=1.1)
    ax.axhline(np.nanmedian(delta), color="tab:red", linestyle="--", linewidth=1.1)
    ax.set_xlabel(r"$D_0/H$ [$10^{28}\,\mathrm{cm^2\,s^{-1}\,kpc^{-1}}$]")
    ax.set_ylabel(r"$\delta$")
    fig.tight_layout()
    return fig


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot MCMC parameter posterior diagnostics")
    p.add_argument("chain_file", help=".npz file from run_mcmc.py")
    p.add_argument("--discard", type=int, default=0, help="discard first N samples")
    p.add_argument("--thin", type=int, default=1, help="thin chain by N")
    p.add_argument("--output", default=None,
                   help="output prefix (default: <chain_file> stem); "
                        "figures are saved as <prefix>_<tag>.pdf")
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    prefix = Path(args.output) if args.output else chain_path.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    chain, param_names, acceptance = _load_chain(chain_path, args.discard, args.thin)
    print(f"Chain: {chain.shape[0]:,} samples, parameters: {param_names}")
    print(f"Acceptance fraction: {acceptance:.3f}")

    use_tex = STYLE.exists() and shutil.which("latex") is not None
    if STYLE.exists():
        try:
            plt.style.use(str(STYLE))
            if not use_tex:
                plt.rcParams["text.usetex"] = False
        except Exception:
            pass

    plots = [
        ("slopes", _plot_slope_posteriors(chain, param_names)),
        ("halo_size", _plot_halo_size_posterior(chain, param_names)),
        ("d0h_delta", _plot_d0h_delta_posterior(chain, param_names)),
    ]
    for tag, fig in plots:
        if fig is None:
            continue
        out = prefix.parent / f"{prefix.name}_{tag}.pdf"
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"Saved: {out.resolve()}")


if __name__ == "__main__":
    main()
