#!/usr/bin/env python3
"""Plot the MCMC best-fit model and posterior band against AMS-02 H and He data.

Usage
-----
    python plot_bestfit.py h_he_chain.npz
    python plot_bestfit.py h_he_chain.npz --nsamples 200 --output bestfit.pdf
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fitting import KISS_DIR, _read_kiss_table, unpack_theta
from run_mcmc import PARAMETERS
from runner import CramsRunner

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"

SPECIES = [
    dict(symbol="H",  filename="AMS-02_H_rigidity.txt",
         color="#0072B2", label=r"H",  R_min=10.0, R_max=1500.0),
    dict(symbol="He", filename="AMS-02_He_rigidity.txt",
         color="#D55E00", label=r"He", R_min=10.0, R_max=2500.0),
]

POWER = 2.7   # multiply flux by R^POWER for display


def _load_chain(path: Path, discard: int = 0, thin: int = 1):
    d = np.load(path, allow_pickle=True)
    chain = d["chain"]
    if discard:
        chain = chain[discard:]
    if thin > 1:
        chain = chain[::thin]
    param_names = list(d["param_names"])
    acceptance  = float(d["acceptance"])
    return chain, param_names, acceptance


def _median_params(chain: np.ndarray, param_names: list[str]) -> dict[str, float]:
    medians = np.median(chain, axis=0)
    active_vals = dict(zip(param_names, medians))
    # merge with fixed parameters from PARAMETERS
    ini = {p.name: p.value for p in PARAMETERS}
    ini.update(active_vals)
    return ini


def _run_sample(runner: CramsRunner, chain: np.ndarray,
                param_names: list[str], n: int, rng) -> list[dict]:
    """Run crams for *n* random samples drawn from *chain*."""
    idx = rng.choice(len(chain), size=n, replace=False)
    spectra_list = []
    for i in idx:
        vals = dict(zip(param_names, chain[i]))
        ini  = {p.name: p.value for p in PARAMETERS}
        ini.update(vals)
        s = runner.run(ini)
        if s is not None:
            spectra_list.append(s)
    return spectra_list


def _plot_species(ax, sp: dict, spectra_median, spectra_samples: list) -> None:
    sym      = sp["symbol"]
    color    = sp["color"]
    R_min    = sp["R_min"]
    R_max    = sp["R_max"]

    # ── data ──────────────────────────────────────────────────────────────────
    x, y, err_lo, err_hi = _read_kiss_table(KISS_DIR / sp["filename"])
    cut = (x >= R_min) & (x >= 1.0)   # show data from 1 GV for context
    x, y = x[cut], y[cut]
    err_lo, err_hi = err_lo[cut], err_hi[cut]

    scale = x ** POWER
    ax.errorbar(
        x, scale * y,
        yerr=[scale * err_lo, scale * err_hi],
        fmt="o", color=color, markersize=4,
        elinewidth=1.2, capsize=0,
        label=r"AMS-02 " + sp["label"],
        zorder=3,
    )

    # ── posterior band ────────────────────────────────────────────────────────
    if spectra_samples:
        R_model = spectra_samples[0]["R"]
        band_mask = (R_model >= R_min) & (R_model <= R_max)
        fluxes = np.array([s[sym][band_mask] for s in spectra_samples])
        R_band = R_model[band_mask]
        lo = np.percentile(fluxes, 16, axis=0)
        hi = np.percentile(fluxes, 84, axis=0)
        ax.fill_between(
            R_band, R_band**POWER * lo, R_band**POWER * hi,
            color=color, alpha=0.25, linewidth=0,
        )

    # ── best-fit (median) ─────────────────────────────────────────────────────
    R_model = spectra_median["R"]
    fit_mask = (R_model >= R_min) & (R_model <= R_max)
    R_fit = R_model[fit_mask]
    ax.plot(
        R_fit, R_fit**POWER * spectra_median[sym][fit_mask],
        color=color, linewidth=2.2, zorder=4,
    )


def _apply_style(ax, sp: dict) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ [GV]")
    ax.set_ylabel(r"$R^{2.7} \times \Phi\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    ax.legend()
    ax.set_xlim(8, sp["R_max"] * 1.1)


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot MCMC best-fit vs AMS-02 H and He")
    p.add_argument("chain_file",  help=".npz file from run_mcmc.py")
    p.add_argument("--nsamples",  type=int, default=100,
                   help="chain samples used for the posterior band (default: 100)")
    p.add_argument("--discard",   type=int, default=0,  help="discard first N samples")
    p.add_argument("--thin",      type=int, default=1,  help="thin chain by N")
    p.add_argument("--output",    default=None,
                   help="output file (default: <chain_file>.pdf)")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--seed",      type=int, default=0)
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    output     = Path(args.output) if args.output else chain_path.with_suffix(".pdf")
    rng        = np.random.default_rng(args.seed)

    chain, param_names, acceptance = _load_chain(chain_path, args.discard, args.thin)
    print(f"Chain: {chain.shape[0]:,} samples, parameters: {param_names}")
    print(f"Acceptance fraction: {acceptance:.3f}")

    runner = CramsRunner(build_dir=args.build_dir)

    print("Running median best-fit…")
    ini_median      = _median_params(chain, param_names)
    spectra_median  = runner.run(ini_median)
    if spectra_median is None:
        sys.exit("Median parameter run failed.")

    nsamples = min(args.nsamples, len(chain))
    print(f"Running {nsamples} posterior samples for band…")
    spectra_samples = _run_sample(runner, chain, param_names, nsamples, rng)
    print(f"  {len(spectra_samples)} successful runs")

    # ── plot ──────────────────────────────────────────────────────────────────
    use_tex = STYLE.exists() and shutil.which("latex") is not None
    if STYLE.exists():
        try:
            plt.style.use(str(STYLE))
            if not use_tex:
                plt.rcParams["text.usetex"] = False
        except Exception:
            pass

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.4))

    for ax, sp in zip(axes, SPECIES):
        _plot_species(ax, sp, spectra_median, spectra_samples)
        _apply_style(ax, sp)

    fig.savefig(output, bbox_inches="tight", dpi=150)
    print(f"Saved: {output.resolve()}")


if __name__ == "__main__":
    main()
