#!/usr/bin/env python3
"""Plot the MCMC best-fit model and posterior band against AMS-02 data.

Produces one figure per observable: absolute fluxes (H, He, C, O, Fe, …) and
ratios (B/C, C/O, …). Each panel shows the AMS-02 data, the posterior 68% band,
and the median best-fit model.

Usage
-----
    python plot_bestfit.py h_he_chain.npz
    python plot_bestfit.py chain.npz --nsamples 200 --output figs/bestfit
"""
from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fitting import KISS_DIR, _read_kiss_table, unpack_theta
from run_mcmc import PARAMETERS
from runner import CramsRunner

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"

POWER = 2.7   # multiply flux by R^POWER for display (fluxes only; ratios use 0)

COLORS = {
    "H":  "#0072B2",
    "He": "#D55E00",
    "C":  "#8D5524",
    "N":  "#E76F9A",
    "O":  "#D62728",
    "B":  "#7E57C2",
    "Be": "#2CA02C",
    "Fe": "#5F6A6A",
}


@dataclass(frozen=True)
class Panel:
    """One observable to plot: a flux (denominator='') or a ratio."""

    numerator: str
    filename: str
    denominator: str = ""
    R_min: float = 10.0
    R_max: float = 2500.0
    color: str = "tab:blue"

    @property
    def is_ratio(self) -> bool:
        return bool(self.denominator)

    @property
    def label(self) -> str:
        return f"{self.numerator}/{self.denominator}" if self.is_ratio else self.numerator

    @property
    def tag(self) -> str:
        return f"{self.numerator}_{self.denominator}" if self.is_ratio else self.numerator


PANELS = [
    Panel("H",  "AMS-02_H_rigidity.txt",  R_max=1500.0, color=COLORS["H"]),
    Panel("He", "AMS-02_He_rigidity.txt",               color=COLORS["He"]),
    Panel("C",  "AMS-02_C_rigidity.txt",                color=COLORS["C"]),
    Panel("O",  "AMS-02_O_rigidity.txt",                color=COLORS["O"]),
    Panel("Fe", "AMS-02_Fe_rigidity.txt",               color=COLORS["Fe"]),
    Panel("B",  "AMS-02_B_C_rigidity.txt", denominator="C", color=COLORS["B"]),
    Panel("C",  "AMS-02_C_O_rigidity.txt", denominator="O", color=COLORS["C"]),
    # Isotope ratio: model-only posterior (no published AMS-02 data table).
    Panel("Be10", "", denominator="Be9", R_min=2.0, R_max=100.0, color=COLORS["Be"]),
]

# Pretty axis labels for ratios whose generic "num/den" form is not ideal.
ISOTOPE_LABELS = {
    "Be10_Be9": r"$^{10}\mathrm{Be}/^{9}\mathrm{Be}$",
}


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


def _model_observable(spectra: dict, panel: Panel) -> tuple[np.ndarray, np.ndarray]:
    """Return (R, y) for the model: flux, or numerator/denominator for a ratio."""
    R = spectra["R"]
    num = spectra[panel.numerator]
    if panel.is_ratio:
        den = spectra[panel.denominator]
        y = np.full_like(num, np.nan)
        np.divide(num, den, out=y, where=den > 0)
    else:
        y = num
    return R, y


def _plot_panel(panel: Panel, spectra_median: dict, spectra_samples: list) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    power = 0.0 if panel.is_ratio else POWER

    # ── data (shown from 1 GV for context, up to R_max) ────────────────────────
    if panel.filename:
        x, y, err_lo, err_hi = _read_kiss_table(KISS_DIR / panel.filename)
        cut = (x >= 1.0) & (x <= panel.R_max)
        x, y = x[cut], y[cut]
        err_lo, err_hi = err_lo[cut], err_hi[cut]

        scale = x ** power
        ax.errorbar(
            x, scale * y,
            yerr=[scale * err_lo, scale * err_hi],
            fmt="o", color=panel.color, markersize=4,
            elinewidth=1.2, capsize=0,
            label=f"AMS-02 {panel.label}",
            zorder=3,
        )

    # ── posterior 68% band ─────────────────────────────────────────────────────
    if spectra_samples:
        R_model = spectra_samples[0]["R"]
        band_mask = (R_model >= panel.R_min) & (R_model <= panel.R_max)
        R_band = R_model[band_mask]
        ys = np.array([_model_observable(s, panel)[1][band_mask] for s in spectra_samples])
        lo = np.nanpercentile(ys, 16, axis=0)
        hi = np.nanpercentile(ys, 84, axis=0)
        ax.fill_between(
            R_band, R_band**power * lo, R_band**power * hi,
            color=panel.color, alpha=0.25, linewidth=0,
        )

    # ── best-fit (median) ──────────────────────────────────────────────────────
    R_model, y_med = _model_observable(spectra_median, panel)
    fit_mask = (R_model >= panel.R_min) & (R_model <= panel.R_max)
    R_fit = R_model[fit_mask]
    ax.plot(
        R_fit, R_fit**power * y_med[fit_mask],
        color=panel.color, linewidth=2.2, zorder=4,
        # label the model only when there is no data series to anchor the legend
        label=f"crams {ISOTOPE_LABELS.get(panel.tag, panel.label)}" if not panel.filename else None,
    )

    # ── style ──────────────────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_xlabel(r"$R$ [GV]")
    ax.set_xlim(1.0, panel.R_max * 1.1)
    if panel.is_ratio:
        ax.set_ylabel(ISOTOPE_LABELS.get(panel.tag, panel.label))
    else:
        ax.set_yscale("log")
        ax.set_ylabel(r"$R^{2.7} \times \Phi\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    ax.legend()
    return fig


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot MCMC best-fit vs AMS-02 fluxes and ratios")
    p.add_argument("chain_file",  help=".npz file from run_mcmc.py")
    p.add_argument("--nsamples",  type=int, default=100,
                   help="chain samples used for the posterior band (default: 100)")
    p.add_argument("--discard",   type=int, default=0,  help="discard first N samples")
    p.add_argument("--thin",      type=int, default=1,  help="thin chain by N")
    p.add_argument("--output",    default=None,
                   help="output prefix (default: <chain_file> stem); "
                        "each panel is saved as <prefix>_<species>.pdf")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--seed",      type=int, default=0)
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    prefix     = Path(args.output) if args.output else chain_path.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rng        = np.random.default_rng(args.seed)

    chain, param_names, acceptance = _load_chain(chain_path, args.discard, args.thin)
    print(f"Chain: {chain.shape[0]:,} samples, parameters: {param_names}")
    print(f"Acceptance fraction: {acceptance:.3f}")

    runner = CramsRunner(build_dir=args.build_dir, read_isotopes=True)

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

    for panel in PANELS:
        fig = _plot_panel(panel, spectra_median, spectra_samples)
        out = prefix.parent / f"{prefix.name}_{panel.tag}.pdf"
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"Saved: {out.resolve()}")


if __name__ == "__main__":
    main()
