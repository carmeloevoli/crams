#!/usr/bin/env python3
"""Plot MCMC posterior predictions against AMS-02 data.

Produces one figure per observable: absolute fluxes (H, He, C, O, Fe, …) and
ratios (B/C, C/O, …). Each panel shows the AMS-02 data, the posterior 3-sigma
band, and the posterior median model.

Usage
-----
    python plot_mcmc_posterior.py h_he_chain.npz
    python plot_mcmc_posterior.py chain.npz --nsamples 200 --output figs/posterior
"""
from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fitting import KISS_DIR, _read_kiss_table
from run_mcmc import PARAMETERS, set_halo_size
from runner import FRAGMENTATION_MODELS, CramsRunner

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"
PRELIMINARY_DIR = Path(__file__).parent / "preliminary"
BE_ISOTOPES_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_isotopes.csv"
BE_RATIOS_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_ratios.csv"

POWER = 2.7   # multiply flux by R^POWER for display (fluxes only; ratios use 0)
BAND_PERCENTILES = (0.1349898, 99.8650102)  # Gaussian-equivalent central 3 sigma
BE_ISOTOPE_R_MIN = 1.0
BE_ISOTOPE_R_MAX = 100.0

COLORS = {
    "H":  "#0072B2",
    "He": "#D55E00",
    "C":  "#8D5524",
    "N":  "#E76F9A",
    "O":  "#D62728",
    "B":  "#7E57C2",
    "Be": "#2CA02C",
    "Be7": "#009E73",
    "Be9": "#E69F00",
    "Si": "#F0E442",
    "Fe": "#5F6A6A",
}


@dataclass(frozen=True)
class Panel:
    """One observable to plot: a flux (denominator='') or a ratio."""

    numerator: str
    filename: str
    denominator: str = ""
    R_min: float = 5.0
    R_max: float = 2500.0
    color: str = "tab:blue"
    csv_path: Path | None = None
    csv_column: str = ""

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
    Panel("N",  "AMS-02_N_rigidity.txt",                color=COLORS["N"]),
    Panel("O",  "AMS-02_O_rigidity.txt",                color=COLORS["O"]),
    Panel("Si", "AMS-02_Si_rigidity.txt",               color=COLORS["Si"]),
    Panel("Fe", "AMS-02_Fe_rigidity.txt",               color=COLORS["Fe"]),
    Panel("B",  "AMS-02_B_C_rigidity.txt", denominator="C", color=COLORS["B"]),
    Panel("B",  "AMS-02_B_O_rigidity.txt", denominator="O", color=COLORS["B"]),
    Panel("C",  "AMS-02_C_O_rigidity.txt", denominator="O", color=COLORS["C"]),
    Panel("Be",  "AMS-02_Be_B_rigidity.txt", denominator="B", color=COLORS["B"]),
    # Isotope ratio: overplot preliminary AMS-02 data from mcmc/preliminary/.
    Panel("Be10", "", denominator="Be9", R_min=5.0, R_max=100.0, color=COLORS["Be"],
          csv_path=BE_RATIOS_CSV, csv_column="Be10_over_Be9"),
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
    param_names = [str(name) for name in d["param_names"]]
    acceptance  = float(d["acceptance"])
    # Fragmentation model the chain was fit with ("" if absent / crams default).
    frag = str(d["fragmentation_model"]) if "fragmentation_model" in d else ""
    # Halo half-height the chain was fit with (None if not recorded).
    halo = float(d["halo_size"]) if "halo_size" in d else None
    return chain, param_names, acceptance, (frag or None), halo


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
            spectra_list.append(s)  # fudge applied inside crams via the .ini
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


def _read_csv_data(path: Path, column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (R, y, err_lo, err_hi) from a preliminary AMS-02 CSV column."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    err_column = column[:-len("_flux_R2p7")] if column.endswith("_flux_R2p7") else column
    return (
        np.asarray(data["R_GV"], dtype=float),
        np.asarray(data[column], dtype=float),
        np.asarray(data[f"{err_column}_err_minus"], dtype=float),
        np.asarray(data[f"{err_column}_err_plus"], dtype=float),
    )


def _plot_panel(panel: Panel, spectra_median: dict, spectra_samples: list) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    power = 0.0 if panel.is_ratio else POWER

    # ── data (shown over the plotted rigidity range) ───────────────────────────
    if panel.filename or panel.csv_column:
        if panel.filename:
            x, y, err_lo, err_hi = _read_kiss_table(KISS_DIR / panel.filename)
            data_label = f"AMS-02 {panel.label}"
        else:
            if panel.csv_path is None:
                raise ValueError(f"panel {panel.tag} has no preliminary CSV path")
            x, y, err_lo, err_hi = _read_csv_data(panel.csv_path, panel.csv_column)
            data_label = "AMS-02 preliminary"

        cut = (
            np.isfinite(x) & np.isfinite(y) & np.isfinite(err_lo) & np.isfinite(err_hi)
            & (err_lo > 0) & (err_hi > 0)
            & (x >= panel.R_min) & (x <= panel.R_max)
        )
        x, y = x[cut], y[cut]
        err_lo, err_hi = err_lo[cut], err_hi[cut]

        scale = x ** power
        ax.errorbar(
            x, scale * y,
            yerr=[scale * err_lo, scale * err_hi],
            fmt="o", color=panel.color, markersize=4,
            elinewidth=1.2, capsize=0,
            label=data_label,
            zorder=3,
        )

    # ── posterior 3-sigma band ─────────────────────────────────────────────────
    if spectra_samples:
        R_model = spectra_samples[0]["R"]
        band_mask = (R_model >= panel.R_min) & (R_model <= panel.R_max)
        R_band = R_model[band_mask]
        ys = np.array([_model_observable(s, panel)[1][band_mask] for s in spectra_samples])
        lo, hi = np.nanpercentile(ys, BAND_PERCENTILES, axis=0)
        ax.fill_between(
            R_band, R_band**power * lo, R_band**power * hi,
            color=panel.color, alpha=0.25, linewidth=0,
        )

    # ── posterior median model ─────────────────────────────────────────────────
    R_model, y_med = _model_observable(spectra_median, panel)
    fit_mask = (R_model >= panel.R_min) & (R_model <= panel.R_max)
    R_fit = R_model[fit_mask]
    ax.plot(
        R_fit, R_fit**power * y_med[fit_mask],
        color=panel.color, linewidth=2.2, zorder=4,
        # Published AMS panels use the data as legend anchor; label model-only
        # and preliminary isotope panels explicitly.
        label=f"crams {ISOTOPE_LABELS.get(panel.tag, panel.label)}" if not panel.filename else None,
    )

    # ── style ──────────────────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_xlabel(r"$R$ [GV]")
    ax.set_xlim(panel.R_min, panel.R_max * 1.1)
    if panel.is_ratio:
        ax.set_ylabel(ISOTOPE_LABELS.get(panel.tag, panel.label))
    else:
        #ax.set_yscale("log")
        ax.set_ylabel(r"$R^{2.7} I\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    ax.legend()
    return fig


def _plot_be7_be9_fluxes(spectra_median: dict, spectra_samples: list) -> plt.Figure | None:
    isotopes = [
        ("Be7", "Be7_flux_R2p7", r"$^{7}\mathrm{Be}$", COLORS["Be7"]),
        ("Be9", "Be9_flux_R2p7", r"$^{9}\mathrm{Be}$", COLORS["Be9"]),
    ]
    needed = {name for name, _, _, _ in isotopes}
    if not needed.issubset(spectra_median):
        missing = sorted(needed - set(spectra_median))
        print(f"Skipping Be7/Be9 isotope posterior plot: missing spectra {missing}")
        return None

    fig, ax = plt.subplots(figsize=(12.5, 8.0))

    for name, csv_column, label, color in isotopes:
        x, y, err_lo, err_hi = _read_csv_data(BE_ISOTOPES_CSV, csv_column)
        data_mask = (
            np.isfinite(x) & np.isfinite(y) & np.isfinite(err_lo) & np.isfinite(err_hi)
            & (err_lo > 0) & (err_hi > 0)
            & (x >= BE_ISOTOPE_R_MIN) & (x <= BE_ISOTOPE_R_MAX)
        )
        ax.errorbar(
            x[data_mask], y[data_mask],
            yerr=[err_lo[data_mask], err_hi[data_mask]],
            fmt="o", color=color, markersize=4,
            elinewidth=1.2, capsize=0,
            label=rf"AMS-02 preliminary {label}",
            zorder=3,
        )

        if spectra_samples:
            R_band = spectra_samples[0]["R"]
            band_mask = (R_band >= BE_ISOTOPE_R_MIN) & (R_band <= BE_ISOTOPE_R_MAX)
            R_band = R_band[band_mask]
            ys = np.array([
                s["R"][band_mask] ** POWER * s[name][band_mask]
                for s in spectra_samples
                if name in s
            ])
            if len(ys):
                lo, hi = np.nanpercentile(ys, BAND_PERCENTILES, axis=0)
                ax.fill_between(
                    R_band, lo, hi,
                    color=color, alpha=0.22, linewidth=0,
                )

        R_model = spectra_median["R"]
        fit_mask = (R_model >= BE_ISOTOPE_R_MIN) & (R_model <= BE_ISOTOPE_R_MAX)
        ax.plot(
            R_model[fit_mask], R_model[fit_mask] ** POWER * spectra_median[name][fit_mask],
            color=color, linewidth=2.2, zorder=4,
            label=rf"crams {label}",
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"$R$ [GV]")
    ax.set_xlim(BE_ISOTOPE_R_MIN, BE_ISOTOPE_R_MAX)
    ax.set_ylabel(r"$R^{2.7} \times \Phi\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    ax.legend()
    return fig


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot MCMC posterior predictions vs AMS-02 fluxes and ratios")
    p.add_argument("chain_file",  help=".npz file from run_mcmc.py")
    p.add_argument("--nsamples",  type=int, default=100,
                   help="chain samples used for the posterior band (default: 100)")
    p.add_argument("--discard",   type=int, default=0,  help="discard first N samples")
    p.add_argument("--thin",      type=int, default=1,  help="thin chain by N")
    p.add_argument("--output",    default=None,
                   help="output prefix (default: <chain_file> stem); "
                        "figures are saved as <prefix>_<tag>.pdf")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--fragmentation-model", default=None, choices=FRAGMENTATION_MODELS,
                   help="override the fragmentation model (default: the one the "
                        "chain was fit with, recorded in the .npz)")
    p.add_argument("--halosize", type=float, default=None,
                   help="override the halo half-height h [kpc] (default: the one "
                        "the chain was fit with, recorded in the .npz)")
    p.add_argument("--seed",      type=int, default=0)
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    prefix     = Path(args.output) if args.output else chain_path.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rng        = np.random.default_rng(args.seed)

    chain, param_names, acceptance, chain_frag, chain_halo = _load_chain(
        chain_path, args.discard, args.thin)
    print(f"Chain: {chain.shape[0]:,} samples, parameters: {param_names}")
    print(f"Acceptance fraction: {acceptance:.3f}")

    # Plot with the model/halo the chain was fit with unless the user overrides.
    frag_model = args.fragmentation_model or chain_frag
    halo = args.halosize if args.halosize is not None else chain_halo
    if halo is not None:
        set_halo_size(halo)
    halo_size = next(p.value for p in PARAMETERS if p.name == "h")
    print(f"Fragmentation model: {frag_model or 'crams default'}   halo h: {halo_size} kpc")
    runner = CramsRunner(build_dir=args.build_dir, read_isotopes=True,
                         fragmentation_model=frag_model)

    print("Running posterior median model…")
    ini_median      = _median_params(chain, param_names)
    spectra_median  = runner.run(ini_median)
    if spectra_median is None:
        sys.exit("Median parameter run failed.")
    # fudge applied inside crams via the .ini; no post-scaling needed

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

    extra_plots = [
        ("Be7_Be9", _plot_be7_be9_fluxes(spectra_median, spectra_samples)),
    ]
    for tag, fig in extra_plots:
        if fig is None:
            continue
        out = prefix.parent / f"{prefix.name}_{tag}.pdf"
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"Saved: {out.resolve()}")


if __name__ == "__main__":
    main()
