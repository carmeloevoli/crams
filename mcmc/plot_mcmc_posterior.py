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
BE_ISOTOPES_ECRS_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_isotopes_ECRS.csv"
BE_RATIOS_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_ratios.csv"

POWER = 2.7   # multiply flux by R^POWER for display (fluxes only; ratios use 0)
BAND_PERCENTILES = (0.1349898, 99.8650102)  # Gaussian-equivalent central 3 sigma
BE_ISOTOPE_R_MIN = 1.0
BE_ISOTOPE_R_MAX = 100.0
DISPLAY_R_MIN = 2.0   # show data + model down to here …
FIT_R_MIN = 5.0       # … but the fit uses only R >= FIT_R_MIN; shade the 2–5 GV gap

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
    "Be10": "#CC79A7",
    "Si": "#F0E442",
    "S":  "#17BECF",
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
    Panel("N",  "AMS-02_N_rigidity.txt",  R_min=30.0,   color=COLORS["N"]),
    Panel("O",  "AMS-02_O_rigidity.txt",                color=COLORS["O"]),
    Panel("Si", "AMS-02_Si_rigidity.txt", R_min=30.0,   color=COLORS["Si"]),
    Panel("S",  "AMS-02_S_rigidity.txt",  R_min=30.0,   color=COLORS["S"]),
    Panel("Fe", "AMS-02_Fe_rigidity.txt", R_min=30.0,   color=COLORS["Fe"]),
    Panel("Be", "AMS-02_Be_rigidity.txt",               color=COLORS["Be"]),
    Panel("B",  "AMS-02_B_C_rigidity.txt", denominator="C", color=COLORS["B"]),
    Panel("B",  "AMS-02_B_O_rigidity.txt", denominator="O", color=COLORS["B"]),
    Panel("C",  "AMS-02_C_O_rigidity.txt", denominator="O", color=COLORS["C"]),
    Panel("Be",  "AMS-02_Be_B_rigidity.txt", denominator="B", R_max=1000.0, color=COLORS["B"]),
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
    # Per-sample log-probability, discarded/thinned in lockstep with the chain so
    # the MAP point can be located (None if the chain file predates log_prob).
    log_prob = d["log_prob"] if "log_prob" in d else None
    if discard:
        chain = chain[discard:]
        if log_prob is not None:
            log_prob = log_prob[discard:]
    if thin > 1:
        chain = chain[::thin]
        if log_prob is not None:
            log_prob = log_prob[::thin]
    param_names = [str(name) for name in d["param_names"]]
    acceptance  = float(d["acceptance"])
    # Fragmentation model the chain was fit with ("" if absent / crams default).
    frag = str(d["fragmentation_model"]) if "fragmentation_model" in d else ""
    # Halo half-height the chain was fit with (None if not recorded).
    halo = float(d["halo_size"]) if "halo_size" in d else None
    return chain, param_names, acceptance, (frag or None), halo, log_prob


def _median_params(chain: np.ndarray, param_names: list[str]) -> dict[str, float]:
    medians = np.median(chain, axis=0)
    active_vals = dict(zip(param_names, medians))
    # merge with fixed parameters from PARAMETERS
    ini = {p.name: p.value for p in PARAMETERS}
    ini.update(active_vals)
    return ini


def _map_params(chain: np.ndarray, log_prob: np.ndarray,
                param_names: list[str]) -> dict[str, float]:
    """Full crams parameter dict at the maximum-a-posteriori chain sample."""
    active_vals = dict(zip(param_names, chain[int(np.argmax(log_prob))]))
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


def _plot_mask(R: np.ndarray, R_min: float, R_max: float) -> np.ndarray:
    """Mask selecting the model nodes to draw between R_min and R_max.

    crams outputs on a fixed log grid (1 GV … 10 TeV, 100 nodes), so the first
    node with R >= R_min falls slightly above R_min (e.g. 5.34 GV for R_min=5),
    which is why a plain ``R >= R_min`` cut makes the model curve start visibly
    to the right of the axis. Include the last node just below R_min as well so
    the drawn curve reaches the left axis edge; set_xlim clips it back at R_min.
    """
    mask = (R >= R_min) & (R <= R_max)
    below = np.flatnonzero(R < R_min)
    if below.size:
        mask[below[-1]] = True
    return mask


def _model_at(x_data: np.ndarray, R_model: np.ndarray, y_model: np.ndarray) -> np.ndarray:
    """Interpolate the model curve (on its R grid) onto the data rigidities."""
    m = np.isfinite(R_model) & np.isfinite(y_model)
    return np.interp(np.log(x_data), np.log(R_model[m]), y_model[m])


def _shade_excluded(axes, disp_min: float, fit_min: float, label=None) -> None:
    """Shade [disp_min, fit_min] on each axis: data shown but excluded from the fit."""
    if fit_min <= disp_min:
        return
    for i, a in enumerate(axes):
        a.axvspan(disp_min, fit_min, color="0.5", alpha=0.15, linewidth=0,
                  zorder=0, label=label if i == 0 else None)


def _plot_panel(panel: Panel, spectra_median: dict, spectra_samples: list,
                spectra_map: dict | None = None) -> plt.Figure:
    power = 0.0 if panel.is_ratio else POWER
    # Display down to DISPLAY_R_MIN but never above the panel's own fit floor.
    disp_min = min(DISPLAY_R_MIN, panel.R_min)

    fig, (ax, axr) = plt.subplots(
        2, 1, sharex=True, figsize=(12.5, 9.5),
        height_ratios=[3, 1], layout="constrained",
    )
    fig.get_layout_engine().set(hspace=0.03)

    # median model on its native grid; reused for the curve and the residuals
    R_model, y_med = _model_observable(spectra_median, panel)

    # ── data (shown over the plotted rigidity range) ───────────────────────────
    x = y = err_lo = err_hi = None
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
            & (x >= disp_min) & (x <= panel.R_max)
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

    # ── posterior 3-sigma band (top panel + fractional residual band) ──────────
    band_mask = med_on_band = R_band = lo = hi = None
    if spectra_samples:
        band_mask = _plot_mask(R_model, disp_min, panel.R_max)
        R_band = R_model[band_mask]
        med_on_band = y_med[band_mask]
        ys = np.array([_model_observable(s, panel)[1][band_mask] for s in spectra_samples])
        lo, hi = np.nanpercentile(ys, BAND_PERCENTILES, axis=0)
        ax.fill_between(
            R_band, R_band**power * lo, R_band**power * hi,
            color=panel.color, alpha=0.25, linewidth=0,
        )

    # ── posterior median (solid) and MAP (dotted) models ───────────────────────
    fit_mask = _plot_mask(R_model, disp_min, panel.R_max)
    R_fit = R_model[fit_mask]
    # Published AMS panels use the data as legend anchor and leave the model
    # unlabelled; label model-only/preliminary panels, and always distinguish
    # median from MAP once both are drawn.
    model_label = ISOTOPE_LABELS.get(panel.tag, panel.label)
    base = "crams" if panel.filename else f"crams {model_label}"
    show_label = (not panel.filename) or (spectra_map is not None)
    ax.plot(
        R_fit, R_fit**power * y_med[fit_mask],
        color=panel.color, linewidth=2.2, zorder=4,
        label=(f"{base} (median)" if spectra_map is not None else base) if show_label else None,
    )
    if spectra_map is not None:
        R_map, y_map = _model_observable(spectra_map, panel)
        map_mask = _plot_mask(R_map, disp_min, panel.R_max)
        ax.plot(
            R_map[map_mask], R_map[map_mask]**power * y_map[map_mask],
            color=panel.color, linewidth=2.0, linestyle=":", zorder=4,
            label=f"{base} (MAP)",
        )

    # ── residual panel: (data - model) / model ──────────────────────────────────
    if med_on_band is not None:
        ok = med_on_band > 0
        axr.fill_between(
            R_band[ok], (lo / med_on_band - 1.0)[ok], (hi / med_on_band - 1.0)[ok],
            color=panel.color, alpha=0.25, linewidth=0,
        )
    if x is not None and len(x):
        m_at = _model_at(x, R_model, y_med)
        good = m_at > 0
        axr.errorbar(
            x[good], y[good] / m_at[good] - 1.0,
            yerr=[err_lo[good] / m_at[good], err_hi[good] / m_at[good]],
            fmt="o", color=panel.color, markersize=4,
            elinewidth=1.2, capsize=0, zorder=3,
        )
    axr.axhline(0.0, color="k", linewidth=0.8, zorder=1)

    # ── shade the sub-threshold region excluded from the fit ────────────────────
    _shade_excluded((ax, axr), disp_min, panel.R_min,
                    label=rf"$R < {panel.R_min:g}$ GV (not fitted)")

    # ── style ──────────────────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_xlim(disp_min, panel.R_max * 1.1)
    if panel.is_ratio:
        ax.set_ylabel(ISOTOPE_LABELS.get(panel.tag, panel.label))
    else:
        #ax.set_yscale("log")
        ax.set_ylabel(r"$R^{2.7} I\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    axr.set_xlabel(r"$R$ [GV]")
    axr.set_ylabel("residual")
    ax.legend()
    return fig


# Be isotope flux panels, drawn from the preliminary ECRS release. Each entry is
# (spectrum key, CSV column, legend label, colour).
BE7 = ("Be7", "Be7_flux_R2p7", r"$^{7}\mathrm{Be}$", COLORS["Be7"])
BE9 = ("Be9", "Be9_flux_R2p7", r"$^{9}\mathrm{Be}$", COLORS["Be9"])
BE10 = ("Be10", "Be10_flux_R2p7", r"$^{10}\mathrm{Be}$", COLORS["Be10"])


def _plot_isotope_fluxes(isotopes: list, spectra_median: dict, spectra_samples: list,
                         spectra_map: dict | None = None,
                         r_min: float = FIT_R_MIN) -> plt.Figure | None:
    """Plot R^2.7-scaled Be-isotope fluxes vs the preliminary AMS-02 ECRS data.

    *r_min* is the fit lower bound: data are shown down to DISPLAY_R_MIN and the
    excluded [DISPLAY_R_MIN, r_min] region is shaded.
    """
    needed = {name for name, _, _, _ in isotopes}
    if not needed.issubset(spectra_median):
        missing = sorted(needed - set(spectra_median))
        print(f"Skipping {sorted(needed)} isotope posterior plot: missing spectra {missing}")
        return None

    disp_min = min(DISPLAY_R_MIN, r_min)
    fig, (ax, axr) = plt.subplots(
        2, 1, sharex=True, figsize=(12.5, 9.5),
        height_ratios=[3, 1], layout="constrained",
    )
    fig.get_layout_engine().set(hspace=0.03)
    R_model = spectra_median["R"]

    for name, csv_column, label, color in isotopes:
        # Previous (non-ECRS) preliminary release, drawn in gray for comparison.
        xo, yo, elo, ehi = _read_csv_data(BE_ISOTOPES_CSV, csv_column)
        old_mask = (
            np.isfinite(xo) & np.isfinite(yo) & np.isfinite(elo) & np.isfinite(ehi)
            & (elo > 0) & (ehi > 0)
            & (xo >= disp_min) & (xo <= BE_ISOTOPE_R_MAX)
        )
        ax.errorbar(
            xo[old_mask], yo[old_mask],
            yerr=[elo[old_mask], ehi[old_mask]],
            fmt="o", color="0.6", markersize=4,
            elinewidth=1.2, capsize=0,
            label="_nolegend_",
            zorder=2,
        )

        x, y, err_lo, err_hi = _read_csv_data(BE_ISOTOPES_ECRS_CSV, csv_column)
        data_mask = (
            np.isfinite(x) & np.isfinite(y) & np.isfinite(err_lo) & np.isfinite(err_hi)
            & (err_lo > 0) & (err_hi > 0)
            & (x >= disp_min) & (x <= BE_ISOTOPE_R_MAX)
        )
        x, y = x[data_mask], y[data_mask]
        err_lo, err_hi = err_lo[data_mask], err_hi[data_mask]
        ax.errorbar(
            x, y, yerr=[err_lo, err_hi],
            fmt="o", color=color, markersize=4,
            elinewidth=1.2, capsize=0,
            label=rf"AMS-02 preliminary {label} (ECRS)",
            zorder=3,
        )

        # median model (R^2.7-scaled) on its native grid, for curve + residuals
        y_med = R_model ** POWER * spectra_median[name]

        band_mask = med_on_band = R_band = lo = hi = None
        if spectra_samples:
            band_mask = _plot_mask(R_model, disp_min, BE_ISOTOPE_R_MAX)
            R_band = R_model[band_mask]
            med_on_band = y_med[band_mask]
            ys = np.array([
                s["R"][band_mask] ** POWER * s[name][band_mask]
                for s in spectra_samples
                if name in s
            ])
            if len(ys):
                lo, hi = np.nanpercentile(ys, BAND_PERCENTILES, axis=0)
                ax.fill_between(R_band, lo, hi, color=color, alpha=0.22, linewidth=0)

        fit_mask = _plot_mask(R_model, disp_min, BE_ISOTOPE_R_MAX)
        ax.plot(
            R_model[fit_mask], y_med[fit_mask],
            color=color, linewidth=2.2, zorder=4,
            label=rf"crams {label}" + (" (median)" if spectra_map is not None else ""),
        )
        if spectra_map is not None and name in spectra_map:
            R_map = spectra_map["R"]
            y_map = R_map ** POWER * spectra_map[name]
            map_mask = _plot_mask(R_map, disp_min, BE_ISOTOPE_R_MAX)
            ax.plot(
                R_map[map_mask], y_map[map_mask],
                color=color, linewidth=2.0, linestyle=":", zorder=4,
                label=rf"crams {label} (MAP)",
            )

        # ── residual: (data - model) / model ──
        if med_on_band is not None and lo is not None:
            okb = med_on_band > 0
            axr.fill_between(
                R_band[okb], (lo / med_on_band - 1.0)[okb], (hi / med_on_band - 1.0)[okb],
                color=color, alpha=0.22, linewidth=0,
            )
        if len(x):
            m_at = _model_at(x, R_model, y_med)
            good = m_at > 0
            axr.errorbar(
                x[good], y[good] / m_at[good] - 1.0,
                yerr=[err_lo[good] / m_at[good], err_hi[good] / m_at[good]],
                fmt="o", color=color, markersize=4,
                elinewidth=1.2, capsize=0, zorder=3,
            )

    axr.axhline(0.0, color="k", linewidth=0.8, zorder=1)
    _shade_excluded((ax, axr), disp_min, r_min,
                    label=rf"$R < {r_min:g}$ GV (not fitted)")

    ax.set_xscale("log")
    ax.set_xlim(disp_min, BE_ISOTOPE_R_MAX)
    ax.set_ylabel(r"$R^{2.7} \times \Phi\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    axr.set_xlabel(r"$R$ [GV]")
    axr.set_ylabel("residual")
    ax.legend(fontsize=12)
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
    p.add_argument("--no-map", dest="show_map", action="store_false",
                   help="do not overlay the MAP (maximum-a-posteriori) model; "
                        "by default it is drawn dotted alongside the solid median")
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    prefix     = Path(args.output) if args.output else chain_path.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rng        = np.random.default_rng(args.seed)

    chain, param_names, acceptance, chain_frag, chain_halo, log_prob = _load_chain(
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

    # MAP (maximum-a-posteriori) model, overlaid dotted unless disabled/unavailable.
    spectra_map = None
    if args.show_map:
        if log_prob is None:
            print("No log_prob in chain file; skipping MAP overlay.")
        else:
            print("Running MAP (maximum-a-posteriori) model…")
            spectra_map = runner.run(_map_params(chain, log_prob, param_names))
            if spectra_map is None:
                print("MAP parameter run failed; drawing median only.")

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
        fig = _plot_panel(panel, spectra_median, spectra_samples, spectra_map)
        out = prefix.parent / f"{prefix.name}_{panel.tag}.pdf"
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"Saved: {out.resolve()}")

    extra_plots = [
        ("Be7_Be9", _plot_isotope_fluxes([BE7, BE9], spectra_median, spectra_samples, spectra_map)),
        ("Be10", _plot_isotope_fluxes([BE10], spectra_median, spectra_samples, spectra_map)),
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
