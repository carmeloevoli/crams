#!/usr/bin/env python3
"""Paper-specific MCMC plots for the PRD manuscript.

This script is intentionally separate from the general plotting utilities. It
can duplicate a little code from plot_mcmc_posterior.py and
plot_mcmc_params_posterior.py so the publication figures can evolve without
making the reusable scripts more specialized.

The first supported figure compares the six selected-scenario injection-slope
posterior density functions: three slopes for each of the Evoli et al. 2026 W93
and ST99 fragmentation models.

Usage
-----
    python plot_mcmc_prd.py
    python plot_mcmc_prd.py --scenario variable_h_beb
    python plot_mcmc_prd.py --scenario variable_h_preliminary_be
    python plot_mcmc_prd.py --scenario variable_h_variable_xsecs_preliminary_be
    python plot_mcmc_prd.py --output figs/mcmc_evoli2026_slopes.pdf
    python plot_mcmc_prd.py --violin-output figs/mcmc_evoli2026_slopes_violin.pdf
    python plot_mcmc_prd.py --d0h-delta-output figs/mcmc_evoli2026_d0h_delta.pdf
    python plot_mcmc_prd.py --halo-output figs/mcmc_evoli2026_halo_size.pdf
    python plot_mcmc_prd.py --be-nuisance-output figs/mcmc_evoli2026_Be7_Be9_nuisance.pdf
    python plot_mcmc_prd.py --be-isotopes-output figs/mcmc_evoli2026_Be7_Be9_fluxes.pdf
    python plot_mcmc_prd.py --h-he-output figs/mcmc_evoli2026_observables.pdf
    python plot_mcmc_prd.py --be10-be9-output figs/mcmc_evoli2026_Be10_Be9.pdf
"""
from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from fitting import KISS_DIR, _read_kiss_table
from runner import CramsRunner

MCMC_DIR = Path(__file__).resolve().parent
STYLE = MCMC_DIR.parent.parent / "crams-plots" / "crams.mplstyle"
PRELIMINARY_DIR = MCMC_DIR / "preliminary"
BE_ISOTOPES_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_isotopes_ECRS.csv"
# NOTE: the _ECRS ratios table is binned in Ek/n (column Ek_over_n_GeV_per_n),
# not rigidity, so it is incompatible with the rigidity-based Be10/Be9 figure.
# Keep the R_GV-binned table here until that figure is reworked for Ek/n.
BE_RATIOS_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_ratios.csv"
DEFAULT_CHAIN_DIR = MCMC_DIR / "mcmc_chains"
DEFAULT_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_slopes.pdf"
DEFAULT_VIOLIN_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_slopes_violin.pdf"
DEFAULT_D0H_DELTA_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_d0h_delta.pdf"
DEFAULT_HALO_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_halo_size.pdf"
DEFAULT_BE_NUISANCE_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_Be7_Be9_nuisance.pdf"
DEFAULT_BE_ISOTOPES_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_Be7_Be9_fluxes.pdf"
DEFAULT_H_HE_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_observables.pdf"
DEFAULT_BE10_BE9_OUTPUT = MCMC_DIR / "figs" / "mcmc_evoli2026_Be10_Be9.pdf"
VIOLIN_YLIM = (4.30, 4.43)
VIOLIN_BINS = 180
D0H_DELTA_XLIM = (0.30, 0.58)
D0H_DELTA_YLIM = (0.50, 0.60)
D0H_DELTA_BINS = 140
POWER = 2.7
BE_ISOTOPE_R_MIN = 5.0
BE_ISOTOPE_R_MAX = 45.0
PREDICTION_INTERVAL = ("95.45%", 2.2750131948, 97.7249868052)
SLOPE_STYLES = [
    ("hslope", r"$\gamma_\mathrm{H}$", "#0072B2"),
    ("heslope", r"$\gamma_\mathrm{He}$", "#D55E00"),
    ("slope", r"$\gamma_{Z\geq 3}$", "#8D5524"),
]

BE_NUISANCE_STYLES = [
    ("fudge_be7", r"$^{7}\mathrm{Be}$", "#009E73"),
    ("fudge_be9", r"$^{9}\mathrm{Be}$", "#CC79A7"),
    ("fudge_be10", r"$^{10}\mathrm{Be}$", "#E69F00"),
]

BE_ISOTOPE_STYLES = [
    ("Be7", "Be7_flux_R2p7", r"$^{7}\mathrm{Be}$", "#009E73", "o"),
    ("Be9", "Be9_flux_R2p7", r"$^{9}\mathrm{Be}$", "#CC79A7", "^"),
    ("Be10", "Be10_flux_R2p7", r"$^{10}\mathrm{Be}$", "#E69F00", "s"),
]

MODEL_STYLES = {
    "w93": {
        "fragmentation_model": "evoli2026w93",
        "label": "W93",
        "color": "#0072B2",
        "linestyle": "-",
        "marker": "o",
    },
    "st99": {
        "fragmentation_model": "evoli2026st99",
        "label": "ST99",
        "color": "#D55E00",
        "linestyle": "--",
        "marker": "^",
    },
}

SIGMA_INTERVALS = [
    ("95.45%", 2.2750131948, 97.7249868052, 4.0, 0.95),
]

CONTOUR_INTERVALS = [
    ("68.27%", 0.6826894921, "-"),
    ("95.45%", 0.9544997361, "--"),
]

SCENARIOS = (
    "baseline",
    "variable_h_beb",
    "variable_h_preliminary_be",
    "variable_h_variable_xsecs_preliminary_be",
)


@dataclass(frozen=True)
class ChainPosterior:
    path: Path
    model_key: str
    label: str
    linestyle: str
    chain: np.ndarray
    param_names: list[str]
    acceptance: float


@dataclass(frozen=True)
class ObservablePanel:
    numerator: str
    filename: str
    label: str
    r_min: float
    r_max: float
    denominator: str = ""

    @property
    def is_ratio(self) -> bool:
        return bool(self.denominator)

    @property
    def power(self) -> float:
        return 0.0 if self.is_ratio else POWER

    @property
    def y_label(self) -> str:
        if self.is_ratio:
            return self.label
        return r"$R^{2.7} I\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$"

    @property
    def tag(self) -> str:
        return f"{self.numerator}_{self.denominator}" if self.is_ratio else self.numerator


OBSERVABLE_PANELS = [
    ObservablePanel("H", "AMS-02_H_rigidity.txt", r"$\mathrm{H}$", 4.0, 1500.0),
    ObservablePanel("He", "AMS-02_He_rigidity.txt", r"$\mathrm{He}$", 4.0, 2500.0),
    ObservablePanel("C", "AMS-02_C_rigidity.txt", r"$\mathrm{C}$", 4.0, 2500.0),
    ObservablePanel("O", "AMS-02_O_rigidity.txt", r"$\mathrm{O}$", 4.0, 2500.0),
    ObservablePanel("B", "AMS-02_B_C_rigidity.txt", r"$\mathrm{B}/\mathrm{C}$",
                    4.0, 2500.0, denominator="C"),
    ObservablePanel("B", "AMS-02_B_O_rigidity.txt", r"$\mathrm{B}/\mathrm{O}$",
                    4.0, 2500.0, denominator="O"),
    ObservablePanel("C", "AMS-02_C_O_rigidity.txt", r"$\mathrm{C}/\mathrm{O}$",
                    4.0, 2500.0, denominator="O"),
    ObservablePanel("Be", "AMS-02_Be_B_rigidity.txt", r"$\mathrm{Be}/\mathrm{B}$",
                    4.0, 1200.0, denominator="B"),
    ObservablePanel("Be", "AMS-02_Be_C_rigidity.txt", r"$\mathrm{Be}/\mathrm{C}$",
                    4.0, 1200.0, denominator="C"),
]

BE10_BE9_PANEL = ObservablePanel(
    "Be10",
    "",
    r"$^{10}\mathrm{Be}/^{9}\mathrm{Be}$",
    5.0,
    100.0,
    denominator="Be9",
)


def _format_halo_tag(halosize: float) -> str:
    if float(halosize).is_integer():
        return str(int(halosize))
    return f"{halosize:g}".replace(".", "p")


def _scenario_chain_path(
    chain_dir: Path,
    model_key: str,
    scenario: str,
    halosize: float,
) -> Path:
    model = MODEL_STYLES[model_key]["fragmentation_model"]
    # run_all_mcmc.sh tags chains with the ECRS best-fit seed suffix (h is a
    # free parameter, so no halo tag is baked into the filename).
    del halosize  # kept for signature compatibility; unused in the ECRS naming
    return chain_dir / f"mcmc_{model}_{scenario}_ecrs.npz"


def _load_chain(
    path: Path,
    model_key: str,
    discard: int = 0,
    thin: int = 1,
) -> ChainPosterior:
    style = MODEL_STYLES[model_key]
    d = np.load(path, allow_pickle=True)
    chain = d["chain"]
    if discard:
        chain = chain[discard:]
    if thin > 1:
        chain = chain[::thin]
    return ChainPosterior(
        path=path,
        model_key=model_key,
        label=str(style["label"]),
        linestyle=str(style["linestyle"]),
        chain=chain,
        param_names=[str(name) for name in d["param_names"]],
        acceptance=float(d["acceptance"]),
    )


def _finite_param(chain: ChainPosterior, param_name: str) -> np.ndarray:
    if param_name not in chain.param_names:
        raise KeyError(f"{chain.path} has no parameter {param_name!r}")
    values = chain.chain[:, chain.param_names.index(param_name)]
    return values[np.isfinite(values)]


def _shared_bins(values: list[np.ndarray], nbins: int) -> np.ndarray:
    all_values = np.concatenate([v for v in values if len(v)])
    lo, hi = float(np.min(all_values)), float(np.max(all_values))
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise ValueError("cannot build histogram bins from non-finite values")
    if lo == hi:
        pad = max(abs(lo) * 1e-3, 1e-3)
    else:
        pad = 0.06 * (hi - lo)
    return np.linspace(lo - pad, hi + pad, nbins + 1)


def _summary(values: np.ndarray) -> tuple[float, float, float]:
    lo, median, hi = np.nanpercentile(values, [15.8655253931, 50, 84.1344746069])
    return float(lo), float(median), float(hi)


def _apply_style() -> None:
    use_tex = STYLE.exists() and shutil.which("latex") is not None
    if STYLE.exists():
        try:
            plt.style.use(str(STYLE))
            if not use_tex:
                plt.rcParams["text.usetex"] = False
        except Exception:
            pass


def _histogram_violin_density(
    values: np.ndarray,
    y_range: tuple[float, float],
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(y_range[0], y_range[1], nbins + 1)
    finite = values[np.isfinite(values)]
    finite = finite[(finite >= y_range[0]) & (finite <= y_range[1])]
    counts, _ = np.histogram(finite, bins=edges)
    density = counts.astype(float)

    # A short symmetric kernel keeps the violin smooth without invoking a KDE.
    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0])
    kernel /= kernel.sum()
    for _ in range(2):
        density = np.convolve(density, kernel, mode="same")

    if np.max(density) > 0.0:
        density /= np.max(density)

    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density


def _smooth_1d_histogram_density(
    values: np.ndarray,
    bins: np.ndarray,
    passes: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    finite = values[np.isfinite(values)]
    density, edges = np.histogram(finite, bins=bins, density=True)
    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0])
    kernel /= kernel.sum()
    for _ in range(passes):
        density = np.convolve(density, kernel, mode="same")
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density


def _draw_histogram_violin(
    ax,
    values: np.ndarray,
    position: float,
    width: float,
    color: str,
    y_range: tuple[float, float],
) -> None:
    y, density = _histogram_violin_density(values, y_range, VIOLIN_BINS)
    half_width = 0.5 * width * density
    ax.fill_betweenx(
        y,
        position - half_width,
        position + half_width,
        facecolor=color,
        edgecolor=color,
        linewidth=2.0,
        alpha=0.32,
    )


def _smooth_grid(grid: np.ndarray, passes: int = 2) -> np.ndarray:
    smoothed = grid.astype(float)
    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0])
    kernel /= kernel.sum()
    for _ in range(passes):
        smoothed = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"),
                                       axis=0, arr=smoothed)
        smoothed = np.apply_along_axis(lambda col: np.convolve(col, kernel, mode="same"),
                                       axis=1, arr=smoothed)
    return smoothed


def _credible_density_levels(
    density: np.ndarray,
    fractions: list[float],
) -> list[float]:
    flat = density.ravel()
    flat = flat[np.isfinite(flat) & (flat > 0)]
    if len(flat) < 2:
        return []

    ordered = np.sort(flat)[::-1]
    cdf = np.cumsum(ordered) / np.sum(ordered)
    levels = []
    for fraction in fractions:
        idx = min(np.searchsorted(cdf, fraction), len(ordered) - 1)
        levels.append(float(ordered[idx]))

    max_density = float(np.max(flat))
    return sorted({level for level in levels if 0.0 < level < max_density})


def _posterior_2d_density(
    x: np.ndarray,
    y: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    nbins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    finite = (
        np.isfinite(x) & np.isfinite(y)
        & (x >= x_range[0]) & (x <= x_range[1])
        & (y >= y_range[0]) & (y <= y_range[1])
    )
    counts, x_edges, y_edges = np.histogram2d(
        x[finite], y[finite], bins=nbins, range=[x_range, y_range])
    density = _smooth_grid(counts, passes=2)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    return x_centers, y_centers, density


def _median_fit_params(chain: ChainPosterior, halo_size: float) -> dict[str, float]:
    medians = np.nanmedian(chain.chain, axis=0)
    params = dict(zip(chain.param_names, medians))
    params.setdefault("h", halo_size)
    return params


def _sample_fit_params(
    chain: ChainPosterior,
    halo_size: float,
    nsamples: int,
    rng: np.random.Generator,
) -> list[dict[str, float]]:
    n = min(nsamples, len(chain.chain))
    idx = rng.choice(len(chain.chain), size=n, replace=False)
    samples = []
    for i in idx:
        params = dict(zip(chain.param_names, chain.chain[i]))
        params.setdefault("h", halo_size)
        samples.append(params)
    return samples


def _run_flux_predictions(
    chain: ChainPosterior,
    halo_size: float,
    nsamples: int,
    rng: np.random.Generator,
    build_dir: Path | None,
) -> tuple[dict[str, np.ndarray], list[dict[str, np.ndarray]]]:
    style = MODEL_STYLES[chain.model_key]
    runner = CramsRunner(
        build_dir=build_dir,
        read_isotopes=True,
        fragmentation_model=style["fragmentation_model"],
    )

    # The Be fudge factors are written to the .ini and applied inside crams
    # (runner.run), so no post-scaling of the returned spectra is needed.
    median_params = _median_fit_params(chain, halo_size)
    median_spectra = runner.run(median_params)
    if median_spectra is None:
        raise RuntimeError(f"median model run failed for {chain.label}")

    sample_spectra = []
    for params in _sample_fit_params(chain, halo_size, nsamples, rng):
        spectra = runner.run(params)
        if spectra is not None:
            sample_spectra.append(spectra)

    print(f"{chain.label}: {len(sample_spectra)}/{min(nsamples, len(chain.chain))} "
          "posterior flux samples succeeded")
    return median_spectra, sample_spectra


def _flux_data(
    filename: str,
    r_min: float,
    r_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x, y, err_lo, err_hi = _read_kiss_table(KISS_DIR / filename)
    cut = (
        np.isfinite(x) & np.isfinite(y) & np.isfinite(err_lo) & np.isfinite(err_hi)
        & (err_lo > 0.0) & (err_hi > 0.0)
        & (x >= r_min) & (x <= r_max)
    )
    return x[cut], y[cut], err_lo[cut], err_hi[cut]


def _preliminary_csv_data(
    path: Path,
    column: str,
    r_min: float,
    r_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    x = np.asarray(data["R_GV"], dtype=float)
    y = np.asarray(data[column], dtype=float)
    err_column = column[:-len("_flux_R2p7")] if column.endswith("_flux_R2p7") else column
    err_lo = np.asarray(data[f"{err_column}_err_minus"], dtype=float)
    err_hi = np.asarray(data[f"{err_column}_err_plus"], dtype=float)
    cut = (
        np.isfinite(x) & np.isfinite(y) & np.isfinite(err_lo) & np.isfinite(err_hi)
        & (err_lo > 0.0) & (err_hi > 0.0)
        & (x >= r_min) & (x <= r_max)
    )
    return x[cut], y[cut], err_lo[cut], err_hi[cut]


def _strip_known_output_tag(stem: str) -> str:
    tags = [
        "slopes_violin",
        "slopes",
        "d0h_delta",
        "halo_size",
        "Be7_Be9_nuisance",
        "Be7_Be9_fluxes",
        "Be10_Be9",
        "H_He_BC_BO",
        "observables",
        BE10_BE9_PANEL.tag,
        *(panel.tag for panel in OBSERVABLE_PANELS),
    ]
    for tag in tags:
        marker = f"_{tag}"
        if stem.endswith(marker):
            return stem[: -len(marker)]
    return stem


def _tagged_output_path(output: Path, scenario: str, tag: str) -> Path:
    suffix = output.suffix or ".pdf"
    stem = _strip_known_output_tag(output.stem)
    if scenario != "baseline":
        scenario_marker = f"_{scenario}"
        if not stem.endswith(scenario_marker):
            stem = f"{stem}{scenario_marker}"
    return output.with_name(f"{stem}_{tag}{suffix}")


def _interpolate_loglog(
    r_model: np.ndarray,
    y_model: np.ndarray,
    r_data: np.ndarray,
) -> np.ndarray:
    mask = (r_model > 0.0) & (y_model > 0.0) & np.isfinite(r_model) & np.isfinite(y_model)
    if mask.sum() < 2:
        return np.full_like(r_data, np.nan)
    log_y = np.interp(
        np.log(r_data),
        np.log(r_model[mask]),
        np.log(y_model[mask]),
        left=np.nan,
        right=np.nan,
    )
    return np.exp(log_y)


def _model_observable(
    spectra: dict[str, np.ndarray],
    panel: ObservablePanel,
) -> tuple[np.ndarray, np.ndarray]:
    r_model = spectra["R"]
    numerator = spectra[panel.numerator]
    if panel.is_ratio:
        denominator = spectra[panel.denominator]
        y_model = np.full_like(numerator, np.nan)
        np.divide(numerator, denominator, out=y_model, where=denominator > 0.0)
    else:
        y_model = numerator
    return r_model, y_model


def _median_residuals(
    spectra: dict[str, np.ndarray],
    panel: ObservablePanel,
    r_data: np.ndarray,
    y_data: np.ndarray,
    err_lo: np.ndarray,
    err_hi: np.ndarray,
) -> np.ndarray:
    r_model, observable = _model_observable(spectra, panel)
    y_model = _interpolate_loglog(r_model, observable, r_data)
    residual = y_model - y_data
    sigma = np.where(residual >= 0.0, err_hi, err_lo)
    out = np.full_like(residual, np.nan)
    np.divide(residual, sigma, out=out, where=sigma > 0.0)
    return out


def plot_slope_comparison(
    chains: list[ChainPosterior],
    output: Path,
    nbins: int,
) -> None:
    """Plot all selected W93/ST99 injection-slope posteriors on one axes."""
    _apply_style()
    fig, ax = plt.subplots(figsize=(12, 10))
    all_values = [
        _finite_param(chain, param_name)
        for param_name, _, _ in SLOPE_STYLES
        for chain in chains
    ]
    bins = _shared_bins(all_values, nbins)

    for param_name, param_label, color in SLOPE_STYLES:
        for chain in chains:
            values = _finite_param(chain, param_name)
            ax.hist(
                values,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=3.0,
                color=color,
                linestyle=chain.linestyle,
            )
            lo, median, hi = _summary(values)
            print(
                f"{chain.label:4s} {param_name:7s}: "
                f"{median:.4f} +{hi - median:.4f} / -{median - lo:.4f}"
            )

    slope_handles = [
        Line2D([0], [0], color=color, lw=3.0, label=label)
        for _, label, color in SLOPE_STYLES
    ]
    model_handles = [
        Line2D([0], [0], color="0.2", lw=3.0,
               linestyle=MODEL_STYLES[key]["linestyle"],
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    slope_legend = ax.legend(handles=slope_handles, frameon=False, loc="upper left",
                             fontsize=29)
    ax.add_artist(slope_legend)
    ax.legend(handles=model_handles, frameon=False, loc="upper right", fontsize=25)

    ax.set_xlim(4.31, 4.43)
    ax.set_xlabel("Injection spectral index")
    ax.set_ylabel("Posterior density")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_d0h_delta_comparison(
    chains: list[ChainPosterior],
    output: Path,
) -> None:
    """Plot overlaid selected W93/ST99 contours in D0/H versus delta."""
    _apply_style()
    fig, ax = plt.subplots(figsize=(11, 9.5))
    fractions = [fraction for _, fraction, _ in CONTOUR_INTERVALS]

    for chain in chains:
        d0_h = _finite_param(chain, "d0_h")
        delta = _finite_param(chain, "delta")
        n = min(len(d0_h), len(delta))
        d0_h, delta = d0_h[:n], delta[:n]
        x_centers, y_centers, density = _posterior_2d_density(
            d0_h, delta, D0H_DELTA_XLIM, D0H_DELTA_YLIM, D0H_DELTA_BINS)
        levels = _credible_density_levels(density, fractions)
        if not levels:
            print(f"Skipping D0/H-delta contours for {chain.label}: no finite density")
            continue

        style = MODEL_STYLES[chain.model_key]
        # _credible_density_levels returns increasing density thresholds; map
        # them back to the requested outer-to-inner visual styles.
        style_by_level = {
            level: linestyle
            for level, (_, _, linestyle) in zip(levels, reversed(CONTOUR_INTERVALS))
        }
        for level in levels:
            ax.contour(
                x_centers,
                y_centers,
                density.T,
                levels=[level],
                colors=[style["color"]],
                linestyles=[style_by_level[level]],
                linewidths=3.0,
            )

        d0_med = float(np.nanmedian(d0_h))
        delta_med = float(np.nanmedian(delta))
        ax.plot(d0_med, delta_med, marker="o", markersize=9.0,
                color=style["color"], markeredgecolor="white", markeredgewidth=1.2)
        print(f"{chain.label:4s} d0_h/delta median: {d0_med:.4f}, {delta_med:.4f}")

    model_handles = [
        Line2D([0], [0], color=MODEL_STYLES[key]["color"], lw=4.0,
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    interval_handles = [
        Line2D([0], [0], color="0.25", lw=3.0, linestyle=linestyle, label=label)
        for label, _, linestyle in CONTOUR_INTERVALS
    ]
    model_legend = ax.legend(handles=model_handles, frameon=False,
                             loc="upper right", fontsize=25)
    ax.add_artist(model_legend)
    ax.legend(handles=interval_handles, frameon=False, loc="lower left",
              fontsize=22)

    ax.set_xlim(*D0H_DELTA_XLIM)
    ax.set_ylim(*D0H_DELTA_YLIM)
    ax.set_xlabel(r"$D_0/H$ [$10^{28}\,\mathrm{cm^2\,s^{-1}\,kpc^{-1}}$]")
    ax.set_ylabel(r"$\delta$")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_halo_size_posterior(
    chains: list[ChainPosterior],
    output: Path,
    nbins: int,
) -> None:
    """Plot the sampled halo-size posterior for W93/ST99, if available."""
    values_by_chain = []
    for chain in chains:
        if "h" not in chain.param_names:
            continue
        values = _finite_param(chain, "h")
        if len(values):
            values_by_chain.append((chain, values))

    if not values_by_chain:
        print("Skipping halo-size posterior: H is fixed in this scenario")
        return

    _apply_style()
    fig, ax = plt.subplots(figsize=(11, 8.5))
    bins = _shared_bins([values for _, values in values_by_chain], nbins)
    summary_lines = []

    for chain, values in values_by_chain:
        style = MODEL_STYLES[chain.model_key]
        x_density, y_density = _smooth_1d_histogram_density(values, bins)
        ax.fill_between(
            x_density,
            0.0,
            y_density,
            color=style["color"],
            alpha=0.1,
            linewidth=0,
        )
        ax.plot(
            x_density,
            y_density,
            linewidth=4.0,
            color=style["color"],
            linestyle=style["linestyle"],
        )
        lower95, median, upper95 = np.nanpercentile(
            values, [2.2750131948, 50.0, 97.7249868052])
        ax.axvline(median, color=style["color"], linestyle=":", linewidth=4.0)
        summary_lines.append((
            style["color"],
            rf"{chain.label}: $H={median:.2f}^{{+{upper95 - median:.2f}}}"
            rf"_{{-{median - lower95:.2f}}}\,\mathrm{{kpc}}$",
        ))
        print(
            f"{chain.label:4s} H: {median:.2f} kpc "
            f"[{lower95:.2f}, {upper95:.2f}] 95.45%"
        )

    model_handles = [
        Line2D([0], [0], color=MODEL_STYLES[key]["color"], lw=4.0,
               linestyle=MODEL_STYLES[key]["linestyle"],
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    median_handle = Line2D([0], [0], color="0.25", lw=4.0, linestyle=":",
                           label="median")
    ax.legend(handles=[*model_handles, median_handle], frameon=False,
              loc="upper right", fontsize=20)
    for i, (color, text) in enumerate(summary_lines):
        ax.text(
            0.04,
            0.95 - 0.075 * i,
            text,
            transform=ax.transAxes,
            color=color,
            fontsize=20,
            ha="left",
            va="top",
        )
    ax.set_xlim(6.0, 16.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$H$ [kpc]")
    ax.set_ylabel("Posterior density")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_be_nuisance_posterior(
    chains: list[ChainPosterior],
    output: Path,
    nbins: int,
) -> None:
    """Plot the Be7 and Be9 nuisance-factor posteriors, if available."""
    entries = []
    for param_name, isotope_label, color in BE_NUISANCE_STYLES:
        for chain in chains:
            if param_name not in chain.param_names:
                continue
            values = _finite_param(chain, param_name)
            if len(values):
                entries.append((chain, param_name, isotope_label, color, values))

    if not entries:
        print("Skipping Be7/Be9 nuisance posterior: nuisance factors are fixed")
        return

    _apply_style()
    fig, ax = plt.subplots(figsize=(11, 8.5))
    bins = _shared_bins([values for *_, values in entries], nbins)

    for chain, param_name, isotope_label, color, values in entries:
        x_density, y_density = _smooth_1d_histogram_density(values, bins)
        ax.fill_between(
            x_density,
            0.0,
            y_density,
            color=color,
            alpha=0.08,
            linewidth=0,
        )
        ax.plot(
            x_density,
            y_density,
            linewidth=3.5,
            color=color,
            linestyle=chain.linestyle,
        )
        lower95, median, upper95 = np.nanpercentile(
            values, [2.2750131948, 50.0, 97.7249868052])
        print(
            f"{chain.label:4s} {param_name:10s}: {median:.3f} "
            f"[{lower95:.3f}, {upper95:.3f}] 95.45%"
        )

    ax.axvline(1.0, color="0.35", linestyle=":", linewidth=2.0)

    isotope_handles = [
        Line2D([0], [0], color=color, lw=4.0, label=label)
        for _, label, color in BE_NUISANCE_STYLES
    ]
    model_handles = [
        Line2D([0], [0], color="0.25", lw=3.5,
               linestyle=MODEL_STYLES[key]["linestyle"],
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    nominal_handle = Line2D([0], [0], color="0.35", lw=2.0,
                            linestyle=":", label="nominal")
    isotope_legend = ax.legend(handles=isotope_handles, frameon=False,
                               loc="upper left", fontsize=22)
    ax.add_artist(isotope_legend)
    ax.legend(handles=[*model_handles, nominal_handle], frameon=False,
              loc="upper right", fontsize=20)

    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Beryllium isotope nuisance factor")
    ax.set_ylabel("Posterior density")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_be7_be9_fluxes(
    chains: list[ChainPosterior],
    predictions: dict[str, tuple[dict[str, np.ndarray], list[dict[str, np.ndarray]]]],
    output: Path,
) -> None:
    """Plot preliminary Be7 and Be9 fluxes with W93/ST99 posterior predictions."""
    _apply_style()
    interval_label, lower_pct, upper_pct = PREDICTION_INTERVAL

    fig = plt.figure(figsize=(11.0, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.0), hspace=0.05)
    ax = fig.add_subplot(grid[0, 0])
    ax_res = fig.add_subplot(grid[1, 0], sharex=ax)
    plotted_model = False

    for isotope, csv_column, isotope_label, color, data_marker in BE_ISOTOPE_STYLES:
        x, y, err_lo, err_hi = _preliminary_csv_data(
            BE_ISOTOPES_CSV, csv_column, BE_ISOTOPE_R_MIN, BE_ISOTOPE_R_MAX)
        ax.errorbar(
            x,
            y,
            yerr=[err_lo, err_hi],
            fmt=data_marker,
            color=color,
            markersize=4.5,
            elinewidth=1.2,
            capsize=0,
            zorder=5,
        )

        for chain in chains:
            style = MODEL_STYLES[chain.model_key]
            median_spectra, sample_spectra = predictions[chain.model_key]
            if isotope not in median_spectra:
                print(f"Skipping {isotope} flux for {chain.label}: isotope spectra missing")
                continue

            valid_samples = [s for s in sample_spectra if isotope in s]
            if valid_samples:
                r_model = valid_samples[0]["R"]
                mask = (r_model >= BE_ISOTOPE_R_MIN) & (r_model <= BE_ISOTOPE_R_MAX)
                band = np.array([
                    s["R"][mask] ** POWER * s[isotope][mask]
                    for s in valid_samples
                ])
                lo, hi = np.nanpercentile(band, [lower_pct, upper_pct], axis=0)
                ax.fill_between(
                    r_model[mask],
                    lo,
                    hi,
                    color=color,
                    alpha=0.10,
                    linewidth=0,
                    zorder=1,
                )

            r_med = median_spectra["R"]
            mask = (r_med >= BE_ISOTOPE_R_MIN) & (r_med <= BE_ISOTOPE_R_MAX)
            y_med = r_med[mask] ** POWER * median_spectra[isotope][mask]
            ax.plot(
                r_med[mask],
                y_med,
                color=color,
                linestyle=style["linestyle"],
                linewidth=3.0,
                zorder=6,
            )
            plotted_model = True

            y_interp = _interpolate_loglog(
                r_med[mask],
                y_med,
                x,
            )
            residual = y_interp - y
            sigma = np.where(residual >= 0.0, err_hi, err_lo)
            residual_sigma = np.full_like(residual, np.nan)
            np.divide(residual, sigma, out=residual_sigma, where=sigma > 0.0)
            ax_res.plot(
                x,
                residual_sigma,
                color=color,
                linestyle="None",
                marker=style["marker"],
                markersize=4.0,
            )

    if not plotted_model:
        print("Skipping Be7/Be9 flux plot: isotope spectra missing")
        plt.close(fig)
        return

    ax.set_xscale("log")
    ax_res.set_xscale("log")
    ax.set_xlim(BE_ISOTOPE_R_MIN, BE_ISOTOPE_R_MAX)
    ax.set_ylabel(r"$R^{2.7} I\ [\mathrm{GV^{1.7}\ m^{-2}\ s^{-1}\ sr^{-1}}]$")
    ax.tick_params(labelbottom=False)
    ax_res.axhline(0.0, color="0.35", linestyle=":", linewidth=1.3)
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.set_xlabel(r"$R$ [GV]")
    ax_res.set_ylabel(r"$\Delta/\sigma$")

    isotope_handles = [
        Line2D([0], [0], color=color, marker=marker, lw=3.0,
               label=isotope_label)
        for _, _, isotope_label, color, marker in BE_ISOTOPE_STYLES
    ]
    model_handles = [
        Line2D([0], [0], color="0.25", lw=3.0,
               linestyle=MODEL_STYLES[key]["linestyle"],
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    data_handle = Line2D([0], [0], color="0.25", marker="o", linestyle="None",
                         label="AMS-02 preliminary")
    band_handle = Line2D([0], [0], color="0.25", lw=7.0, alpha=0.18,
                         label=interval_label)
    isotope_legend = ax.legend(handles=isotope_handles, frameon=False,
                               loc="upper left", fontsize=19)
    ax.add_artist(isotope_legend)
    ax.legend(handles=[*model_handles, data_handle, band_handle], frameon=False,
              loc="upper right", fontsize=17)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_h_he_flux_posteriors(
    chains: list[ChainPosterior],
    output: Path,
    scenario: str,
    halo_size: float,
    nsamples: int,
    rng: np.random.Generator,
    build_dir: Path | None,
) -> dict[str, tuple[dict[str, np.ndarray], list[dict[str, np.ndarray]]]]:
    """Plot posterior bands plus median residuals for key observables."""
    _apply_style()
    predictions = {
        chain.model_key: _run_flux_predictions(chain, halo_size, nsamples, rng, build_dir)
        for chain in chains
    }

    interval_label, lower_pct, upper_pct = PREDICTION_INTERVAL

    for panel in OBSERVABLE_PANELS:
        fig = plt.figure(figsize=(11.0, 8.5), constrained_layout=True)
        grid = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.0), hspace=0.05)
        ax = fig.add_subplot(grid[0, 0])
        ax_res = fig.add_subplot(grid[1, 0], sharex=ax)
        x, y, err_lo, err_hi = _flux_data(panel.filename, panel.r_min, panel.r_max)
        scale = x**panel.power

        ax.errorbar(
            x,
            scale * y,
            yerr=[scale * err_lo, scale * err_hi],
            fmt="o",
            color="0.18",
            markersize=4.0,
            elinewidth=1.2,
            capsize=0,
            label="AMS-02",
            zorder=5,
        )

        for chain in chains:
            style = MODEL_STYLES[chain.model_key]
            median_spectra, sample_spectra = predictions[chain.model_key]

            if sample_spectra:
                r_model, _ = _model_observable(sample_spectra[0], panel)
                mask = (r_model >= panel.r_min) & (r_model <= panel.r_max)
                band = np.array([
                    _model_observable(s, panel)[1][mask]
                    for s in sample_spectra
                ])
                lo, hi = np.nanpercentile(band, [lower_pct, upper_pct], axis=0)
                ax.fill_between(
                    r_model[mask],
                    r_model[mask] ** panel.power * lo,
                    r_model[mask] ** panel.power * hi,
                    color=style["color"],
                    alpha=0.18,
                    linewidth=0,
                    label=f"{style['label']} {interval_label}",
                )

            r_med, median_observable = _model_observable(median_spectra, panel)
            mask = (r_med >= panel.r_min) & (r_med <= panel.r_max)
            ax.plot(
                r_med[mask],
                r_med[mask] ** panel.power * median_observable[mask],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=3.0,
                label=f"{style['label']} median",
                zorder=6,
            )

            residual = _median_residuals(
                median_spectra, panel, x, y, err_lo, err_hi)
            ax_res.plot(
                x,
                residual,
                color=style["color"],
                linestyle="None",
                marker=style["marker"],
                markersize=4.0,
            )

        ax.set_xscale("log")
        ax_res.set_xscale("log")
        ax.set_xlim(panel.r_min, panel.r_max)
        ax.set_ylabel(panel.y_label)
        ax.tick_params(labelbottom=False)
        ax_res.axhline(0.0, color="0.35", linestyle=":", linewidth=1.3)
        ax_res.set_ylim(-5.0, 5.0)
        ax_res.set_xlabel(r"$R$ [GV]")
        ax_res.set_ylabel(r"$\Delta/\sigma$")

        ax.legend(frameon=False, fontsize=15, loc="best")
        out = _tagged_output_path(output, scenario, panel.tag)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight", dpi=300)
        plt.close(fig)
        print(f"Saved: {out.resolve()}")

    return predictions


def plot_be10_be9_posterior(
    chains: list[ChainPosterior],
    predictions: dict[str, tuple[dict[str, np.ndarray], list[dict[str, np.ndarray]]]],
    output: Path,
) -> None:
    """Plot selected-scenario predictions against preliminary AMS-02 Be10/Be9."""
    _apply_style()
    panel = BE10_BE9_PANEL
    interval_label, lower_pct, upper_pct = PREDICTION_INTERVAL
    x, y, err_lo, err_hi = _preliminary_csv_data(
        BE_RATIOS_CSV, "Be10_over_Be9", panel.r_min, panel.r_max)

    fig = plt.figure(figsize=(11.0, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.0), hspace=0.05)
    ax = fig.add_subplot(grid[0, 0])
    ax_res = fig.add_subplot(grid[1, 0], sharex=ax)

    ax.errorbar(
        x,
        y,
        yerr=[err_lo, err_hi],
        fmt="o",
        color="0.18",
        markersize=4.0,
        elinewidth=1.2,
        capsize=0,
        label="AMS-02 preliminary",
        zorder=5,
    )

    for chain in chains:
        style = MODEL_STYLES[chain.model_key]
        median_spectra, sample_spectra = predictions[chain.model_key]
        if panel.numerator not in median_spectra or panel.denominator not in median_spectra:
            print(f"Skipping {panel.label} for {chain.label}: isotope spectra missing")
            continue

        valid_samples = [
            s for s in sample_spectra
            if panel.numerator in s and panel.denominator in s
        ]
        if valid_samples:
            r_model, _ = _model_observable(valid_samples[0], panel)
            mask = (r_model >= panel.r_min) & (r_model <= panel.r_max)
            band = np.array([
                _model_observable(s, panel)[1][mask]
                for s in valid_samples
            ])
            lo, hi = np.nanpercentile(band, [lower_pct, upper_pct], axis=0)
            ax.fill_between(
                r_model[mask],
                lo,
                hi,
                color=style["color"],
                alpha=0.18,
                linewidth=0,
                label=f"{style['label']} {interval_label}",
            )

        r_med, median_observable = _model_observable(median_spectra, panel)
        mask = (r_med >= panel.r_min) & (r_med <= panel.r_max)
        ax.plot(
            r_med[mask],
            median_observable[mask],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=3.0,
            label=f"{style['label']} median",
            zorder=6,
        )

        residual = _median_residuals(median_spectra, panel, x, y, err_lo, err_hi)
        ax_res.plot(
            x,
            residual,
            color=style["color"],
            linestyle="None",
            marker=style["marker"],
            markersize=4.0,
        )

    ax.set_xscale("log")
    ax_res.set_xscale("log")
    ax.set_xlim(panel.r_min, panel.r_max)
    ax.set_ylabel(panel.label)
    ax.tick_params(labelbottom=False)
    ax_res.axhline(0.0, color="0.35", linestyle=":", linewidth=1.3)
    ax_res.set_ylim(-5.0, 5.0)
    ax_res.set_xlabel(r"$R$ [GV]")
    ax_res.set_ylabel(r"$\Delta/\sigma$")
    ax.legend(frameon=False, fontsize=15, loc="best")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def plot_slope_violin_comparison(
    chains: list[ChainPosterior],
    output: Path,
) -> None:
    """Plot grouped W93/ST99 violins for the three selected injection slopes."""
    _apply_style()
    fig, ax = plt.subplots(figsize=(11, 9.5))
    group_positions = np.arange(1, len(SLOPE_STYLES) + 1, dtype=float)
    offsets = {"w93": -0.17, "st99": 0.17}
    width = 0.28

    for i, (param_name, _, _) in enumerate(SLOPE_STYLES):
        x0 = group_positions[i]
        for chain in chains:
            values = _finite_param(chain, param_name)
            position = x0 + offsets[chain.model_key]
            style = MODEL_STYLES[chain.model_key]
            _draw_histogram_violin(ax, values, position, width, style["color"],
                                   VIOLIN_YLIM)

            for _, lower_pct, upper_pct, linewidth, alpha in SIGMA_INTERVALS:
                lo, hi = np.nanpercentile(values, [lower_pct, upper_pct])
                ax.vlines(position, lo, hi, color=style["color"],
                          linewidth=linewidth, alpha=alpha)

            _, median, _ = _summary(values)
            ax.plot(
                [position - 0.11, position + 0.11],
                [median, median],
                color=style["color"],
                linewidth=5.0,
                linestyle=":",
                solid_capstyle="round",
            )

    ax.set_xticks(group_positions)
    ax.set_xticklabels([label for _, label, _ in SLOPE_STYLES])
    ax.tick_params(axis="x", labelsize=35)
    ax.set_xlim(0.45, len(SLOPE_STYLES) + 0.55)
    ax.set_ylim(*VIOLIN_YLIM)
    ax.set_ylabel("Injection spectral index")

    model_handles = [
        Line2D([0], [0], color=MODEL_STYLES[key]["color"], lw=5.0,
               label=MODEL_STYLES[key]["label"])
        for key in ("w93", "st99")
    ]
    interval_handles = [
        Line2D([0], [0], color="0.25", lw=linewidth, alpha=alpha, label=label)
        for label, _, _, linewidth, alpha in SIGMA_INTERVALS
    ]
    interval_handles.append(Line2D([0], [0], color="0.25", lw=5.0,
                                   linestyle=":", label="median"))

    model_legend = ax.legend(handles=model_handles, frameon=False,
                             loc="upper right", fontsize=25)
    ax.add_artist(model_legend)
    ax.legend(handles=interval_handles, frameon=False, loc="lower right",
              fontsize=20)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {output.resolve()}")


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description="Create PRD-specific MCMC figures")
    p.add_argument("--scenario", choices=SCENARIOS, default="baseline",
                   help="MCMC scenario to plot")
    p.add_argument("--halosize", type=float, default=7.0,
                   help="halo tag used in default chain names")
    p.add_argument("--chain-dir", type=Path, default=DEFAULT_CHAIN_DIR,
                   help="directory containing the MCMC .npz chains")
    p.add_argument("--w93-chain", type=Path, default=None,
                   help="override the selected Evoli 2026 W93 chain path")
    p.add_argument("--st99-chain", type=Path, default=None,
                   help="override the selected Evoli 2026 ST99 chain path")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                   help="output PDF for the slope-comparison figure")
    p.add_argument("--violin-output", type=Path, default=DEFAULT_VIOLIN_OUTPUT,
                   help="output PDF for the slope violin-comparison figure")
    p.add_argument("--d0h-delta-output", type=Path, default=DEFAULT_D0H_DELTA_OUTPUT,
                   help="output PDF for the D0/H versus delta contour figure")
    p.add_argument("--halo-output", type=Path, default=DEFAULT_HALO_OUTPUT,
                   help="output PDF for the halo-size posterior figure")
    p.add_argument("--be-nuisance-output", type=Path, default=DEFAULT_BE_NUISANCE_OUTPUT,
                   help="output PDF for the Be7/Be9 nuisance-factor posterior figure")
    p.add_argument("--be-isotopes-output", type=Path, default=DEFAULT_BE_ISOTOPES_OUTPUT,
                   help="output PDF for the preliminary Be7/Be9 flux comparison figure")
    p.add_argument("--h-he-output", type=Path, default=DEFAULT_H_HE_OUTPUT,
                   help="output filename stem for the separate H, He, B/C, "
                        "B/O, C/O, and Be/B posterior PDFs")
    p.add_argument("--be10-be9-output", type=Path, default=DEFAULT_BE10_BE9_OUTPUT,
                   help="output PDF for predictions against preliminary Be10/Be9")
    p.add_argument("--flux-samples", type=int, default=100,
                   help="posterior samples per model for observable prediction bands")
    p.add_argument("--build-dir", type=Path, default=None,
                   help="path to crams build/ for model prediction figures")
    p.add_argument("--seed", type=int, default=0,
                   help="random seed for posterior prediction samples")
    p.add_argument("--discard", type=int, default=0,
                   help="discard first N flattened samples before plotting")
    p.add_argument("--thin", type=int, default=1,
                   help="thin the flattened chain by this factor before plotting")
    p.add_argument("--bins", type=int, default=70,
                   help="number of histogram bins")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)
    rng = np.random.default_rng(args.seed)
    scenario = args.scenario
    chain_paths = {
        "w93": args.w93_chain or _scenario_chain_path(
            args.chain_dir, "w93", scenario, args.halosize),
        "st99": args.st99_chain or _scenario_chain_path(
            args.chain_dir, "st99", scenario, args.halosize),
    }

    print(f"Scenario: {scenario}")
    chains = []
    for model_key in ("w93", "st99"):
        chain = _load_chain(chain_paths[model_key], model_key, args.discard, args.thin)
        chains.append(chain)
        print(
            f"{chain.label}: {chain.path}  "
            f"samples={chain.chain.shape[0]:,}  acceptance={chain.acceptance:.3f}"
        )

    plot_slope_comparison(
        chains, _tagged_output_path(args.output, scenario, "slopes"), args.bins)
    plot_slope_violin_comparison(
        chains, _tagged_output_path(args.violin_output, scenario, "slopes_violin"))
    plot_d0h_delta_comparison(
        chains, _tagged_output_path(args.d0h_delta_output, scenario, "d0h_delta"))
    plot_halo_size_posterior(
        chains, _tagged_output_path(args.halo_output, scenario, "halo_size"), args.bins)
    plot_be_nuisance_posterior(
        chains,
        _tagged_output_path(args.be_nuisance_output, scenario, "Be7_Be9_nuisance"),
        args.bins,
    )
    predictions = plot_h_he_flux_posteriors(
        chains,
        args.h_he_output,
        scenario,
        args.halosize,
        args.flux_samples,
        rng,
        args.build_dir,
    )
    plot_be7_be9_fluxes(
        chains,
        predictions,
        _tagged_output_path(args.be_isotopes_output, scenario, "Be7_Be9_fluxes"),
    )
    plot_be10_be9_posterior(
        chains,
        predictions,
        _tagged_output_path(args.be10_be9_output, scenario, BE10_BE9_PANEL.tag),
    )


if __name__ == "__main__":
    main()
