#!/usr/bin/env python3
"""Plot one crams model against Be/B and preliminary AMS-02 Be isotope data.

Parameters come from (in increasing precedence): the run_mcmc defaults, a crams
.ini (--ini), a chain median (--chain), and finally any key=value overrides on
the command line.

Usage
-----
    python plot_model_Be.py --ini bestfits/bestfit_evoli2026w93_h5.ini
    python plot_model_Be.py --ini bestfit.ini d0=2.5 h=5
    python plot_model_Be.py --chain crams_chain.npz --output figs/model
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
from run_mcmc import PARAMETERS
from runner import FRAGMENTATION_MODELS, CramsRunner

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"
PRELIMINARY_DIR = Path(__file__).parent / "preliminary"
BE_ISOTOPES_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_isotopes.csv"
BE_RATIOS_CSV = PRELIMINARY_DIR / "AMS-02_preliminary_Be_ratios.csv"

POWER = 2.7
DATA_COLOR = "tab:blue"
MODEL_COLOR = "tab:red"


@dataclass(frozen=True)
class Panel:
    numerator: str
    label: str
    tag: str
    denominator: str = ""
    filename: str = ""
    csv_path: Path | None = None
    csv_column: str = ""
    power: float = 0.0
    R_min: float | None = None
    R_max: float | None = None


PANELS = [
    Panel("Be", r"$\mathrm{Be}/\mathrm{B}$", "Be_B", denominator="B",
          filename="AMS-02_Be_B_rigidity.txt", R_min=2.0, R_max=2500.0),
    Panel("Be7", r"$R^{2.7}\,I(^{7}\mathrm{Be})$", "Be7",
          csv_path=BE_ISOTOPES_CSV, csv_column="Be7_flux_R2p7", power=POWER),
    Panel("Be9", r"$R^{2.7}\,I(^{9}\mathrm{Be})$", "Be9",
          csv_path=BE_ISOTOPES_CSV, csv_column="Be9_flux_R2p7", power=POWER),
    Panel("Be10", r"$R^{2.7}\,I(^{10}\mathrm{Be})$", "Be10",
          csv_path=BE_ISOTOPES_CSV, csv_column="Be10_flux_R2p7", power=POWER),
    Panel("Be9", r"$^{9}\mathrm{Be}/^{7}\mathrm{Be}$", "Be9_over_Be7",
          denominator="Be7", csv_path=BE_RATIOS_CSV, csv_column="Be9_over_Be7"),
    Panel("Be10", r"$^{10}\mathrm{Be}/^{9}\mathrm{Be}$", "Be10_over_Be9",
          denominator="Be9", csv_path=BE_RATIOS_CSV, csv_column="Be10_over_Be9"),
]


def read_ini(path: Path) -> dict[str, float]:
    """Read numeric 'key value' lines from a crams .ini (skips comments/strings)."""
    vals: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            vals[parts[0]] = float(parts[1])
        except ValueError:
            continue
    return vals


def read_ini_str(path: Path, key: str) -> str | None:
    """Return the value of a non-numeric 'key value' line."""
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0] == key:
            return parts[1]
    return None


def median_params(chain_file: Path) -> dict[str, float]:
    d = np.load(chain_file, allow_pickle=True)
    names = [str(x) for x in d["param_names"]]
    return dict(zip(names, np.median(d["chain"], axis=0)))


_CONJUGATE = {"d0": "d0_h", "d0_h": "d0", "rb": "rb_log", "rb_log": "rb"}


def _merge(dst: dict[str, float], src: dict[str, float]) -> None:
    """Update dst with src; a newly-set key removes its stale conjugate."""
    for key, value in src.items():
        dst.pop(_CONJUGATE.get(key, ""), None)
        dst[key] = value


def resolve_params(args) -> dict[str, float]:
    ini = {p.name: p.value for p in PARAMETERS}
    if args.ini:
        _merge(ini, read_ini(Path(args.ini)))
    if args.chain:
        _merge(ini, median_params(Path(args.chain)))
    for ov in args.overrides:
        if "=" not in ov:
            sys.exit(f"bad override '{ov}', expected key=value")
        key, value = ov.split("=", 1)
        _merge(ini, {key.strip(): float(value)})
    return ini


def _ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    y = np.full_like(num, np.nan)
    np.divide(num, den, out=y, where=den > 0)
    return y


def model_observable(spectra: dict[str, np.ndarray], panel: Panel) -> np.ndarray:
    num = spectra[panel.numerator]
    if panel.denominator:
        return _ratio(num, spectra[panel.denominator])
    return num


def _interpolate_loglog(
    R_model: np.ndarray,
    y_model: np.ndarray,
    R_data: np.ndarray,
) -> np.ndarray:
    """Log-log linear interpolation; returns NaN outside model range."""
    mask = (y_model > 0) & np.isfinite(y_model)
    if mask.sum() < 2:
        return np.full_like(R_data, np.nan)
    log_y = np.interp(
        np.log(R_data),
        np.log(R_model[mask]),
        np.log(y_model[mask]),
        left=np.nan,
        right=np.nan,
    )
    return np.exp(log_y)


def read_panel_data(panel: Panel) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (R, y, err_lo, err_hi) for one Be panel."""
    if panel.filename:
        return _read_kiss_table(KISS_DIR / panel.filename)

    if panel.csv_path is None:
        raise ValueError(f"panel {panel.tag} has no data source")

    data = np.genfromtxt(panel.csv_path, delimiter=",", names=True)
    column = panel.csv_column
    err_column = column[:-len("_flux_R2p7")] if column.endswith("_flux_R2p7") else column
    return (
        np.asarray(data["R_GV"], dtype=float),
        np.asarray(data[column], dtype=float),
        np.asarray(data[f"{err_column}_err_minus"], dtype=float),
        np.asarray(data[f"{err_column}_err_plus"], dtype=float),
    )


def plot_panel(ax, panel: Panel, spectra: dict[str, np.ndarray]) -> tuple[float, int]:
    x, y, lo, hi = read_panel_data(panel)
    data_mask = np.isfinite(x) & np.isfinite(y) & (lo > 0) & (hi > 0)
    x, y, lo, hi = x[data_mask], y[data_mask], lo[data_mask], hi[data_mask]

    ax.errorbar(
        x, y, yerr=[lo, hi],
        fmt="o", ms=4, color=DATA_COLOR,
        elinewidth=1.0, capsize=0, zorder=3,
        label="AMS-02 preliminary" if panel.csv_column else "AMS-02",
    )

    R = spectra["R"]
    ym = R**panel.power * model_observable(spectra, panel)
    y_at_data = _interpolate_loglog(R, ym, x)
    valid = np.isfinite(y_at_data)
    residual = y_at_data[valid] - y[valid]
    sigma = np.where(residual > 0, hi[valid], lo[valid])
    chi2 = float(np.sum((residual / sigma) ** 2)) if valid.any() else float("nan")
    n = int(np.sum(valid))

    R_min = panel.R_min if panel.R_min is not None else x.min() * 0.8
    R_max = panel.R_max if panel.R_max is not None else x.max() * 1.2
    model_mask = (R >= R_min) & (R <= R_max) & np.isfinite(ym)
    ax.plot(R[model_mask], ym[model_mask], color=MODEL_COLOR, lw=2, zorder=4, label="crams")

    ax.set_xscale("log")
    ax.set_xlim(max(0.8, x.min() * 0.8), x.max() * 1.2)
    ax.set_xlabel(r"$R$ [GV]")
    ax.set_ylabel(panel.label)
    if panel.power:
        ax.set_ylim(bottom=0.0)
    ax.set_title(f"{panel.label}   " + rf"$\chi^2$={chi2:.1f}/{n}")
    ax.legend()
    return chi2, n


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot one crams model against Be data")
    p.add_argument("overrides", nargs="*", help="parameter overrides as key=value (e.g. d0=2.5 qne=2e-3)")
    p.add_argument("--ini", default=None, help="crams .ini with parameter values (e.g. bestfit.ini)")
    p.add_argument("--chain", default=None, help="chain .npz; use its posterior median")
    p.add_argument("--output", default="figs/model",
                   help="output prefix; figures are saved as <prefix>_<panel>.pdf")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--fragmentation-model", default=None, choices=FRAGMENTATION_MODELS,
                   help="crams fragmentation model; defaults to the one in --ini, else crams' built-in default")
    args = p.parse_args(argv)

    ini = resolve_params(args)
    frag_model = args.fragmentation_model
    if frag_model is None and args.ini:
        frag_model = read_ini_str(Path(args.ini), "fragmentation_model")
    print(f"Fragmentation model: {frag_model or 'crams default'}")

    runner = CramsRunner(
        build_dir=args.build_dir,
        fragmentation_model=frag_model,
        read_isotopes=True,
    )
    spectra = runner.run(ini)
    if spectra is None:
        sys.exit("crams run failed for the requested parameters.")
    missing = [key for key in ("Be", "B", "Be7", "Be9", "Be10") if key not in spectra]
    if missing:
        sys.exit(f"crams output is missing required spectra: {', '.join(missing)}")

    if STYLE.exists():
        try:
            plt.style.use(str(STYLE))
            if shutil.which("latex") is None:
                plt.rcParams["text.usetex"] = False
        except Exception:
            pass

    prefix = Path(args.output)
    if prefix.suffix:
        prefix = prefix.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    for panel in PANELS:
        fig, ax = plt.subplots(figsize=(11.5, 8.0))
        chi2, n = plot_panel(ax, panel, spectra)
        out = prefix.parent / f"{prefix.name}_{panel.tag}.pdf"
        fig.tight_layout()
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  {panel.tag:<15s} chi2={chi2:8.1f}/{n:<3d} -> {out.name}")

    print(f"Saved {len(PANELS)} figures under {prefix.parent.resolve()}/")


if __name__ == "__main__":
    main()
