#!/usr/bin/env python3
"""Plot one crams model against ALL datasets in run_mcmc.DATASETS.

A quick inspector for a best-fit (or any parameter set): runs crams once and
overlays the model on every dataset used in the fit — fluxes and ratios — in a
single grid, with the fit range highlighted and the chi^2 per dataset. Points
outside the fit range are shown faded for context.

Parameters come from (in increasing precedence): the run_mcmc defaults, a crams
.ini (--ini, e.g. bestfit.ini), a chain median (--chain), and finally any
key=value overrides on the command line — so you can tweak values by hand and
re-plot.

Usage
-----
    python plot_model.py --ini bestfit.ini
    python plot_model.py --ini bestfit.ini qne=2.0e-3 d0=2.5
    python plot_model.py --chain crams_chain.npz --output figs/model.pdf
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fitting import KISS_DIR, _chi2_dataset, _read_kiss_table
from run_mcmc import DATASETS, PARAMETERS
from runner import CramsRunner

STYLE = Path(__file__).parent.parent.parent / "crams-plots" / "crams.mplstyle"
POWER = 2.7
DATA_COLOR = "tab:blue"
MODEL_COLOR = "tab:red"


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
            continue  # e.g. "inelastic_model tripathi99"
    return vals


def read_ini_str(path: Path, key: str) -> str | None:
    """Return the value of a non-numeric 'key value' line (e.g. fragmentation_model)."""
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


# Some quantities have two equivalent representations: the physical key crams
# writes (d0, rb) and the fit-space key run_mcmc samples (d0_h = d0/h,
# rb_log = log10 rb). They must never coexist in the param dict, or
# runner.to_ini_params would let the fit-space value override the physical one.
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


def model_observable(spectra: dict, dataset) -> np.ndarray:
    num = spectra[dataset.numerator]
    if dataset.denominator:
        den = spectra[dataset.denominator]
        y = np.full_like(num, np.nan)
        np.divide(num, den, out=y, where=den > 0)
        return y
    return num


def plot_panel(ax, dataset, spectra, cache) -> float:
    is_ratio = bool(dataset.denominator)
    power = 0.0 if is_ratio else POWER
    label = dataset.numerator + ("/" + dataset.denominator if is_ratio else "")

    # ── data: in-fit points solid, excluded points faded ───────────────────────
    x, y, lo, hi = _read_kiss_table(KISS_DIR / dataset.filename)
    in_fit = (x >= dataset.R_min) & (x <= dataset.R_max)
    sc = x ** power
    ax.errorbar(x[in_fit], (sc * y)[in_fit], yerr=[(sc * lo)[in_fit], (sc * hi)[in_fit]],
                fmt="o", ms=4, color=DATA_COLOR, elinewidth=1.0, capsize=0, zorder=3)
    if (~in_fit).any():
        ax.errorbar(x[~in_fit], (sc * y)[~in_fit], yerr=[(sc * lo)[~in_fit], (sc * hi)[~in_fit]],
                    fmt="o", ms=3, color=DATA_COLOR, alpha=0.25, elinewidth=0.8, capsize=0, zorder=2)

    # ── model on its own grid, over the data span ──────────────────────────────
    R = spectra["R"]
    ym = model_observable(spectra, dataset)
    span = (R >= x.min() * 0.9) & (R <= x.max() * 1.1)
    ax.plot(R[span], (R ** power * ym)[span], color=MODEL_COLOR, lw=2, zorder=4)

    if dataset.R_min > x.min():
        ax.axvline(dataset.R_min, color="gray", ls=":", lw=1)

    chi2 = _chi2_dataset(spectra, dataset, cache)
    n = int(np.sum(in_fit & (lo > 0) & (hi > 0)))
    ax.set_xscale("log")
    if is_ratio:
        ax.set_ylabel(label)
    else:
        ax.set_yscale("log")
        ax.set_ylabel(r"$R^{2.7}\,\Phi$  [GV$^{1.7}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax.set_title(rf"{label}   $\chi^2$={chi2:.0f}/{n}   $w$={dataset.weight:g}")
    ax.set_xlabel(r"$R$ [GV]")
    return chi2


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Plot one crams model against all run_mcmc datasets")
    p.add_argument("overrides", nargs="*", help="parameter overrides as key=value (e.g. d0=2.5 qne=2e-3)")
    p.add_argument("--ini", default=None, help="crams .ini with parameter values (e.g. bestfit.ini)")
    p.add_argument("--chain", default=None, help="chain .npz; use its posterior median")
    p.add_argument("--output", default="figs/model", help="output prefix; each panel is saved as <prefix>_<species>.pdf")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--fragmentation-model", default=None,
                   help="crams fragmentation model; defaults to the one in --ini, else crams' built-in default")
    args = p.parse_args(argv)

    ini = resolve_params(args)
    # The fragmentation model is a string (not a fit parameter), so it isn't in
    # `ini`; read it from the .ini unless overridden, and pass it to the runner.
    # Otherwise crams falls back to its built-in default and the chi^2 won't match
    # the best-fit.
    frag_model = args.fragmentation_model
    if frag_model is None and args.ini:
        frag_model = read_ini_str(Path(args.ini), "fragmentation_model")
    print(f"Fragmentation model: {frag_model or 'crams default'}")
    runner = CramsRunner(build_dir=args.build_dir, fragmentation_model=frag_model)
    spectra = runner.run(ini)
    if spectra is None:
        sys.exit("crams run failed for the requested parameters.")

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

    cache: dict = {}
    total = 0.0
    print("Per-dataset chi^2 (unweighted):")
    for dataset in DATASETS:
        fig, ax = plt.subplots(figsize=(11.5, 8.0))
        chi2 = plot_panel(ax, dataset, spectra, cache)
        total += dataset.weight * chi2

        tag = dataset.numerator + ("_" + dataset.denominator if dataset.denominator else "")
        label = dataset.numerator + ("/" + dataset.denominator if dataset.denominator else "")
        out = prefix.parent / f"{prefix.name}_{tag}.pdf"
        fig.tight_layout()
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  {label:<6s} chi2={chi2:8.1f}  w={dataset.weight:g}  -> {out.name}")

    print(f"Total (weighted) chi^2 = {total:.1f}")
    print(f"Saved {len(DATASETS)} figures under {prefix.parent.resolve()}/")


if __name__ == "__main__":
    main()
