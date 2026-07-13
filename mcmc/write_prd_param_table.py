#!/usr/bin/env python3
"""Write PRD best-fit parameter tables (plain text + LaTeX) per model.

For each of the Evoli et al. 2026 W93 and ST99 chains of a scenario this writes
one table per model listing, for every fitted parameter (the per-species source
normalisations ``q*`` are omitted):

    * MAP        -- the maximum-a-posteriori sample, i.e. the chain point with
                    the largest log-probability (same definition as
                    write_bestfit_ini.py --map);
    * Median     -- the marginalised posterior median;
    * 95.45% CI  -- the [2.275, 97.725] percentile credible interval, reported
                    as asymmetric offsets around the median.

Only the chain .npz is read, so no crams build is required.

Usage
-----
    python write_prd_param_table.py
    python write_prd_param_table.py --scenario free_h_preliminary_be
    python write_prd_param_table.py --scenario free_h_beb --discard 2000 --thin 5
    python write_prd_param_table.py --output-dir prd_figs
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import floor, log10
from pathlib import Path

import numpy as np

MCMC_DIR = Path(__file__).resolve().parent
DEFAULT_CHAIN_DIR = MCMC_DIR / "mcmc_chains"
DEFAULT_OUTPUT_DIR = MCMC_DIR / "prd_figs"

# 95.45% (2 sigma) equal-tailed credible interval.
CI_LABEL = "95.45%"
CI_LABEL_TEX = CI_LABEL.replace("%", r"\%")  # '%' is the LaTeX comment char
CI_LOWER_PCT = 2.2750131948
CI_UPPER_PCT = 97.7249868052

MODELS = {
    "w93": {"fragmentation_model": "evoli2026w93", "label": "W93"},
    "st99": {"fragmentation_model": "evoli2026st99", "label": "ST99"},
}

SCENARIOS = (
    "baseline",
    "free_h_beb",
    "free_h_preliminary_be",
)


@dataclass(frozen=True)
class ParamStyle:
    name: str          # chain parameter name
    plain: str         # plain-text label
    latex: str         # LaTeX label
    unit_plain: str    # plain-text unit ("" if none)
    unit_latex: str    # LaTeX unit ("" if none)


# Ordered list of the parameters we report; the per-species source
# normalisations (qh, qhe, ...) are intentionally excluded. Parameters absent
# from a given chain (e.g. h and the Be fudges in the baseline scenario) are
# silently skipped.
PARAM_STYLES = [
    ParamStyle("hslope", "gamma_H", r"$\gamma_{\rm H}$", "", ""),
    ParamStyle("heslope", "gamma_He", r"$\gamma_{\rm He}$", "", ""),
    ParamStyle("slope", "gamma_Z>=3", r"$\gamma_{Z\geq3}$", "", ""),
    ParamStyle("d0_h", "D0/H", r"$D_0/H$",
               "1e28 cm2/s/kpc", r"$10^{28}\,{\rm cm^2\,s^{-1}\,kpc^{-1}}$"),
    ParamStyle("delta", "delta", r"$\delta$", "", ""),
    ParamStyle("ddelta", "Delta_delta", r"$\Delta\delta$", "", ""),
    ParamStyle("rb_log", "log10(Rb/GV)", r"$\log_{10}(R_b/{\rm GV})$", "", ""),
    ParamStyle("va", "vA", r"$v_A$", "km/s", r"${\rm km\,s^{-1}}$"),
    ParamStyle("phi", "phi", r"$\phi$", "GV", r"${\rm GV}$"),
    ParamStyle("h", "H", r"$H$", "kpc", r"${\rm kpc}$"),
    ParamStyle("fudge_be7", "f(7Be)", r"$f_{^{7}{\rm Be}}$", "", ""),
    ParamStyle("fudge_be9", "f(9Be)", r"$f_{^{9}{\rm Be}}$", "", ""),
    ParamStyle("fudge_be10", "f(10Be)", r"$f_{^{10}{\rm Be}}$", "", ""),
]


@dataclass(frozen=True)
class ParamRow:
    style: ParamStyle
    map: float
    median: float
    err_lo: float   # median - lower bound (>= 0)
    err_hi: float   # upper bound - median (>= 0)
    decimals: int


@dataclass(frozen=True)
class ModelTable:
    model_key: str
    label: str
    fragmentation_model: str
    scenario: str
    halo_size: float
    n_samples: int
    acceptance: float
    rows: list[ParamRow]


def _format_halo_tag(halosize: float) -> str:
    if float(halosize).is_integer():
        return str(int(halosize))
    return f"{halosize:g}".replace(".", "p")


def _scenario_chain_path(chain_dir: Path, model_key: str, scenario: str,
                         halosize: float) -> Path:
    model = MODELS[model_key]["fragmentation_model"]
    halo_tag = _format_halo_tag(halosize)
    return chain_dir / f"mcmc_{model}_{scenario}_h{halo_tag}.npz"


def _uncertainty_decimals(err_lo: float, err_hi: float, sig: int = 2) -> int:
    """Decimal places so the smaller error bar shows *sig* significant figures."""
    err = min(err_lo, err_hi)
    if not np.isfinite(err) or err <= 0.0:
        return 4
    return max(0, -(int(floor(log10(err))) - (sig - 1)))


def _build_model_table(path: Path, model_key: str, scenario: str,
                       discard: int, thin: int) -> ModelTable:
    data = np.load(path, allow_pickle=True)
    chain = data["chain"]
    param_names = [str(name) for name in data["param_names"]]
    log_prob = data["log_prob"]

    if discard:
        chain = chain[discard:]
        log_prob = log_prob[discard:]
    if thin > 1:
        chain = chain[::thin]
        log_prob = log_prob[::thin]
    if chain.size == 0:
        raise ValueError(f"{path}: chain empty after discard/thin")

    map_point = chain[int(np.argmax(log_prob))]
    index = {name: i for i, name in enumerate(param_names)}

    rows = []
    for style in PARAM_STYLES:
        if style.name not in index:
            continue
        col = index[style.name]
        values = chain[:, col]
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue
        lo, median, hi = np.percentile(
            finite, [CI_LOWER_PCT, 50.0, CI_UPPER_PCT])
        err_lo = float(median - lo)
        err_hi = float(hi - median)
        rows.append(ParamRow(
            style=style,
            map=float(map_point[col]),
            median=float(median),
            err_lo=err_lo,
            err_hi=err_hi,
            decimals=_uncertainty_decimals(err_lo, err_hi),
        ))

    halo_size = float(data["halo_size"]) if "halo_size" in data else float("nan")
    return ModelTable(
        model_key=model_key,
        label=MODELS[model_key]["label"],
        fragmentation_model=MODELS[model_key]["fragmentation_model"],
        scenario=scenario,
        halo_size=halo_size,
        n_samples=int(chain.shape[0]),
        acceptance=float(data["acceptance"]) if "acceptance" in data else float("nan"),
        rows=rows,
    )


def _fmt(value: float, decimals: int) -> str:
    return f"{value:.{decimals}f}"


def _plain_table(table: ModelTable) -> str:
    lines = []
    lines.append(
        f"Model: {table.label}  ({table.fragmentation_model}, "
        f"scenario={table.scenario}, H={table.halo_size:g} kpc)"
    )
    lines.append(
        f"  samples={table.n_samples:,}  acceptance={table.acceptance:.3f}  "
        f"central value = posterior median; interval = {CI_LABEL} CI"
    )
    header = f"{'Parameter':<14} {'Unit':<26} {'MAP':>12} {'Median':>12} {'-'+CI_LABEL:>11} {'+'+CI_LABEL:>11}"
    lines.append(header)
    lines.append("-" * len(header))
    for row in table.rows:
        d = row.decimals
        lines.append(
            f"{row.style.plain:<14} {row.style.unit_plain:<26} "
            f"{_fmt(row.map, d):>12} {_fmt(row.median, d):>12} "
            f"{_fmt(row.err_lo, d):>11} {_fmt(row.err_hi, d):>11}"
        )
    return "\n".join(lines)


def _latex_table(table: ModelTable) -> str:
    scenario_tex = table.scenario.replace("_", r"\_")
    lines = []
    lines.append(r"\begin{table}")
    lines.append(r"  \centering")
    lines.append(
        r"  \caption{Best-fit parameters for the "
        f"{table.label} "
        r"fragmentation model (scenario \texttt{"
        f"{scenario_tex}"
        r"}, $H="
        f"{table.halo_size:g}"
        r"\,\mathrm{kpc}$). Central values are posterior medians with "
        f"{CI_LABEL_TEX} "
        r"credible intervals; the MAP column is the maximum-a-posteriori "
        r"sample.}"
    )
    lines.append(r"  \begin{tabular}{lccc}")
    lines.append(r"    \hline\hline")
    lines.append(
        r"    Parameter & Unit & MAP & Median (" + CI_LABEL_TEX + r" CI) \\")
    lines.append(r"    \hline")
    for row in table.rows:
        d = row.decimals
        unit = row.style.unit_latex if row.style.unit_latex else "--"
        median_cell = (
            f"${_fmt(row.median, d)}"
            f"^{{+{_fmt(row.err_hi, d)}}}_{{-{_fmt(row.err_lo, d)}}}$"
        )
        lines.append(
            f"    {row.style.latex} & {unit} & "
            f"${_fmt(row.map, d)}$ & {median_cell} \\\\"
        )
    lines.append(r"    \hline\hline")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Write PRD best-fit parameter tables (plain text + LaTeX)")
    p.add_argument("--scenario", choices=SCENARIOS, default="baseline",
                   help="MCMC scenario to tabulate")
    p.add_argument("--halosize", type=float, default=7.0,
                   help="halo tag used in default chain names")
    p.add_argument("--chain-dir", type=Path, default=DEFAULT_CHAIN_DIR,
                   help="directory containing the MCMC .npz chains")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                   help="directory for the output .txt and .tex tables")
    p.add_argument("--discard", type=int, default=0,
                   help="discard first N flattened samples before tabulating")
    p.add_argument("--thin", type=int, default=1,
                   help="thin the flattened chain by this factor")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)
    scenario = args.scenario
    tables = []
    for model_key in ("w93", "st99"):
        path = _scenario_chain_path(args.chain_dir, model_key, scenario, args.halosize)
        table = _build_model_table(path, model_key, scenario, args.discard, args.thin)
        tables.append(table)

    halo_tag = _format_halo_tag(args.halosize)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    txt_path = args.output_dir / f"params_{scenario}_h{halo_tag}.txt"
    tex_path = args.output_dir / f"params_{scenario}_h{halo_tag}.tex"

    plain_blocks = [_plain_table(t) for t in tables]
    latex_blocks = [_latex_table(t) for t in tables]

    txt_content = ("\n\n".join(plain_blocks)) + "\n"
    tex_content = (
        f"% PRD best-fit parameter tables -- scenario {scenario}, H{halo_tag}\n"
        + "\n\n".join(latex_blocks) + "\n"
    )
    txt_path.write_text(txt_content)
    tex_path.write_text(tex_content)

    print(txt_content)
    print(f"Saved: {txt_path.resolve()}")
    print(f"Saved: {tex_path.resolve()}")


if __name__ == "__main__":
    main()
