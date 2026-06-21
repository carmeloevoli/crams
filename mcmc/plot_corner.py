#!/usr/bin/env python3
"""Corner plot of MCMC posterior from a chain .npz file.

Usage
-----
    python plot_corner.py h_he_chain.npz
    python plot_corner.py h_he_chain.npz --output corner.pdf
    python plot_corner.py h_he_chain.npz --thin 10 --discard 50
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import corner
import matplotlib.pyplot as plt

# LaTeX-style labels keyed by parameter name
LABELS: dict[str, str] = {
    "qh":      r"$q_\mathrm{H}$",
    "qhe":     r"$q_\mathrm{He}$",
    "hslope":  r"$\gamma_\mathrm{H}$",
    "heslope": r"$\gamma_\mathrm{He}$",
    "phi":     r"$\phi$ [GV]",
    "d0":      r"$D_0\ [10^{28}\ \mathrm{cm^2/s}]$",
    "delta":   r"$\delta$",
    "ddelta":  r"$\Delta\delta$",
    "va":      r"$v_A$ [km/s]",
    "rb":      r"$R_b$ [GV]",
    "h":       r"$H$ [kpc]",
    "xs":      r"$X_s$ [g/cm$^2$]",
}


def load_chain(path: Path, thin: int = 1, discard: int = 0):
    d = np.load(path, allow_pickle=True)
    chain = d["chain"]           # shape: (nsamples, ndim)
    param_names = list(d["param_names"])
    acceptance = float(d["acceptance"])

    if discard > 0:
        chain = chain[discard:]
    if thin > 1:
        chain = chain[::thin]

    return chain, param_names, acceptance


def make_corner(chain, param_names, acceptance, output: Path) -> None:
    labels = [LABELS.get(n, n) for n in param_names]
    ndim = chain.shape[1]

    quantiles = [0.16, 0.50, 0.84]
    fig = corner.corner(
        chain,
        labels=labels,
        quantiles=quantiles,
        show_titles=True,
        title_fmt=".4g",
        title_kwargs={"fontsize": 11},
        label_kwargs={"fontsize": 12},
        truths=None,
    )

    fig.suptitle(
        f"{chain.shape[0]:,} samples · {ndim}D · acceptance {acceptance:.2f}",
        fontsize=11,
        y=1.01,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=150)
    print(f"Saved: {output.resolve()}")

    # Print medians and 1-sigma intervals
    print("\nPosterior summary (median ± 1σ):")
    for i, name in enumerate(param_names):
        lo, mid, hi = np.percentile(chain[:, i], [16, 50, 84])
        print(f"  {name:<10s}: {mid:.5g}  +{hi - mid:.3g} / -{mid - lo:.3g}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Corner plot from MCMC chain .npz")
    p.add_argument("chain_file", help=".npz file produced by run_mcmc.py")
    p.add_argument("--output",  default=None, help="output file (default: <chain_file>.pdf)")
    p.add_argument("--thin",    type=int, default=1,  help="keep every N-th sample")
    p.add_argument("--discard", type=int, default=0,  help="discard first N samples")
    args = p.parse_args(argv)

    chain_path = Path(args.chain_file)
    output = Path(args.output) if args.output else chain_path.with_suffix(".pdf")

    chain, param_names, acceptance = load_chain(chain_path, args.thin, args.discard)
    print(f"Loaded {chain.shape[0]:,} samples, {chain.shape[1]} parameters: {param_names}")

    make_corner(chain, param_names, acceptance, output)


if __name__ == "__main__":
    main()
