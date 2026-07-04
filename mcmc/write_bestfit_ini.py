#!/usr/bin/env python3
"""Write a crams .ini holding the best-fit parameters from an MCMC chain.

The best-fit point is either the posterior median (default) or the maximum
a-posteriori sample (``--map``, the chain point with the highest log-probability).
Active (fitted) parameters take their value from the chain; all remaining crams
parameters take their fixed value from ``run_mcmc.PARAMETERS``, so the resulting
.ini reproduces exactly the model the fit assumed.

Usage
-----
    python write_bestfit_ini.py crams_chain.npz                 # -> crams_chain.ini
    python write_bestfit_ini.py crams_chain.npz -o bestfit.ini
    python write_bestfit_ini.py crams_chain.npz --map --discard 2000
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from run_mcmc import PARAMETERS
from runner import to_ini_params


def best_fit_params(
    chain_file: str | Path,
    *,
    estimator: str = "median",
    discard: int = 0,
    thin: int = 1,
) -> dict[str, float]:
    """Return the full crams parameter dict at the chain best-fit point.

    *estimator* is ``"median"`` (per-parameter posterior median) or ``"map"``
    (the single sample with the largest log-probability).
    """
    data = np.load(Path(chain_file), allow_pickle=True)
    chain = data["chain"]
    param_names = [str(name) for name in data["param_names"]]

    if discard:
        chain = chain[discard:]
    if thin > 1:
        chain = chain[::thin]
    if chain.size == 0:
        raise ValueError("chain is empty after discard/thin")

    if estimator == "median":
        point = np.median(chain, axis=0)
    elif estimator == "map":
        if "log_prob" not in data:
            raise ValueError("chain file has no 'log_prob'; cannot use estimator='map'")
        log_prob = data["log_prob"]
        if discard:
            log_prob = log_prob[discard:]
        if thin > 1:
            log_prob = log_prob[::thin]
        point = chain[int(np.argmax(log_prob))]
    else:
        raise ValueError(f"unknown estimator {estimator!r}; choose 'median' or 'map'")

    # Start from the full crams parameter set (fixed values), then overwrite the
    # fitted ones with the chain best-fit.
    ini = {p.name: p.value for p in PARAMETERS}
    ini.update(dict(zip(param_names, point)))
    return ini


def write_bestfit_ini(
    chain_file: str | Path,
    output: str | Path | None = None,
    *,
    estimator: str = "median",
    discard: int = 0,
    thin: int = 1,
) -> Path:
    """Write the best-fit crams .ini and return its path."""
    chain_path = Path(chain_file)
    out_path = Path(output) if output else chain_path.with_suffix(".ini")

    ini = best_fit_params(chain_file, estimator=estimator, discard=discard, thin=thin)

    with open(out_path, "w") as f:
        f.write(f"# crams best-fit parameters from {chain_path.name}\n")
        f.write(f"# estimator: {estimator}  discard: {discard}  thin: {thin}\n")
        for key, value in to_ini_params(ini).items():
            f.write(f"{key} {value:.6e}\n")
        f.write("id 0\n")
    return out_path


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Write a crams .ini with the MCMC best-fit parameters")
    p.add_argument("chain_file", help=".npz file from run_mcmc.py")
    p.add_argument("-o", "--output", default=None, help="output .ini (default: <chain_file>.ini)")
    p.add_argument("--map", action="store_true", help="use the maximum-a-posteriori sample instead of the median")
    p.add_argument("--discard", type=int, default=0, help="discard first N samples")
    p.add_argument("--thin", type=int, default=1, help="thin chain by N")
    args = p.parse_args(argv)

    estimator = "map" if args.map else "median"
    out_path = write_bestfit_ini(
        args.chain_file,
        args.output,
        estimator=estimator,
        discard=args.discard,
        thin=args.thin,
    )

    ini = best_fit_params(args.chain_file, estimator=estimator, discard=args.discard, thin=args.thin)
    print(f"Best-fit ({estimator}) written to {out_path.resolve()}")
    for key, value in ini.items():
        print(f"  {key:<8s} {value:.6g}")


if __name__ == "__main__":
    main()
