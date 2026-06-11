#!/usr/bin/env python3
"""MCMC fit of the crams propagation model to cosmic-ray data.

Usage
-----
    python run_mcmc.py                          # default settings
    python run_mcmc.py --nwalkers 64 --nsteps 1000 --nburn 300
    python run_mcmc.py --output my_chain.npz

Customisation
-------------
Edit the PARAMETERS and DATASETS lists below to add, remove, or freeze
parameters and to choose which datasets enter the fit.

Output
------
An .npz file with keys:
  chain       – flat chain, shape (nwalkers * nsteps, ndim)
  log_prob    – log-posterior for each sample
  param_names – active parameter names
  acceptance  – mean acceptance fraction

To plot results afterwards use corner.py:
  import corner, numpy as np
  d = np.load("mcmc_chain.npz", allow_pickle=True)
  corner.corner(d["chain"], labels=d["param_names"])
"""
from __future__ import annotations

import argparse
import multiprocessing
import sys
from pathlib import Path

import numpy as np

try:
    import emcee
except ImportError:
    sys.exit("emcee is required. Install with:  pip install emcee")

from runner import CramsRunner
from fitting import (
    Dataset,
    Parameter,
    log_posterior,
    pack_theta,
)

# ── Parameters ─────────────────────────────────────────────────────────────────
# name   : key used in params.ini
# value  : initial value (in the units crams expects)
# prior_lo / prior_hi : flat prior bounds
# active : set False to hold fixed at *value*
#
# crams .ini units:
#   qh, qhe, qc, qo   – injection abundance (dimensionless)
#   hslope, heslope   – injection spectral index for H and He
#   slope             – common injection spectral index for nuclei Z≥3 (C, O, …)
#   phi               – solar modulation potential [GV]
#   d0                – diffusion coefficient [units of 1e28 cm²/s]
#   delta, ddelta     – diffusion spectral index and break amplitude
#   va                – Alfvén speed [km/s]
#   rb                – diffusion break rigidity [GV]
#   h                 – halo half-height [kpc]
PARAMETERS: list[Parameter] = [
    # --- free parameters ---
    Parameter("qh",      5.07e-2,    1e-2,  2e-1,  active=True),   # H injection abundance
    Parameter("qhe",     2.54e-2,    5e-3,  1e-1,  active=True),   # He injection abundance
    Parameter("qc",      3.98879e-3, 1e-3,  2e-2,  active=True),   # C injection abundance
    Parameter("qo",      7.15129e-3, 1e-3,  3e-2,  active=True),   # O injection abundance
    Parameter("qfe",     7.15129e-3, 1e-3,  3e-2,  active=True),   # Fe injection abundance
    Parameter("hslope",  4.37,       4.0,   4.8,   active=True),   # H spectral index
    Parameter("heslope", 4.31,       4.0,   4.8,   active=True),   # He spectral index
    Parameter("slope",   4.33,       4.0,   4.8,   active=True),   # nuclei common spectral index
    Parameter("phi",     0.488,      0.1,   1.0,   active=True),   # solar modulation [GV]
    Parameter("rb",      290.0,      100.,  600.,  active=True),   # diffusion break rigidity [GV]
    # --- fixed propagation parameters ---
    Parameter("d0",      2.48,    0.5,   6.0,   active=False),
    Parameter("delta",   0.565,   0.3,   0.8,   active=False),
    Parameter("va",      4.41,    1.0,   15.0,  active=False),
    Parameter("h",       7.0,     1.0,   15.0,  active=False),
    Parameter("ddelta",  0.22,    0.0,   0.5,   active=False),
]

# ── Datasets ───────────────────────────────────────────────────────────────────
# filename    : file inside mcmc/kiss_tables/
# numerator   : element symbol (e.g. 'H')
# denominator : element symbol for ratio, or '' for absolute flux
# R_min/R_max : rigidity range [GV] included in the chi²
# weight      : relative weight of this dataset in the total chi²
DATASETS: list[Dataset] = [
    Dataset("AMS-02_H_rigidity.txt",  "H",  "", R_min=10.0, R_max=1500.0, weight=1.0),
    Dataset("AMS-02_He_rigidity.txt", "He", "", R_min=10.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_C_rigidity.txt",  "C",  "", R_min=10.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_O_rigidity.txt",  "O",  "", R_min=10.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_Fe_rigidity.txt", "Fe", "", R_min=10.0, R_max=2500.0, weight=1.0),
]

# ── MCMC defaults ──────────────────────────────────────────────────────────────
N_WALKERS = 32
N_BURN    = 200
N_STEPS   = 500


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="MCMC fit of crams propagation model to cosmic-ray data"
    )
    p.add_argument("--nwalkers",  type=int, default=N_WALKERS, help="number of emcee walkers")
    p.add_argument("--nburn",     type=int, default=N_BURN,    help="burn-in steps (discarded)")
    p.add_argument("--nsteps",    type=int, default=N_STEPS,   help="production steps")
    p.add_argument("--output",    default="mcmc_chain.npz",    help="output .npz file")
    p.add_argument("--build-dir", default=None,                help="path to crams build/")
    p.add_argument("--seed",      type=int, default=42,        help="random seed")
    p.add_argument("--ncores",    type=int, default=1,
                   help="parallel worker processes (default: 1 = serial)")
    return p.parse_args(argv)


def _print_summary(flat_chain: np.ndarray, param_names: list[str]) -> None:
    print("\nPosterior summary (median ± 1σ):")
    for i, name in enumerate(param_names):
        lo, mid, hi = np.percentile(flat_chain[:, i], [16, 50, 84])
        print(f"  {name:<10s}: {mid:.5g}  +{hi - mid:.3g} / -{mid - lo:.3g}")


def main(argv=None) -> None:
    args = _parse_args(argv)
    rng = np.random.default_rng(args.seed)

    runner = CramsRunner(build_dir=args.build_dir)
    data_cache: dict = {}

    active = [p for p in PARAMETERS if p.active]
    ndim   = len(active)
    theta0 = pack_theta(PARAMETERS)

    if ndim == 0:
        sys.exit("No active parameters. Set active=True for at least one parameter.")
    if args.nwalkers < 2 * ndim:
        sys.exit(f"nwalkers ({args.nwalkers}) must be at least 2 × ndim ({2 * ndim}).")

    ncores = min(args.ncores, args.nwalkers)
    print(f"Active parameters ({ndim}): {[p.name for p in active]}")
    print(f"Fixed parameters: {[p.name for p in PARAMETERS if not p.active]}")
    print(f"Datasets ({len(DATASETS)}): {[d.filename for d in DATASETS]}")
    print(f"Walkers: {args.nwalkers}  Burn-in: {args.nburn}  Production: {args.nsteps}  Cores: {ncores}")

    # Initialise walkers as a tight Gaussian ball around the starting point
    spread = np.array([0.01 * abs(p.prior_hi - p.prior_lo) for p in active])
    p0 = theta0[None, :] + spread[None, :] * rng.standard_normal((args.nwalkers, ndim))

    pool = multiprocessing.Pool(ncores) if ncores > 1 else None
    sampler = emcee.EnsembleSampler(
        args.nwalkers,
        ndim,
        log_posterior,
        args=(PARAMETERS, DATASETS, runner, data_cache),
        pool=pool,
    )

    print(f"\nBurning in ({args.nburn} steps)…")
    try:
        state = sampler.run_mcmc(p0, args.nburn, progress=True, skip_initial_state_check=True)
    except Exception:
        sampler.reset()
        state = sampler.run_mcmc(p0, args.nburn, progress=False, skip_initial_state_check=True)
    sampler.reset()

    print(f"Production run ({args.nsteps} steps)…")
    try:
        sampler.run_mcmc(state, args.nsteps, progress=True)
    except Exception:
        sampler.reset()
        sampler.run_mcmc(state, args.nsteps, progress=False)

    if pool is not None:
        pool.close()
        pool.join()

    flat_chain = sampler.get_chain(flat=True)
    log_probs  = sampler.get_log_prob(flat=True)
    acceptance = float(np.mean(sampler.acceptance_fraction))
    print(f"\nMean acceptance fraction: {acceptance:.3f}")

    param_names = [p.name for p in active]
    _print_summary(flat_chain, param_names)

    output_path = Path(args.output)
    np.savez(
        output_path,
        chain=flat_chain,
        log_prob=log_probs,
        param_names=np.array(param_names),
        acceptance=acceptance,
    )
    print(f"\nChain saved to {output_path.resolve()}")
    print(
        f"Samples: {flat_chain.shape[0]}  "
        f"({args.nwalkers} walkers × {args.nsteps} steps)"
    )


if __name__ == "__main__":
    main()
