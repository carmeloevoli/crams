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

from runner import FRAGMENTATION_MODELS, CramsRunner, from_ini_params
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
#   d0_h              – diffusion / halo height, d0/h [1e28 cm²/s / kpc]; the
#                       runner reconstructs d0 = d0_h * h when writing the .ini,
#                       so h can be changed without re-tuning this prior
#   delta, ddelta     – diffusion spectral index and break amplitude
#   va                – Alfvén speed [km/s]
#   rb_log            – log10(rb/GV); the runner reconstructs rb = 10**rb_log
#                       when writing the .ini (break rigidity is a scale param,
#                       better sampled in log)
#   h                 – halo half-height [kpc]
PARAMETERS: list[Parameter] = [
    # --- free parameters (initial values = best-fit point) ---
    Parameter("qh",      4.14e-2,   1e-2,  2e-1,  active=True),   # H injection abundance
    Parameter("qhe",     2.04e-2,   5e-3,  1e-1,  active=True),   # He injection abundance
    Parameter("qc",      4.00e-3,   1e-3,  2e-2,  active=True),   # C injection abundance
    Parameter("qn",      3.83e-4,   1e-4,  2e-3,  active=True),   # N injection abundance
    Parameter("qo",      7.21e-3,   1e-3,  3e-2,  active=True),   # O injection abundance
    Parameter("qne",     1.35e-3,   3e-4,  5e-3,  active=True),   # Ne injection abundance
    Parameter("qmg",     2.38e-3,   5e-4,  8e-3,  active=True),   # Mg injection abundance
    Parameter("qsi",     2.83e-3,   5e-4,  8e-3,  active=True),   # Si injection abundance
    Parameter("qs",      5.06e-4,   1e-4,  4e-3,  active=True),   # S injection abundance
    Parameter("qfe",     7.53e-3,   1e-3,  3e-2,  active=True),   # Fe injection abundance
    Parameter("hslope",  4.37,      4.1,   4.7,   active=True),   # H spectral index
    Parameter("heslope", 4.30,      4.1,   4.7,   active=True),   # He spectral index
    Parameter("slope",   4.36,      4.1,   4.7,   active=True),   # nuclei common spectral index
    Parameter("phi",     0.47,      0.1,   1.0,   active=True),   # solar modulation [GV]
    # --- diffusion parameters (constrained by B/C) ---
    Parameter("d0_h",    0.34,      0.05,  1.0,   active=True),    # d0/h [1e28 cm²/s / kpc]; d0 = d0_h * h (= 2.376 at h=7)
    Parameter("delta",   0.54,      0.2,   0.8,   active=True),    # diffusion spectral index
    Parameter("ddelta",  0.27,      0.0,   0.5,   active=True),    # low-rigidity diffusion break amplitude
    Parameter("rb_log",  2.26,      2.0,   3.0,   active=True),   # log10(rb/GV); rb = 10**rb_log (= 316.9 GV), prior 100–600 GV
    Parameter("va",      3.32,      1.0,   15.0,  active=True),    # Alfvén speed [km/s]
    # --- fixed propagation parameters ---
    Parameter("h",       5.0,       1.0,   15.0,  active=False),
]

# ── Datasets ───────────────────────────────────────────────────────────────────
# filename    : file inside mcmc/kiss_tables/
# numerator   : element symbol (e.g. 'H')
# denominator : element symbol for ratio, or '' for absolute flux
# R_min/R_max : rigidity range [GV] included in the chi²
# weight      : relative weight of this dataset in the total chi²
DATASETS: list[Dataset] = [
    # Fluxes
    Dataset("AMS-02_H_rigidity.txt",  "H",  "", R_min=5.0, R_max=1500.0, weight=1.0),
    Dataset("AMS-02_He_rigidity.txt", "He", "", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_B_rigidity.txt", "B", "", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_C_rigidity.txt",  "C",  "", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_O_rigidity.txt",  "O",  "", R_min=5.0, R_max=2500.0, weight=1.0),
    # Ratios
    Dataset("AMS-02_H_He_rigidity.txt", "H", "He", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_He_O_rigidity.txt", "He", "O", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_B_C_rigidity.txt", "B", "C", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_B_O_rigidity.txt", "B", "O", R_min=5.0, R_max=2500.0, weight=1.0),
    Dataset("AMS-02_C_O_rigidity.txt", "C", "O", R_min=5.0, R_max=2500.0, weight=1.0),
    # N and the heavier primaries: only above 40 GV and down-weighted — slightly
    # less relevant, and their small error bars would otherwise dominate the fit.
    Dataset("AMS-02_N_rigidity.txt", "N", "", R_min=40.0, R_max=2500.0, weight=0.25),
    Dataset("AMS-02_Ne_rigidity.txt", "Ne", "", R_min=40.0, R_max=2500.0, weight=0.25),
    Dataset("AMS-02_Mg_rigidity.txt", "Mg", "", R_min=40.0, R_max=2500.0, weight=0.25),
    Dataset("AMS-02_Si_rigidity.txt", "Si", "", R_min=40.0, R_max=2500.0, weight=0.25),
    Dataset("AMS-02_S_rigidity.txt", "S", "", R_min=40.0, R_max=2500.0, weight=0.25),
    Dataset("AMS-02_Fe_rigidity.txt", "Fe", "", R_min=40.0, R_max=2500.0, weight=0.25),
]

# ── MCMC defaults ──────────────────────────────────────────────────────────────
N_WALKERS = 96     # ~5x ndim; pilot acceptance ~0.32
N_BURN    = 300    # ~3.5x tau (tau_max ~56 from pilot); walkers start at the best-fit
N_STEPS   = 4000   # ~53x tau -> ~5000 independent samples


def set_halo_size(h_kpc: float) -> None:
    """Override the fixed halo half-height h (kpc) in PARAMETERS, in place.

    Because the diffusion parameter is sampled as d0/h, changing h only rescales
    the reconstructed d0 (d0 = d0_h * h) and leaves the d0_h prior untouched.
    """
    for p in PARAMETERS:
        if p.name == "h":
            p.value = h_kpc
            return
    raise KeyError("no fixed 'h' parameter found in PARAMETERS")


def _read_ini_values(path: Path) -> dict[str, float]:
    """Read numeric 'key value' lines from a crams .ini (skips comments/strings)."""
    vals: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                vals[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return vals


def set_start_from_ini(path: str) -> None:
    """Override the initial value of each active parameter from a crams .ini.

    The .ini stores physical crams keys (d0, rb); from_ini_params maps them back
    to the fit-space parameters (d0_h, rb_log). Parameters absent from the file
    keep their default value. The walkers are seeded in a tight ball around these
    values, so this sets where the MCMC starts (e.g. a MINUIT best-fit point).
    """
    start = from_ini_params(_read_ini_values(Path(path)))
    for p in PARAMETERS:
        if p.active and p.name in start:
            p.value = start[p.name]


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="MCMC fit of crams propagation model to cosmic-ray data"
    )
    p.add_argument("--nwalkers",  type=int, default=N_WALKERS, help="number of emcee walkers")
    p.add_argument("--nburn",     type=int, default=N_BURN,    help="burn-in steps (discarded)")
    p.add_argument("--nsteps",    type=int, default=N_STEPS,   help="production steps")
    p.add_argument("--output",    default="mcmc_chain.npz",    help="output .npz file")
    p.add_argument("--fragmentation-model", default=None, choices=FRAGMENTATION_MODELS,
                   help="crams fragmentation cross-section model "
                        "(default: None = crams built-in default)")
    p.add_argument("--halosize",  type=float, default=None,
                   help="halo half-height h [kpc] (default: PARAMETERS value, 7)")
    p.add_argument("--start",     default=None,
                   help="crams .ini (e.g. a MINUIT bestfit) whose active values "
                        "seed the walkers (default: PARAMETERS values)")
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

    if args.halosize is not None:
        set_halo_size(args.halosize)
    halo_size = next(p.value for p in PARAMETERS if p.name == "h")

    if args.start is not None:
        set_start_from_ini(args.start)

    runner = CramsRunner(build_dir=args.build_dir,
                         fragmentation_model=args.fragmentation_model)
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
    print(f"Fragmentation model: {args.fragmentation_model or 'crams default'}")
    print(f"Halo half-height h: {halo_size} kpc")
    print(f"Start point: {args.start or 'PARAMETERS defaults'}")
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
    print(f"\nMean acceptance fraction: {acceptance:.3f}  (healthy: ~0.2–0.5)")

    # Autocorrelation time: production should be >> tau (rule of thumb: N_STEPS >= 50*tau)
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        tau_max = float(np.nanmax(tau))
        print(f"Autocorrelation time: mean τ = {np.nanmean(tau):.1f}, max τ = {tau_max:.1f} steps")
        print(f"  N_STEPS / max τ = {args.nsteps / tau_max:.0f}  (aim ≳ 50);  "
              f"suggested burn ≈ {2 * tau_max:.0f}, steps ≈ {50 * tau_max:.0f}")
    except Exception as exc:  # chain too short to estimate reliably
        print(f"Autocorrelation time: unavailable ({exc})")

    param_names = [p.name for p in active]
    _print_summary(flat_chain, param_names)

    output_path = Path(args.output)
    np.savez(
        output_path,
        chain=flat_chain,
        log_prob=log_probs,
        param_names=np.array(param_names),
        acceptance=acceptance,
        fragmentation_model=np.array(args.fragmentation_model or ""),
        halo_size=np.array(halo_size),
    )
    print(f"\nChain saved to {output_path.resolve()}")
    print(
        f"Samples: {flat_chain.shape[0]}  "
        f"({args.nwalkers} walkers × {args.nsteps} steps)"
    )


if __name__ == "__main__":
    main()
