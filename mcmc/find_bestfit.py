#!/usr/bin/env python3
"""Fast pre-fit: locate the best-fit point to seed the MCMC walkers.

Reuses PARAMETERS and DATASETS from run_mcmc and minimises the total chi^2
(= -2 log L) over the active parameters. The resulting point is a good starting
position for run_mcmc.py (much faster than letting the chain burn in from a
guess) and is also written out as a ready-to-run crams .ini.

Backends (--method)
-------------------
    nelder-mead   pure-numpy downhill simplex (default, no extra dependencies)
    iminuit       MINUIT / MIGRAD  (pip install iminuit); also reports errors
    scipy         scipy.optimize.minimize (pip install scipy); Powell by default

Each chi^2 evaluation runs crams once, so the cost is dominated by the number of
function evaluations; --maxfev caps it.

Usage
-----
    python find_bestfit.py
    python find_bestfit.py --method iminuit --output bestfit.ini
    python find_bestfit.py --maxfev 400 --seed 1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from fitting import _chi2_dataset, log_likelihood, unpack_theta
from run_mcmc import (
    DATASETS,
    PARAMETERS,
    SCENARIOS,
    datasets_need_isotopes,
    make_datasets,
    make_parameters,
    set_halo_size,
)
from runner import FRAGMENTATION_MODELS, CramsRunner, from_ini_params, to_ini_params

# Fit configuration; defaults to the baseline scenario and rebuilt in main() once
# the --scenario argument is known (see _set_scenario).
ACTIVE = [p for p in PARAMETERS if p.active]
NAMES = [p.name for p in ACTIVE]
LO = np.array([p.prior_lo for p in ACTIVE])
HI = np.array([p.prior_hi for p in ACTIVE])
X0 = np.array([p.value for p in ACTIVE])

BIG = 1e12  # penalty returned for failed / out-of-bounds evaluations


def _set_scenario(scenario: str) -> None:
    """Rebuild the module-level PARAMETERS/DATASETS and fit arrays for *scenario*."""
    global PARAMETERS, DATASETS, ACTIVE, NAMES, LO, HI, X0
    PARAMETERS = make_parameters(scenario)
    DATASETS = make_datasets(scenario)
    ACTIVE = [p for p in PARAMETERS if p.active]
    NAMES = [p.name for p in ACTIVE]
    LO = np.array([p.prior_lo for p in ACTIVE])
    HI = np.array([p.prior_hi for p in ACTIVE])
    X0 = np.array([p.value for p in ACTIVE])


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


def make_chi2(runner: CramsRunner, cache: dict, bounded: bool = True):
    """Return f(theta) -> total chi^2.

    If *bounded*, parameters are clipped to the run_mcmc priors [LO, HI];
    otherwise they are passed to crams unchanged, so the minimiser can leave the
    assumed range (crams failures on invalid values are penalised by BIG).
    """

    def chi2(theta: np.ndarray) -> float:
        theta = np.asarray(theta, float)
        if bounded:
            theta = np.clip(theta, LO, HI)
        ll = log_likelihood(theta, PARAMETERS, DATASETS, runner, cache)
        if not np.isfinite(ll):
            return BIG
        return -2.0 * ll

    return chi2


# ── Optimisers ──────────────────────────────────────────────────────────────

def _nelder_mead(f, x0, step, maxfev=500, tol=1e-4):
    """Minimal Nelder-Mead downhill simplex (no external dependencies)."""
    n = len(x0)
    sim = np.array([x0] + [x0 + step[i] * np.eye(n)[i] for i in range(n)], float)
    fval = np.array([f(x) for x in sim])
    nfev = n + 1
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    while nfev < maxfev:
        order = np.argsort(fval)
        sim, fval = sim[order], fval[order]
        if fval[-1] - fval[0] <= tol * (abs(fval[0]) + tol):
            break
        centroid = sim[:-1].mean(axis=0)
        xr = centroid + alpha * (centroid - sim[-1])
        fr = f(xr); nfev += 1
        if fval[0] <= fr < fval[-2]:
            sim[-1], fval[-1] = xr, fr
        elif fr < fval[0]:
            xe = centroid + gamma * (xr - centroid)
            fe = f(xe); nfev += 1
            sim[-1], fval[-1] = (xe, fe) if fe < fr else (xr, fr)
        else:
            xc = centroid + rho * (sim[-1] - centroid)
            fc = f(xc); nfev += 1
            if fc < fval[-1]:
                sim[-1], fval[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    sim[i] = sim[0] + sigma * (sim[i] - sim[0])
                    fval[i] = f(sim[i]); nfev += 1
    order = np.argsort(fval)
    return sim[order][0], float(fval[order][0]), nfev, None


def _fit_nelder_mead(f, x0, maxfev, rng, bounded=True):
    step = 0.05 * (HI - LO)
    return _nelder_mead(f, np.array(x0, float), step, maxfev=maxfev)


def _fit_iminuit(f, x0, maxfev, rng, bounded=True):
    try:
        from iminuit import Minuit
    except ImportError:
        sys.exit("iminuit not installed.  pip install iminuit  (or use --method nelder-mead)")

    m = Minuit(lambda *a: f(np.array(a)), *x0, name=NAMES)
    m.errordef = 1.0  # objective is a chi^2
    if bounded:
        for name, lo, hi in zip(NAMES, LO, HI):
            m.limits[name] = (lo, hi)
    m.migrad(ncall=maxfev)
    errors = np.array(m.errors)
    return np.array(m.values), float(m.fval), int(m.nfcn), errors


def _fit_scipy(f, x0, maxfev, rng, bounded=True, scipy_method="Powell"):
    try:
        from scipy.optimize import minimize
    except ImportError:
        sys.exit("scipy not installed.  pip install scipy  (or use --method nelder-mead)")

    res = minimize(
        f, x0, method=scipy_method,
        bounds=list(zip(LO, HI)) if bounded else None,
        options={"maxfev": maxfev} if scipy_method in ("Powell", "Nelder-Mead") else {"maxiter": maxfev},
    )
    return np.asarray(res.x, float), float(res.fun), int(res.nfev), None


OPTIMISERS = {"nelder-mead": _fit_nelder_mead, "iminuit": _fit_iminuit, "scipy": _fit_scipy}


# ── Reporting / output ──────────────────────────────────────────────────────

def _full_ini(theta: np.ndarray) -> dict[str, float]:
    """Full crams parameter dict: fixed defaults overwritten with best-fit actives."""
    return unpack_theta(theta, PARAMETERS)


def _write_ini(ini: dict[str, float], path: Path, method: str, chi2: float,
               fragmentation_model: str | None = None) -> None:
    with open(path, "w") as f:
        f.write(f"# crams best-fit parameters (pre-fit, method={method}, chi2={chi2:.1f})\n")
        for key, value in to_ini_params(ini).items():
            f.write(f"{key} {value:.6e}\n")
        if fragmentation_model is not None:
            f.write(f"fragmentation_model {fragmentation_model}\n")
        f.write("id 0\n")


def _report(theta, chi2, errors, runner, cache):
    print("\nBest-fit parameters (* = outside the run_mcmc prior range):")
    any_outside = False
    for i, name in enumerate(NAMES):
        outside = theta[i] < LO[i] or theta[i] > HI[i]
        any_outside = any_outside or outside
        flag = " *" if outside else "  "
        rng = f"[{LO[i]:.3g}, {HI[i]:.3g}]"
        err = f"  ± {errors[i]:.3g}" if errors is not None else ""
        print(f" {flag}{name:<8s} {theta[i]:.6g}{err}   prior {rng}")
    if any_outside:
        print("  -> some best-fit values lie OUTSIDE the assumed prior range; widen run_mcmc.PARAMETERS.")
    else:
        print("  -> all best-fit values lie inside the assumed prior range.")

    spectra = runner.run(_full_ini(theta))
    ndata = 0
    print("\nPer-dataset chi^2:")
    for d in DATASETS:
        c = _chi2_dataset(spectra, d, cache)
        x, _, lo, hi = cache[(d.source, d.filename, d.csv_column, d.error_mode)]
        n = int(np.sum((x >= d.R_min) & (x <= d.R_max) & (lo > 0) & (hi > 0)))
        ndata += n
        label = d.numerator + ("/" + d.denominator if d.denominator else "")
        print(f"  {label:<6s} chi2={c:9.1f}  ({n} points)")

    dof = max(ndata - len(ACTIVE), 1)
    print(f"\nTotal chi^2 = {chi2:.1f}   dof = {dof}   chi^2/dof = {chi2 / dof:.2f}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Locate the best-fit point to seed the MCMC")
    p.add_argument("--method", choices=list(OPTIMISERS), default="nelder-mead",
                   help="optimiser backend (default: nelder-mead)")
    p.add_argument("--maxfev", type=int, default=600, help="max crams evaluations (default: 600)")
    p.add_argument("--unbounded", action="store_true",
                   help="do not restrict to the run_mcmc prior ranges (find the unconstrained best-fit)")
    p.add_argument("--start", default=None,
                   help="crams .ini whose active values are used as the starting point (e.g. bestfit.ini)")
    p.add_argument("--output", default="bestfit.ini", help="best-fit crams .ini (default: bestfit.ini)")
    p.add_argument("--fragmentation-model", default=None, choices=FRAGMENTATION_MODELS,
                   help="crams fragmentation cross-section model "
                        "(default: None = crams built-in default)")
    p.add_argument("--halosize", type=float, default=None,
                   help="halo half-height h [kpc] (default: PARAMETERS value, 7)")
    p.add_argument("--scenario", default="baseline", choices=SCENARIOS,
                   help="fit scenario (parameters + datasets); see run_mcmc.py "
                        "(default: baseline)")
    p.add_argument("--build-dir", default=None, help="path to crams build/")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)

    _set_scenario(args.scenario)
    if not ACTIVE:
        sys.exit("No active parameters to fit.")

    if args.halosize is not None:
        set_halo_size(args.halosize, PARAMETERS)
    halo_size = next(p.value for p in PARAMETERS if p.name == "h")

    x0 = X0.copy()
    if args.start:
        # crams .ini stores d0/rb; map them back to the fit-space d0_h/rb_log.
        start = from_ini_params(_read_ini_values(Path(args.start)))
        x0 = np.array([start.get(name, x0[i]) for i, name in enumerate(NAMES)], float)

    rng = np.random.default_rng(args.seed)
    runner = CramsRunner(build_dir=args.build_dir,
                         fragmentation_model=args.fragmentation_model,
                         read_isotopes=datasets_need_isotopes(DATASETS))
    cache: dict = {}
    bounded = not args.unbounded
    chi2_fn = make_chi2(runner, cache, bounded=bounded)

    chi2_start = chi2_fn(x0)
    print(f"Active parameters ({len(ACTIVE)}): {NAMES}")
    print(f"Datasets ({len(DATASETS)}): {[d.numerator + ('/' + d.denominator if d.denominator else '') for d in DATASETS]}")
    print(f"Method: {args.method}   maxfev: {args.maxfev}   bounds: {'priors' if bounded else 'NONE (unbounded)'}")
    print(f"Fragmentation model: {args.fragmentation_model or 'crams default'}   halo h: {halo_size} kpc")
    print(f"Start: {args.start or 'run_mcmc defaults'}   chi^2 = {chi2_start:.1f}")

    theta, chi2, nfev, errors = OPTIMISERS[args.method](chi2_fn, x0, args.maxfev, rng, bounded=bounded)
    if bounded:
        # the simplex/optimiser may report a coordinate just outside the limits;
        # report the value that was actually evaluated (clipped into the priors).
        theta = np.clip(theta, LO, HI)
    print(f"\nDone after {nfev} evaluations.  chi^2: {chi2_start:.1f} -> {chi2:.1f}")

    _report(theta, chi2, errors, runner, cache)

    out_path = Path(args.output)
    _write_ini(_full_ini(theta), out_path, args.method, chi2,
               fragmentation_model=args.fragmentation_model)
    print(f"\nBest-fit crams .ini written to {out_path.resolve()}")
    print("Use it to seed run_mcmc.py (update the PARAMETERS initial values) "
          "or run it directly:  ./crams " + str(out_path))


if __name__ == "__main__":
    main()
