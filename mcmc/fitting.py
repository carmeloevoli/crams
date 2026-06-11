"""Parameter and dataset definitions, kiss-table I/O, and likelihood functions."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from runner import CramsRunner

KISS_DIR = Path(__file__).parent / "kiss_tables"


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass
class Parameter:
    """A model parameter that may be free or fixed.

    *name* is the key written to the .ini file (e.g. ``'delta'``, ``'d0'``).
    *value* is in the units that crams expects:
      - d0    : 1e28 cm²/s   (e.g. 2.48)
      - h     : kpc           (e.g. 7.0)
      - va    : km/s          (e.g. 4.41)
      - rb    : GV            (e.g. 290.0)
      - phi   : GV            (e.g. 0.488)
      - delta, ddelta, xs : dimensionless / g cm⁻²
    *prior_lo*, *prior_hi* define the flat (uniform) prior.
    Set *active = False* to hold the parameter fixed at *value*.
    """
    name: str
    value: float
    prior_lo: float
    prior_hi: float
    active: bool = True


@dataclass
class Dataset:
    """One cosmic-ray dataset from kiss_tables to include in the fit.

    *filename* : file inside kiss_tables/ (e.g. ``'AMS-02_B_C_rigidity.txt'``).
    *numerator* : element symbol for the model quantity (e.g. ``'B'``).
    *denominator* : element symbol for the denominator; ``''`` = absolute flux.
    *R_min*, *R_max* : rigidity range [GV] used in the chi².
    *weight* : relative weight of this dataset in the total chi².
    *error_mode* : how stat and sys errors are combined —
        ``'quadrature'`` (default), ``'stat'``, ``'sys'``, or ``'linear'``.
    """
    filename: str
    numerator: str
    denominator: str = ""
    R_min: float = 3.0
    R_max: float = 1000.0
    weight: float = 1.0
    error_mode: str = "quadrature"


# ── Data loading ───────────────────────────────────────────────────────────────

def _read_kiss_table(
    path: Path,
    error_mode: str = "quadrature",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (x, y, err_lo, err_hi) from a CRDB/KISS ASCII table.

    Columns: x  y  stat_lo  stat_hi  sys_lo  sys_hi

    error_mode controls how stat and sys are combined:
      'quadrature' – err = sqrt(stat² + sys²)  [default]
      'stat'       – err = stat only
      'sys'        – err = sys only
      'linear'     – err = stat + sys

    err_lo is the downward uncertainty, err_hi the upward uncertainty.
    The caller selects which to use based on the sign of (model - data).
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    x, y = data[:, 0], data[:, 1]
    if data.shape[1] >= 6:
        stat_lo, stat_hi = np.abs(data[:, 2]), np.abs(data[:, 3])
        sys_lo,  sys_hi  = np.abs(data[:, 4]), np.abs(data[:, 5])
        if error_mode == "quadrature":
            err_lo = np.hypot(stat_lo, sys_lo)
            err_hi = np.hypot(stat_hi, sys_hi)
        elif error_mode == "stat":
            err_lo, err_hi = stat_lo, stat_hi
        elif error_mode == "sys":
            err_lo, err_hi = sys_lo, sys_hi
        elif error_mode == "linear":
            err_lo = stat_lo + sys_lo
            err_hi = stat_hi + sys_hi
        else:
            raise ValueError(f"unknown error_mode '{error_mode}'; "
                             "choose 'quadrature', 'stat', 'sys', or 'linear'")
    else:
        err_lo, err_hi = np.abs(data[:, 2]), np.abs(data[:, 3])
    return x, y, err_lo, err_hi


# ── Model evaluation ───────────────────────────────────────────────────────────

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


def compute_observable(
    spectra: dict[str, np.ndarray],
    dataset: Dataset,
    R_data: np.ndarray,
) -> np.ndarray:
    """Return model prediction for *dataset* interpolated at *R_data*."""
    R_model = spectra["R"]
    num = spectra[dataset.numerator]
    if dataset.denominator:
        denom = spectra[dataset.denominator]
        with np.errstate(divide="ignore", invalid="ignore"):
            y_model = np.where(denom > 0, num / denom, np.nan)
    else:
        y_model = num
    return _interpolate_loglog(R_model, y_model, R_data)


# ── Parameter helpers ──────────────────────────────────────────────────────────

def pack_theta(params: list[Parameter]) -> np.ndarray:
    """Collect active parameter values into a 1-D array."""
    return np.array([p.value for p in params if p.active])


def unpack_theta(
    theta: np.ndarray, params: list[Parameter]
) -> dict[str, float]:
    """Map active *theta* values back onto the full .ini parameter dict."""
    ini: dict[str, float] = {p.name: p.value for p in params}
    active = [p for p in params if p.active]
    for val, param in zip(theta, active):
        ini[param.name] = val
    return ini


# ── Likelihood ─────────────────────────────────────────────────────────────────

def log_prior(theta: np.ndarray, params: list[Parameter]) -> float:
    """Flat prior: 0 inside bounds, -inf outside."""
    active = [p for p in params if p.active]
    for val, param in zip(theta, active):
        if not (param.prior_lo <= val <= param.prior_hi):
            return -np.inf
    return 0.0


def _chi2_dataset(
    spectra: dict[str, np.ndarray],
    dataset: Dataset,
    data_cache: dict,
) -> float:
    """Chi² contribution from one dataset."""
    cache_key = (dataset.filename, dataset.error_mode)
    if cache_key not in data_cache:
        data_cache[cache_key] = _read_kiss_table(KISS_DIR / dataset.filename, dataset.error_mode)

    x_all, y_all, err_lo_all, err_hi_all = data_cache[cache_key]
    cut = (x_all >= dataset.R_min) & (x_all <= dataset.R_max) & (err_lo_all > 0) & (err_hi_all > 0)
    x, y = x_all[cut], y_all[cut]
    err_lo, err_hi = err_lo_all[cut], err_hi_all[cut]

    if len(x) == 0:
        return 0.0

    y_model = compute_observable(spectra, dataset, x)
    valid = np.isfinite(y_model)
    if not valid.any():
        return 1e10

    # Use upper error when model > data, lower error when model < data
    residual = y_model[valid] - y[valid]
    sigma = np.where(residual > 0, err_hi[valid], err_lo[valid])
    return float(np.sum((residual / sigma) ** 2))


def log_likelihood(
    theta: np.ndarray,
    params: list[Parameter],
    datasets: list[Dataset],
    runner: CramsRunner,
    data_cache: dict,
) -> float:
    ini_params = unpack_theta(theta, params)
    spectra = runner.run(ini_params)
    if spectra is None:
        return -np.inf
    total_chi2 = sum(
        ds.weight * _chi2_dataset(spectra, ds, data_cache) for ds in datasets
    )
    return -0.5 * total_chi2


def log_posterior(
    theta: np.ndarray,
    params: list[Parameter],
    datasets: list[Dataset],
    runner: CramsRunner,
    data_cache: dict,
) -> float:
    lp = log_prior(theta, params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, params, datasets, runner, data_cache)
