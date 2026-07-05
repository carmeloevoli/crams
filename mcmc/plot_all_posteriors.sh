#!/usr/bin/env bash
#
# Produce all posterior figures for every crams MCMC chain by running the
# plotting scripts on each chain .npz saved by run_all_mcmc.sh:
#   plot_mcmc_posterior.py         — posterior predictions vs AMS-02 fluxes/ratios
#   plot_mcmc_params_posterior.py  — parameter posterior diagnostics
#   plot_corner.py                 — full corner plot of the posterior
#   write_bestfit_ini.py           — crams .ini of the posterior-median best fit
#
# Chains live in mcmc_chains/mcmc_<model>_<scenario>_h<HALOSIZE>.npz. Figures are
# written to mcmc_figs/ (kept separate from the chains to avoid confusion),
# producing mcmc_figs/mcmc_<model>_<scenario>_h<HALOSIZE>_<tag>.pdf.
#
# Edit the configuration below to match run_all_mcmc.sh.
#
set -euo pipefail

cd "$(dirname "$0")"

# ── Configuration ─────────────────────────────────────────────────────────────
HALOSIZE=7        # halo half-height h [kpc]; must match run_all_mcmc.sh
# The chains are already flat (nwalkers × nsteps) with burn-in removed at
# sampling time, so DISCARD/THIN act on the post-burn-in flat chain:
#   DISCARD  extra samples trimmed from the front (0 — burn-in already stripped)
#   THIN     keep every THIN-th sample to decorrelate (fewer, cleaner samples)
#   NSAMPLES posterior draws used for the prediction band (plot_mcmc_posterior)
DISCARD=0
THIN=10
NSAMPLES=500
SEED=42
ESTIMATOR=map     # best-fit .ini point: "median" (posterior median) or "map"
                  # (maximum-a-posteriori, the highest-log-prob chain sample)
# ──────────────────────────────────────────────────────────────────────────────

# Pass --map to write_bestfit_ini.py only when ESTIMATOR=map.
BESTFIT_FLAGS=()
[[ "${ESTIMATOR}" == "map" ]] && BESTFIT_FLAGS+=(--map)

MODELS=(
  evoli2026w93
  evoli2026st99
)

SCENARIOS=(
  baseline
  free_h_beb
  free_h_preliminary_be
)

CHAINDIR="mcmc_chains"
FIGDIR="mcmc_figs"
mkdir -p "$FIGDIR"

for model in "${MODELS[@]}"; do
  for SCENARIO in "${SCENARIOS[@]}"; do
    run_tag="${model}_${SCENARIO}_h${HALOSIZE}"
    chain="$CHAINDIR/mcmc_${run_tag}.npz"
    prefix="$FIGDIR/mcmc_${run_tag}"

    echo "=============================================================="
    echo ">>> Posterior plots for model=${model}  scenario=${SCENARIO}  h=${HALOSIZE}"
    echo "    chain:  ${chain}"
    echo "    output: ${prefix}_<tag>.pdf"
    echo "=============================================================="

    if [[ ! -f "${chain}" ]]; then
      echo "!!! missing chain ${chain} — run run_all_mcmc.sh first; skipping" >&2
      continue
    fi

    python3 plot_mcmc_posterior.py "${chain}" \
      --nsamples "${NSAMPLES}" \
      --discard "${DISCARD}" \
      --thin "${THIN}" \
      --seed "${SEED}" \
      --output "${prefix}"

    python3 plot_mcmc_params_posterior.py "${chain}" \
      --discard "${DISCARD}" \
      --thin "${THIN}" \
      --output "${prefix}"

    python3 plot_corner.py "${chain}" \
      --discard "${DISCARD}" \
      --thin "${THIN}" \
      --output "${prefix}_corner.pdf"

    python3 write_bestfit_ini.py "${chain}" \
      --discard "${DISCARD}" \
      --thin "${THIN}" \
      ${BESTFIT_FLAGS[@]+"${BESTFIT_FLAGS[@]}"} \
      --output "${prefix}.ini"
  done
done
