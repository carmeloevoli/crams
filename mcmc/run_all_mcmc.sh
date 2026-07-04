#!/usr/bin/env bash
#
# Run run_mcmc.py for every crams fragmentation cross-section model, one after
# another, starting each chain from the per-model baseline best-fit seed
# produced by run_all_bestfits.sh:
#   python find_bestfit.py --scenario baseline --fragmentation-model <model> ...
# saved as bestfits/bestfit_<model>_h<HALOSIZE>.ini.
#
# Scenario: baseline — published AMS-02 fluxes/ratios only, halo half-height h
# fixed at HALOSIZE kpc (no free h, no Be isotope data, no fudge_be7/9/10
# nuisance factors). HALOSIZE is passed explicitly via --halosize so the chain
# matches the fixed-h best-fit seed.
#
# Edit the configuration below to change the MCMC settings.
#
set -euo pipefail

cd "$(dirname "$0")"

# ── Configuration ─────────────────────────────────────────────────────────────
HALOSIZE=7        # halo half-height h [kpc]; must match run_all_bestfits.sh
NWALKERS=96
NBURN=400
NSTEPS=6000
NCORES=16
SEED=42
# ──────────────────────────────────────────────────────────────────────────────

MODELS=(
  evoli2026w93
  evoli2026st99
)

SCENARIO="baseline"

BESTFITDIR="bestfits"   # per-model baseline best-fit seeds (bestfit_<model>_h<HALOSIZE>.ini)
OUTDIR="mcmc_chains"
mkdir -p "$OUTDIR"

for model in "${MODELS[@]}"; do
  run_tag="${model}_${SCENARIO}_h${HALOSIZE}"
  start="$BESTFITDIR/bestfit_${model}_h${HALOSIZE}.ini"   # baseline best-fit seed (find_bestfit)
  out="$OUTDIR/mcmc_${run_tag}.npz"                        # MCMC chain
  log="$OUTDIR/mcmc_${run_tag}.log"
  echo "=============================================================="
  echo ">>> MCMC with fragmentation model: ${model}"
  echo "    scenario=${SCENARIO}  h=${HALOSIZE}  nwalkers=${NWALKERS}  nburn=${NBURN}  nsteps=${NSTEPS}  ncores=${NCORES}"
  echo "    start:  ${start}"
  echo "    output: ${out}"
  echo "    log:    ${log}"
  echo "=============================================================="

  if [[ ! -f "${start}" ]]; then
    echo "!!! missing best-fit seed ${start} — run run_all_bestfits.sh first; skipping ${model}" >&2
    continue
  fi

  python3 run_mcmc.py \
    --scenario "${SCENARIO}" \
    --fragmentation-model "${model}" \
    --halosize "${HALOSIZE}" \
    --start "${start}" \
    --nwalkers "${NWALKERS}" \
    --nburn "${NBURN}" \
    --nsteps "${NSTEPS}" \
    --ncores "${NCORES}" \
    --seed "${SEED}" \
    --output "${out}" \
    2>&1 | tee "${log}"
done

echo
echo "=============================================================="
echo "Summary (mean acceptance fraction per model, scenario=${SCENARIO}, h=${HALOSIZE}):"
for model in "${MODELS[@]}"; do
  log="$OUTDIR/mcmc_${model}_${SCENARIO}_h${HALOSIZE}.log"
  line=$(grep -E "^Mean acceptance fraction" "${log}" 2>/dev/null | tail -1 || true)
  printf "  %-26s %s\n" "${model}" "${line:-<no chain found, check ${log}>}"
done
