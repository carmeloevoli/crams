#!/usr/bin/env bash
#
# Run run_mcmc.py for every crams fragmentation cross-section model, one after
# another, starting each chain from the MINUIT best-fit produced by
# run_all_bestfits.sh (bestfits/bestfit_<model>_h<H>.ini).
#
# Edit the configuration below to change the MCMC settings.
#
set -euo pipefail

cd "$(dirname "$0")"

# ── Configuration ─────────────────────────────────────────────────────────────
HALOSIZE=5
NWALKERS=96
NBURN=300
NSTEPS=4000
NCORES=16
SEED=42
# ──────────────────────────────────────────────────────────────────────────────

MODELS=(
  evoli2026w93
  evoli2026st99
)

BESTFITDIR="bestfits"   # where run_all_bestfits.sh wrote the seed .ini files
OUTDIR="mcmc_chains"
mkdir -p "$OUTDIR"

for model in "${MODELS[@]}"; do
  tag="${model}_h${HALOSIZE}"
  start="$BESTFITDIR/bestfit_${tag}.ini"   # MINUIT best-fit from run_all_bestfits.sh
  out="$OUTDIR/mcmc_${tag}.npz"            # MCMC chain
  log="$OUTDIR/mcmc_${tag}.log"
  echo "=============================================================="
  echo ">>> MCMC with fragmentation model: ${model}"
  echo "    halosize=${HALOSIZE}  nwalkers=${NWALKERS}  nburn=${NBURN}  nsteps=${NSTEPS}  ncores=${NCORES}"
  echo "    start:  ${start}"
  echo "    output: ${out}"
  echo "    log:    ${log}"
  echo "=============================================================="

  if [[ ! -f "${start}" ]]; then
    echo "!!! missing best-fit seed ${start} — run run_all_bestfits.sh first; skipping ${model}" >&2
    continue
  fi

  python3 run_mcmc.py \
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
echo "Summary (mean acceptance fraction per model, h=${HALOSIZE}):"
for model in "${MODELS[@]}"; do
  log="$OUTDIR/mcmc_${model}_h${HALOSIZE}.log"
  line=$(grep -E "^Mean acceptance fraction" "${log}" 2>/dev/null | tail -1 || true)
  printf "  %-26s %s\n" "${model}" "${line:-<no chain found, check ${log}>}"
done
