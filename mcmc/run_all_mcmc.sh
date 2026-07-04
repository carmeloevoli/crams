#!/usr/bin/env bash
#
# Run run_mcmc.py for every crams fragmentation cross-section model, one after
# another, starting each chain from the per-model best-fit seed produced by
#   python find_bestfit.py --scenario "$SCENARIO" --fragmentation-model <model> ...
# and saved as bestfits/bestfit_<model>_ecrs.ini (these already contain the
# fitted h and fudge_be7/9/10, so no --halosize is passed: h is a free parameter
# seeded from the .ini).
#
# Edit the configuration below to change the MCMC settings.
#
set -euo pipefail

cd "$(dirname "$0")"

# ── Configuration ─────────────────────────────────────────────────────────────
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

SCENARIO="variable_h_variable_xsecs_preliminary_be"

BESTFITDIR="bestfits"   # per-model ECRS best-fit seeds (bestfit_<model>_ecrs.ini)
OUTDIR="mcmc_chains"
mkdir -p "$OUTDIR"

for model in "${MODELS[@]}"; do
  run_tag="${model}_${SCENARIO}_ecrs"
  start="$BESTFITDIR/bestfit_${model}_ecrs.ini"   # ECRS best-fit seed (find_bestfit)
  out="$OUTDIR/mcmc_${run_tag}.npz"               # MCMC chain
  log="$OUTDIR/mcmc_${run_tag}.log"
  echo "=============================================================="
  echo ">>> MCMC with fragmentation model: ${model}"
  echo "    scenario=${SCENARIO}  nwalkers=${NWALKERS}  nburn=${NBURN}  nsteps=${NSTEPS}  ncores=${NCORES}"
  echo "    start:  ${start}"
  echo "    output: ${out}"
  echo "    log:    ${log}"
  echo "=============================================================="

  if [[ ! -f "${start}" ]]; then
    echo "!!! missing best-fit seed ${start} — run find_bestfit.py --scenario ${SCENARIO} first; skipping ${model}" >&2
    continue
  fi

  python3 run_mcmc.py \
    --scenario "${SCENARIO}" \
    --fragmentation-model "${model}" \
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
echo "Summary (mean acceptance fraction per model, scenario=${SCENARIO}):"
for model in "${MODELS[@]}"; do
  log="$OUTDIR/mcmc_${model}_${SCENARIO}_ecrs.log"
  line=$(grep -E "^Mean acceptance fraction" "${log}" 2>/dev/null | tail -1 || true)
  printf "  %-26s %s\n" "${model}" "${line:-<no chain found, check ${log}>}"
done
