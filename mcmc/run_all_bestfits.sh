#!/usr/bin/env bash
#
# Run find_bestfit.py for every crams fragmentation cross-section model, one
# after another. Each model writes its own bestfit_<model>.ini and a .log.
#
# Edit the configuration below to change the fit settings.
#
set -euo pipefail

cd "$(dirname "$0")"

# ── Configuration ─────────────────────────────────────────────────────────────
HALOSIZE=5
NM_MAXFEV=300     # quick Nelder-Mead pre-fit, used to seed MINUIT
MAXFEV=3000       # MINUIT (iminuit) evaluations
# ──────────────────────────────────────────────────────────────────────────────

MODELS=(
  evoli2026w93
  evoli2026st99
)

OUTDIR="bestfits"
mkdir -p "$OUTDIR"

for model in "${MODELS[@]}"; do
  tag="bestfit_${model}_h${HALOSIZE}"
  seed="$OUTDIR/${tag}_nm.ini"   # Nelder-Mead pre-fit (seed for MINUIT)
  out="$OUTDIR/${tag}.ini"       # final MINUIT best-fit
  log="$OUTDIR/${tag}.log"
  echo "=============================================================="
  echo ">>> Fitting with fragmentation model: ${model}"
  echo "    halosize=${HALOSIZE}  nm_maxfev=${NM_MAXFEV}  maxfev=${MAXFEV}"
  echo "    output: ${out}"
  echo "    log:    ${log}"
  echo "=============================================================="

  echo "--- stage 1/2: Nelder-Mead pre-fit -> ${seed}"
  python3 find_bestfit.py \
    --fragmentation-model "${model}" \
    --method nelder-mead \
    --unbounded \
    --halosize "${HALOSIZE}" \
    --maxfev "${NM_MAXFEV}" \
    --output "${seed}" \
    2>&1 | tee "${log}"

  echo "--- stage 2/2: MINUIT from the Nelder-Mead seed -> ${out}"
  python3 find_bestfit.py \
    --fragmentation-model "${model}" \
    --method iminuit \
    --unbounded \
    --halosize "${HALOSIZE}" \
    --maxfev "${MAXFEV}" \
    --start "${seed}" \
    --output "${out}" \
    2>&1 | tee -a "${log}"
done

echo
echo "=============================================================="
echo "Summary (final MINUIT chi^2 per model, h=${HALOSIZE}):"
for model in "${MODELS[@]}"; do
  log="$OUTDIR/bestfit_${model}_h${HALOSIZE}.log"
  # find_bestfit prints e.g.  "Total chi^2 = 123.4   dof = ..."
  # tail -1 picks the MINUIT (stage 2) value over the Nelder-Mead one.
  line=$(grep -E "^Total chi\^2" "${log}" | tail -1 || true)
  printf "  %-26s %s\n" "${model}" "${line:-<no chi^2 found, check ${log}>}"
done
