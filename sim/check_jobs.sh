#!/bin/bash
# =============================================================================
# check_jobs.sh - status and conditions of every job of one production
#
# Usage: sim/check_jobs.sh OUTDIR
# Per job log (OUTDIR/logs/job_NNNNNN.log) it checks:
#   finished with rc=0, CDB tag MDC2, zero tracking alignment (MDC2 payload
#   8c7a..._tracking_alignment_parameters_empty.txt), release with the polyseed
#   fix, non-distorted MC timing (t0 = 1.5 time bins, ADC clock 53.3262 ns),
#   and that no data-only condition was requested ("Setting reconstruction for data").
# Prints one line per failing job and a summary; exit code 1 if any check fails.
# =============================================================================
OUTDIR="${1:?usage: check_jobs.sh OUTDIR}"
shopt -s nullglob
LOGS=("${OUTDIR}"/logs/job_*.log)
[[ ${#LOGS[@]} -gt 0 ]] || { echo "no job logs in ${OUTDIR}/logs"; exit 1; }
NOK=0
NBAD=0
NRUN=0
for L in "${LOGS[@]}"; do
  J=$(basename "${L}" .log)
  END=$(grep -m1 "^end=" "${L}" || true)
  if [[ -z "${END}" ]]; then
    NRUN=$((NRUN + 1))
    continue
  fi
  P=()
  [[ "${END}" == *"rc=0" ]] || P+=("${END##* }")
  grep -q " CDB MDC2 " "${L}" || P+=("tag!=MDC2")
  grep -q "8c7a4bd674c9ed3fd85fc5a9e216cf1d_tracking_alignment_parameters_empty" "${L}" || P+=("alignment!=zero")
  grep -q "polyseed fix (PR #4473) in release: yes" "${L}" || P+=("release-without-fix")
  grep -q "non-distorted MC: t0=1.5 tbins, ADC clock=53.3262 ns" "${L}" || P+=("TPC-timing")
  grep -q "Setting reconstruction for data" "${L}" && P+=("data-reco-mode")
  if [[ ${#P[@]} -eq 0 ]]; then
    NOK=$((NOK + 1))
  else
    NBAD=$((NBAD + 1))
    echo "${J}: ${P[*]}"
  fi
done
# beam spot (since 2026-10-15): reconstruction uses the RELEASE defaults - the macro calls
# Tracking_Reco_SiliconSeed_run2pp() and Tracking_Reco_Vertex_run2pp() and sets no beam
# numbers of its own (with MDC2: seeding beam spot (0,0), vertex box |x|,|y| < 0.2 cm around (0,0)).
# The simulated collisions use Input::ApplysPHENIXBeamParameter (release G4_Input.C).
MACRO="$(cd "$(dirname "$0")" && pwd)/Fun4All_MBonly_polyseed.C"
CODE=$(sed 's#//.*##' "${MACRO}")   # code without comments
if grep -q "Tracking_Reco_SiliconSeed_run2pp();" <<< "${CODE}" && grep -q "Tracking_Reco_Vertex_run2pp();" <<< "${CODE}" \
   && ! grep -q -E "set_beamSpotXY|setBeamSpotCut|beam_x|beam_y" <<< "${CODE}"; then
  echo "beam spot: release defaults (Tracking_Reco_SiliconSeed_run2pp + Tracking_Reco_Vertex_run2pp, no beam numbers in the macro)"
else
  echo "beam spot: macro does NOT use the release defaults (custom beam-spot code found or release functions missing)"
  NBAD=$((NBAD + 1))
fi
NEVAL=$(ls "${OUTDIR}"/evaluator/outputeval_*.root 2>/dev/null | wc -l)
echo "${OUTDIR}: ${#LOGS[@]} logs: ${NOK} OK, ${NBAD} failed checks, ${NRUN} not finished; ${NEVAL} evaluator files"
[[ ${NBAD} -eq 0 ]]
