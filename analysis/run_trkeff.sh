#!/bin/bash
# =============================================================================
# run_trkeff.sh - single-particle tracking efficiency from EvalSkim files,
# optionally followed by the independent validation against the evaluator.
#
# Usage: analysis/run_trkeff.sh SKIMLIST OUTDIR [CONFIG] [EVALLIST]
#   SKIMLIST  text file with skim_NNNNNN.root files (one per line)
#   OUTDIR    writes OUTDIR/TrkEff_<config>.root, _summary.txt, .log
#   CONFIG    TrkEffAnalysis configuration (default analysis/config/trkeff_baseline.cfg)
#   EVALLIST  optional: the evaluator files the skims were made from; if given,
#             macros/ValidateTrkEff.C recomputes everything from them
#             (OUTDIR/validation_<config>.log, must end with VALIDATION PASSED)
# Needs analysis/bin/TrkEffAnalysis (analysis/build.sh).
# =============================================================================
SKIMLIST="${1:?usage: run_trkeff.sh SKIMLIST OUTDIR [CONFIG] [EVALLIST]}"
OUTDIR="${2:?usage: run_trkeff.sh SKIMLIST OUTDIR [CONFIG] [EVALLIST]}"
ANA="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${3:-${ANA}/config/trkeff_baseline.cfg}"
EVALLIST="${4:-}"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
[[ -x "${ANA}/bin/TrkEffAnalysis" ]] || { echo "ERROR: run ${ANA}/build.sh first" >&2; exit 2; }
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
TAG="$(basename "${CONFIG}" .cfg)"; TAG="${TAG#trkeff_}"
OUT="${OUTDIR}/TrkEff_${TAG}"
"${ANA}/bin/TrkEffAnalysis" -c "${CONFIG}" -l "${SKIMLIST}" -o "${OUT}.root" -t "${OUT}_summary.txt" > "${OUT}.log" 2>&1
grep -E "selected events|errors" "${OUT}.log"
if [[ -n "${EVALLIST}" ]]; then
  EVALLIST="$(cd "$(dirname "${EVALLIST}")" && pwd)/$(basename "${EVALLIST}")"
  ( cd "${ANA}/macros" && root -l -b -q "ValidateTrkEff.C+(\"${EVALLIST}\",\"${OUT}.root\",\"${OUTDIR}/validation_${TAG}.log\",\"${OUT}_summary.txt\",\"\")" \
      > "${OUTDIR}/validation_${TAG}_root.out" 2>&1 )
  tail -1 "${OUTDIR}/validation_${TAG}.log"
fi
