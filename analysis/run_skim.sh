#!/bin/bash
# =============================================================================
# run_skim.sh - skim a list of SvtxEvaluator files (one skim file per job)
#
# Usage: analysis/run_skim.sh EVALLIST OUTDIR [MAXEVENTS]
#   EVALLIST   text file with evaluator files (absolute paths)
#   OUTDIR     writes OUTDIR/skim_NNNNNN.root, OUTDIR/skim.log, OUTDIR/skim_files.list
#   MAXEVENTS  optional: first N events per file (tests)
# Needs analysis/bin/EvalSkim (analysis/build.sh).
# =============================================================================
LIST="${1:?usage: run_skim.sh EVALLIST OUTDIR [MAXEVENTS]}"
OUTDIR="${2:?usage: run_skim.sh EVALLIST OUTDIR [MAXEVENTS]}"
MAXEV="${3:-0}"
ANA="$(cd "$(dirname "$0")" && pwd)"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
[[ -x "${ANA}/bin/EvalSkim" ]] || { echo "ERROR: run ${ANA}/build.sh first" >&2; exit 2; }
mkdir -p "${OUTDIR}"
cp "${LIST}" "${OUTDIR}/evaluator_files.list"
rc=0
"${ANA}/bin/EvalSkim" -c "${ANA}/config/skim.cfg" -l "${LIST}" -o "${OUTDIR}" -n "${MAXEV}" > "${OUTDIR}/skim.log" 2>&1 || rc=$?
ls "${OUTDIR}"/skim_*.root > "${OUTDIR}/skim_files.list"
grep -E "^\[total\]|structural errors" "${OUTDIR}/skim.log"
exit ${rc}
