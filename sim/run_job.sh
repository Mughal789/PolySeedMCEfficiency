#!/bin/bash
# =============================================================================
# run_job.sh - run Fun4All_MBonly_polyseed.C for one job
#
# Usage: run_job.sh NEVENTS JOBID OUTDIR RELEASE [CLOCKFIX]
#   RELEASE is mandatory and must be an explicit build (e.g. new.1), never
#   "new" (the link moves). Release names cycle new.1 ... new.20, so the
#   job log identifies the build by its commits and records whether it
#   contains the polyseed fix (coresoftware PR #4473: PHGarfield.h has
#   GetNonDistortedMCT0Tbins - the same check the macro uses).
#   Local builds (e.g. a modified evaluator) are added only by setting
#   LOCAL_INSTALL=/path/to/install before calling.
#   CLOCKFIX: 0/1 (default 0) -> macro switch manualClockFix (only used on
#   releases without the fix, e.g. new.20 of 2026-10-13).
#   Environment RANDOM_SEED=1 -> macro switch randomSeed (random seeds as in Peters'
#   production, every seed printed in the job log); default 0 = fixed seed 100000 + JOBID.
# The job log records the release, its coresoftware/macros commits and
# LD_LIBRARY_PATH, so it is always known which libraries were used.
# =============================================================================
NEVENTS="${1:?usage: run_job.sh NEVENTS JOBID OUTDIR RELEASE [CLOCKFIX]}"
JOBID="${2:?}"
OUTDIR="${3:?}"
RELEASE="${4:?RELEASE missing (explicit build, e.g. new.1)}"
CLOCKFIX="${5:-0}"
b() { [[ "$1" == 1 ]] && echo true || echo false; }
MACRODIR="$(cd "$(dirname "$0")" && pwd)"
RELDIR=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/${RELEASE}

[[ "${RELEASE}" != "new" ]] || { echo "ERROR: give an explicit release (e.g. new.1), not \"new\"" >&2; exit 2; }
[[ -d "${RELDIR}" ]] || { echo "ERROR: release ${RELEASE} not found at ${RELDIR}" >&2; exit 2; }
if grep -q GetNonDistortedMCT0Tbins "${RELDIR}/include/phgarfield/PHGarfield.h" 2>/dev/null; then
  POLYFIX="yes (shared PHGarfield, clock fix built in)"
else
  POLYFIX="no (legacy polyseed configuration)"
fi
# Environment first: the sPHENIX setup script is not compatible with "set -u".
source /opt/sphenix/core/bin/sphenix_setup.sh -n "${RELEASE}" >/dev/null 2>&1
if [[ -n "${LOCAL_INSTALL:-}" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "${LOCAL_INSTALL}" >/dev/null 2>&1
fi
set -euo pipefail
mkdir -p "${OUTDIR}/logs"
LOG="${OUTDIR}/logs/job_$(printf %06d "${JOBID}").log"
{
  echo "host=$(hostname) start=$(date '+%F %T') nevents=${NEVENTS} job=${JOBID}"
  echo "release=${RELEASE} OFFLINE_MAIN=${OFFLINE_MAIN}"
  grep "tag:" "${RELDIR}/rebuild.info" || true
  echo "polyseed fix (PR #4473) in release: ${POLYFIX}"
  echo "LOCAL_INSTALL=${LOCAL_INSTALL:-none} CLOCKFIX=${CLOCKFIX} RANDOM_SEED=${RANDOM_SEED:-0}"
  echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
} > "${LOG}"
cd "${MACRODIR}"
set +e
root -l -b -q "Fun4All_MBonly_polyseed.C(${NEVENTS},${JOBID},\"${OUTDIR}\",$(b "${CLOCKFIX}"),$(b "${RANDOM_SEED:-0}"))" >> "${LOG}" 2>&1
RC=$?
set -e
echo "end=$(date '+%F %T') rc=${RC}" >> "${LOG}"
exit ${RC}
