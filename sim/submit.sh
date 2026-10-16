#!/bin/bash
# =============================================================================
# submit.sh - write (and submit) the Condor job file for one production
#
# Usage: sim/submit.sh OUTDIR NJOBS NEVENTS RELEASE [--dry-run]
#   OUTDIR   production directory (created); jobs write OUTDIR/{DST,evaluator,logs}
#   NJOBS    number of jobs, JOBID = 0 .. NJOBS-1, random seed = 100000 + JOBID
#            (the same JOBID always gives the same events)
#   NEVENTS  events per job (50 events take ~2-3 h)
#   RELEASE  explicit build, e.g. new.1 (never "new"); must contain coresoftware
#            PR #4473 (checked by run_job.sh and written to every job log)
#   --dry-run  only write OUTDIR/condor.job, do not submit
# The Condor log goes to /tmp (site policy); job stdout/stderr go to OUTDIR/logs.
# =============================================================================
OUTDIR="${1:?usage: submit.sh OUTDIR NJOBS NEVENTS RELEASE [--dry-run]}"
NJOBS="${2:?}"
NEVENTS="${3:?}"
RELEASE="${4:?RELEASE missing (explicit build, e.g. new.1)}"
DRY="${5:-}"
set -euo pipefail
SIMDIR="$(cd "$(dirname "$0")" && pwd)"
[[ "${RELEASE}" != "new" ]] || { echo "ERROR: give an explicit release (e.g. new.1), not \"new\"" >&2; exit 2; }
mkdir -p "${OUTDIR}/logs"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
cat > "${OUTDIR}/condor.job" <<EOF
# written by ${SIMDIR}/submit.sh on $(date '+%F %T')
# ${NJOBS} jobs x ${NEVENTS} events, release ${RELEASE}, seeds 100000 + JOBID (JOBID 0..$((NJOBS - 1)))
Universe           = vanilla
Executable         = ${SIMDIR}/run_job.sh
should_transfer_files = NO
request_cpus       = 1
request_memory     = 6000MB
PeriodicHold       = (JobStatus == 2) && ((time() - EnteredCurrentStatus) > 6*3600)
Priority           = 20
STEPDIR            = ${OUTDIR}
Log                = /tmp/condor-vtxeffpoly-\$(Cluster).log
Arguments          = ${NEVENTS} \$(Process) \$(STEPDIR) ${RELEASE}
Output             = \$(STEPDIR)/logs/condor-\$(Process).out
Error              = \$(STEPDIR)/logs/condor-\$(Process).err
Queue ${NJOBS}
EOF
echo "wrote ${OUTDIR}/condor.job"
if [[ "${DRY}" == "--dry-run" ]]; then
  exit 0
fi
condor_submit "${OUTDIR}/condor.job"
