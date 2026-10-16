#!/bin/bash
# =============================================================================
# run_trkeff_plots.sh - plots of a TrkEffAnalysis output + validation of the plots
#
# Usage: analysis/run_trkeff_plots.sh TRKEFF_ROOT OUTDIR LABEL [SKIMLIST] [PDFNAME]
#   TRKEFF_ROOT  TrkEffAnalysis output (e.g. TrkEff_baseline.root); its text summary
#                <name>_summary.txt must sit next to it (written by run_trkeff.sh)
#   OUTDIR       plot directory (PDF, figures/, TrackEfficiency_plots.root, numbers)
#   LABEL        sample label printed on every plot
#   SKIMLIST     optional: skim files -> adds the low-pT diagnostics page
#   PDFNAME      default TrackEfficiency.pdf
# The input file is opened read-only. Validation: OUTDIR/validation_plots.log
# (must end with VALIDATION PASSED).
# =============================================================================
IN="${1:?usage: run_trkeff_plots.sh TRKEFF_ROOT OUTDIR LABEL [SKIMLIST] [PDFNAME]}"
OUTDIR="${2:?usage: run_trkeff_plots.sh TRKEFF_ROOT OUTDIR LABEL [SKIMLIST] [PDFNAME]}"
LABEL="${3:?usage: run_trkeff_plots.sh TRKEFF_ROOT OUTDIR LABEL [SKIMLIST] [PDFNAME]}"
SKIMLIST="${4:-}"
PDF="${5:-TrackEfficiency.pdf}"
MACROS="$(cd "$(dirname "$0")" && pwd)/macros"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
IN="$(abs "${IN}")"
SUMMARY="${IN%.root}_summary.txt"
[[ -f "${SUMMARY}" ]] || { echo "ERROR: ${SUMMARY} not found" >&2; exit 2; }
[[ -z "${SKIMLIST}" ]] || SKIMLIST="$(abs "${SKIMLIST}")"
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
cd "${MACROS}"
root -l -b -q "PlotTrkEff.C+(\"${IN}\",\"${OUTDIR}\",\"${LABEL}\",\"${PDF}\",\"${SKIMLIST}\")" > "${OUTDIR}/PlotTrkEff_root.log" 2>&1
tail -1 "${OUTDIR}/PlotTrkEff_root.log"
if [[ -n "${SKIMLIST}" ]]; then
  root -l -b -q "ValidateTrkEffPlots.C+(\"${IN}\",\"${SUMMARY}\",\"${OUTDIR}\",\"${OUTDIR}/validation_plots.log\")" \
      > "${OUTDIR}/validation_plots_root.out" 2>&1
  tail -1 "${OUTDIR}/validation_plots.log"
else
  echo "(no skim list: low-pT page and its validation skipped; run ValidateTrkEffPlots.C after adding it)"
fi
