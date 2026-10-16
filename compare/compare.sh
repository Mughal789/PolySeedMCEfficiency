#!/bin/bash
# =============================================================================
# compare.sh - compare samples analysed with analysis/run_analysis.sh
#
# Usage: compare/compare.sh OUTDIR "SAMPLEDIR_A,SAMPLEDIR_B,..." "label A,label B,..."
#   SAMPLEDIR_x  the directory given to run_analysis.sh (its analysis/ subdirectory is read)
#   up to 8 samples, drawn in the given order (first page = epsilon_vtx per sample)
# Writes OUTDIR/:
#   CompareVtxEff.pdf, png/, CompareVtxEff_table.txt   physics comparison (CompareVtxEff.C)
#   qa/CompareSteps.pdf, png/, CompareSteps_table.txt  every evaluator-QA histogram overlaid
#                                                      (CompareSteps.C)
# Labels must not contain commas.
# =============================================================================
OUTDIR="${1:?usage: compare.sh OUTDIR \"dirA,dirB\" \"labelA,labelB\"}"
DIRS="${2:?}"
LABELS="${3:?}"
CMP="$(cd "$(dirname "$0")" && pwd)"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
mkdir -p "${OUTDIR}/qa"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
ANADIRS=""
QAFILES=""
IFS=',' read -ra DS <<< "${DIRS}"
for D in "${DS[@]}"; do
  A="$(cd "${D}/analysis" && pwd)" || { echo "ERROR: ${D}/analysis missing - run analysis/run_analysis.sh first" >&2; exit 2; }
  ANADIRS+="${ANADIRS:+,}${A}"
  QAFILES+="${QAFILES:+,}${A}/qa/EvalQA_hists.root"
done
cd "${CMP}"
root -l -b -q "CompareVtxEff.C+(\"${ANADIRS}\",\"${LABELS}\",\"${OUTDIR}\")" > "${OUTDIR}/CompareVtxEff_root.log" 2>&1
grep -E "^wrote|WARNING" "${OUTDIR}/CompareVtxEff_root.log" || true
root -l -b -q "CompareSteps.C+(\"${QAFILES}\",\"${LABELS}\",\"${OUTDIR}/qa\")" > "${OUTDIR}/qa/CompareSteps_root.log" 2>&1
echo "wrote ${OUTDIR}/qa/CompareSteps.pdf"
cat "${OUTDIR}/CompareVtxEff_table.txt"
