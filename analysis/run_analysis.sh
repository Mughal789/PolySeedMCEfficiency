#!/bin/bash
# =============================================================================
# run_analysis.sh - full analysis of one sample (one production, or any list of
# SvtxEvaluator files)
#
# Usage: analysis/run_analysis.sh SAMPLEDIR [LABEL] [EVALLIST]
#   SAMPLEDIR  production directory from sim/submit.sh (evaluator/outputeval_*.root);
#              results go to SAMPLEDIR/analysis/
#   LABEL      text for plot titles (default: basename of SAMPLEDIR)
#   EVALLIST   optional list of evaluator files (absolute paths) to use instead of
#              SAMPLEDIR/evaluator/ - e.g. somebody else's production (read only);
#              SAMPLEDIR is then only the output location
#
# Writes SAMPLEDIR/analysis/:
#   evaluator_files.list
#   qa/                 PlotEvalQA: vertex, track, truth-matching QA (EvalQA_hists.root)
#   z{all,8,10}/        VtxEffAnalysis for good tracks TPC-only (VtxEff_TPC.root) and
#                       >= 3 MVTX clusters (VtxEff_MVTX3.root, the analysis definition),
#                       validation_*.log, plots_*/, tables_*.txt,
#                       breakdown/ (VtxEffBreakdown: efficiencies, no-vertex reasons, VtxEffBreakdown_summary.txt)
#   audit/              ResidPhi.txt (MVTX truth-reco distance vs phi),
#                       MvtxAudit.txt, LostTracks.txt/.pdf
#   SUMMARY.txt         main numbers
# Good tracks: pT > 0.1 GeV/c, |eta| < 1.1, >= 18 TPC clusters (+ >= 3 MVTX for MVTX3).
# Needs analysis/bin/VtxEffAnalysis (analysis/build.sh). ROOT macros are compiled
# with ACLiC next to the sources (the package directory must be writable).
# =============================================================================
SAMPLE_IN="${1:?usage: run_analysis.sh SAMPLEDIR [LABEL] [EVALLIST]}"
mkdir -p "${SAMPLE_IN}"
SAMPLE="$(cd "${SAMPLE_IN}" && pwd)"
LABEL="${2:-$(basename "${SAMPLE}")}"
EVALLIST="${3:-}"
ANA="$(cd "$(dirname "$0")" && pwd)"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
[[ -x "${ANA}/bin/VtxEffAnalysis" ]] || { echo "ERROR: run ${ANA}/build.sh first" >&2; exit 2; }
OUT="${SAMPLE}/analysis"
mkdir -p "${OUT}/qa" "${OUT}/audit"
LIST="${OUT}/evaluator_files.list"
if [[ -n "${EVALLIST}" ]]; then
  cp "${EVALLIST}" "${LIST}"
else
  ls "${SAMPLE}"/evaluator/outputeval_*.root > "${LIST}"
fi
echo "${LABEL}: $(wc -l < "${LIST}") evaluator file(s) -> ${OUT}"
cd "${ANA}/macros"

# 1. evaluator QA
root -l -b -q "PlotEvalQA.C+(\"${LIST}\",\"${OUT}/qa\",\"${LABEL}\")" > "${OUT}/qa/PlotEvalQA_root.log" 2>&1
echo "qa: $(tail -1 "${OUT}/qa/EvalQA_summary.txt" | cut -c1-120)"

# 2. vertex efficiency, three z windows, two good-track selections + VtxEffBreakdown
for Z in all 8 10; do
  ZD="${OUT}/z${Z}"
  ZC=$([[ ${Z} == all ]] && echo 1000 || echo ${Z})
  mkdir -p "${ZD}/breakdown"
  for SEL in MVTX3 TPC; do
    F="${ZD}/VtxEff_${SEL}.root"
    "${ANA}/bin/VtxEffAnalysis" -c "${ANA}/config/vtxeff_${SEL}_z${Z}.cfg" -l "${LIST}" -o "${F}" > "${ZD}/analysis_${SEL}.log" 2>&1
    root -l -b -q "ValidateVtxEff.C+(\"${F}\",\"${ZD}/validation_${SEL}.log\")" > "${ZD}/validation_root_${SEL}.out" 2>&1
    mkdir -p "${ZD}/plots_${SEL}"
    root -l -b -q "PlotVtxEff.C+(\"${F}\",\"${ZD}/plots_${SEL}\",\"${ZD}/tables_${SEL}.txt\")" > "${ZD}/plot_${SEL}.log" 2>&1
    echo "z${Z} ${SEL}: $(tail -1 "${ZD}/validation_${SEL}.log")"
  done
  root -l -b -q "VtxEffBreakdown.C+(\"${ZD}/VtxEff_TPC.root\",\"${LIST}\",\"${ZD}/breakdown\",\"${LABEL}\",${ZC},3)" > "${ZD}/breakdown/VtxEffBreakdown_root.log" 2>&1
done

# 3. MVTX audit: truth-reco cluster distance vs phi, association, lost tracks
root -l -b -q "ResidPhi.C+(\"${LIST}\")" 2>&1 | grep "median" > "${OUT}/audit/ResidPhi.txt"
root -l -b -q "MvtxAudit.C+(\"${LIST}\",\"${LABEL}\")" > "${OUT}/audit/MvtxAudit.txt" 2>&1
root -l -b -q "LostTracks.C+(\"${LIST}\",\"${OUT}/audit/LostTracks.pdf\",\"${LABEL}\")" > "${OUT}/audit/LostTracks.txt" 2>&1

# 4. summary
{
  echo "sample: ${LABEL}"
  echo "files : $(wc -l < "${LIST}") (${LIST})"
  grep -E "^  all collisions +:" "${OUT}/zall/breakdown/VtxEffBreakdown_summary.txt" | sed "s/^ */eps_vtx /"
  grep -E "^  \|z_true\| < 8 cm +:" "${OUT}/z8/breakdown/VtxEffBreakdown_summary.txt" | sed "s/^ */eps_vtx /"
  echo "MVTX truth-reco distance (median per 20 deg phi bin, um): $(awk '{print $(NF-2)}' "${OUT}/audit/ResidPhi.txt" | sort -n | sed -n '1p;$p' | paste -sd'-')"
} | tee "${OUT}/SUMMARY.txt"
