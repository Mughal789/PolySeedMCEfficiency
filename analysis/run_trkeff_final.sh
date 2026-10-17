#!/bin/bash
# =============================================================================
# run_trkeff_final.sh - FINAL tracking-efficiency selection (MVTX >= 2) on one
# sample: analysis, independent validation, step checks and per-species maps
#
# Usage: analysis/run_trkeff_final.sh SKIMLIST OUTDIR VARIANTSDIR LABEL [EVALLIST]
#   SKIMLIST     skim files of the sample
#   OUTDIR       writes (all TrkEffAnalysis v4, config/trkeff.cfg):
#                  TrkEff_final.root            FINAL (TPC 20, MVTX 2, INTT 1, DCA 1.5, reco vertex, no truth eta) = the
#                                               data track selection (Collect + Corral); the correction maps
#                  TrkEff_final_mvtx3.root      the same with MVTX >= 3 (variant / systematic, not a correction map)
#                  TrkEff_step{0,1,2,3}.root    one-cut steps from MVTX 2 + TPC 18 to FINAL
#                  TrkEff_final_samedenom.root  FINAL cuts, truth |eta| < 1.1 and the Step-2 binning (for the table)
#                  validation_*.log, final_validation.log, maps/ (TrkEffFinalMaps.root/.pdf, figures/)
#   VARIANTSDIR  the Step-4A variants of the same sample (validation table)
#   LABEL        sample label on the plots
#   EVALLIST     optional evaluator files: independent validation of FINAL and the MVTX3 variant (ValidateTrkEff.C)
# Environment: TRUTH_ETA_CUT=1 runs every output with use_truth_eta_cut = 1 (|eta_true| < 1.1,
#              truth-eta mode); default 0 = the FINAL definition (no truth-eta cut). Use a
#              separate OUTDIR per mode.
# =============================================================================
LIST="${1:?usage: run_trkeff_final.sh SKIMLIST OUTDIR VARIANTSDIR LABEL [EVALLIST]}"
OUTDIR="${2:?}"; VARDIR="${3:?}"; LABEL="${4:?}"; EVALLIST="${5:-}"
ANA="$(cd "$(dirname "$0")" && pwd)"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n "${RELEASE:-new.1}" >/dev/null 2>&1
set -euo pipefail
mkdir -p "${OUTDIR}"; OUTDIR="$(cd "${OUTDIR}" && pwd)"; VARDIR="$(cd "${VARDIR}" && pwd)"
LIST="$(cd "$(dirname "${LIST}")" && pwd)/$(basename "${LIST}")"
CFG="${ANA}/config/trkeff.cfg"
MODE=(-s "use_truth_eta_cut=${TRUTH_ETA_CUT:-0}")   # truth-eta mode of every output of this run
BASECFG="${ANA}/config/trkeff_baseline.cfg"
YB=$(grep '^y_edges' "${BASECFG}" | cut -d= -f2 | tr -d ' '); PB=$(grep '^pt_edges' "${BASECFG}" | cut -d= -f2 | tr -d ' ')
run() { # name [overrides...]
  local name=$1; shift
  "${ANA}/bin/TrkEffAnalysis" -c "${CFG}" "${MODE[@]}" "$@" -l "${LIST}" -o "${OUTDIR}/${name}.root" -t "${OUTDIR}/${name}_summary.txt" \
      > "${OUTDIR}/${name}.log" 2>&1 || { echo "ERROR: ${name} (see ${name}.log)"; exit 3; }
  echo "${name}: $(grep -o 'selected true particles.*' "${OUTDIR}/${name}.log")"
}
run TrkEff_final
run TrkEff_final_mvtx3 -s trk_nmvtx_min=3
run TrkEff_step0 -s trk_ntpc_min=18 -s trk_nintt_min=0 -s trk_dca_xy_max=0 -s trk_dca_z_max=0
run TrkEff_step1 -s trk_nintt_min=0 -s trk_dca_xy_max=0 -s trk_dca_z_max=0
run TrkEff_step2 -s trk_dca_xy_max=0 -s trk_dca_z_max=0
run TrkEff_step3
run TrkEff_final_samedenom -s use_truth_eta_cut=1 -s "y_edges=${YB}" -s "pt_edges=${PB}"  # truth |eta| < 1.1 (later -s wins)
# absolute path before changing directory
[[ -z "${EVALLIST}" ]] || EVALLIST="$(cd "$(dirname "${EVALLIST}")" && pwd)/$(basename "${EVALLIST}")"
cd "${ANA}/macros"
if [[ -n "${EVALLIST}" ]]; then
  for n in TrkEff_final TrkEff_final_mvtx3; do
    root -l -b -q "ValidateTrkEff.C+(\"${EVALLIST}\",\"${OUTDIR}/${n}.root\",\"${OUTDIR}/validation_${n}.log\",\"${OUTDIR}/${n}_summary.txt\")" > "${OUTDIR}/validation_${n}_root.out" 2>&1 || true
    echo "validation ${n}: $(tail -1 "${OUTDIR}/validation_${n}.log")"
  done
fi
root -l -b -q "ValidateFinalTrkEff.C+(\"${OUTDIR}\",\"${VARDIR}\",\"${OUTDIR}/final_validation.log\")" > "${OUTDIR}/final_validation_root.out" 2>&1 || true
echo "final checks: $(tail -1 "${OUTDIR}/final_validation.log")"
root -l -b -q "MakeFinalTrkEffMaps.C+(\"${OUTDIR}/TrkEff_final.root\",\"${OUTDIR}/maps\",\"${LABEL}\",\"${OUTDIR}/TrkEff_final_mvtx3.root\",\"MVTX#geq3\")" > "${OUTDIR}/maps_root.out" 2>&1 || true
echo "maps: $(tail -2 "${OUTDIR}/maps/maps_validation.log" | head -1)"
