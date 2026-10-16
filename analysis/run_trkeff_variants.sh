#!/bin/bash
# =============================================================================
# run_trkeff_variants.sh - tracking efficiency for all selection variants
#
# Runs TrkEffAnalysis once per variant on the same skim (no re-simulation),
# changing only these keys of the baseline configuration (via -s key=value):
#   trk_nmvtx_min         0, 2, 3      (TPC >= 18 in all variants)
#   use_reco_zvtx         0 (ztrue), 1 (zreco)
#   truth_min_tpc_layers  0 (incl = all accepted G4 primaries),
#                         ${TRACKABLE_LAYERS:-18} (track = trackable primaries)
# Output: OUTDIR/TrkEff_mvtx<M>_<ztrue|zreco>_<incl|track>.root + _summary.txt + .log
#
# Usage: analysis/run_trkeff_variants.sh SKIMLIST OUTDIR [BASECFG] [EVALLIST [STEP2_BASELINE_ROOT]]
# Comparison plots (OUTDIR/plots, CompareTrkEffVariants.C) are always made; the
# sample label printed on them is taken from the environment variable LABEL.
#   EVALLIST  optional: validate every variant against the original evaluator
#             (ValidateTrkEff.C, including a split-input + hadd test), then
#             run the cross-variant checks (ValidateTrkEffVariants.C; needs
#             STEP2_BASELINE_ROOT = TrkEff_baseline.root from run_trkeff.sh, for the
#             "variant mvtx0_ztrue_incl reproduces the baseline" check)
# =============================================================================
SKIMLIST="${1:?usage: run_trkeff_variants.sh SKIMLIST OUTDIR [BASECFG] [EVALLIST [STEP2_BASELINE_ROOT]]}"
OUTDIR="${2:?usage: run_trkeff_variants.sh SKIMLIST OUTDIR [BASECFG] [EVALLIST [STEP2_BASELINE_ROOT]]}"
ANA="$(cd "$(dirname "$0")" && pwd)"
BASECFG="${3:-${ANA}/config/trkeff_baseline.cfg}"
EVALLIST="${4:-}"
STEP2="${5:-}"
LAYERS="${TRACKABLE_LAYERS:-18}"
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
[[ -x "${ANA}/bin/TrkEffAnalysis" ]] || { echo "ERROR: run ${ANA}/build.sh first" >&2; exit 2; }
mkdir -p "${OUTDIR}/split"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
SKIMLIST="$(abs "${SKIMLIST}")"
[[ -z "${EVALLIST}" ]] || EVALLIST="$(abs "${EVALLIST}")"
[[ -z "${STEP2}" ]] || STEP2="$(abs "${STEP2}")"
# split input lists (first / second half of the skim files) for the hadd test
NF=$(grep -c . "${SKIMLIST}")
head -n $((NF / 2)) "${SKIMLIST}" > "${OUTDIR}/split/listA"
tail -n +$((NF / 2 + 1)) "${SKIMLIST}" > "${OUTDIR}/split/listB"

for M in 0 2 3; do
  for Z in ztrue zreco; do
    for D in incl track; do
      TAG="mvtx${M}_${Z}_${D}"
      OPTS=(-s "trk_nmvtx_min=${M}" -s "use_reco_zvtx=$([[ ${Z} == zreco ]] && echo 1 || echo 0)"
            -s "truth_min_tpc_layers=$([[ ${D} == track ]] && echo "${LAYERS}" || echo 0)")
      OUT="${OUTDIR}/TrkEff_${TAG}"
      "${ANA}/bin/TrkEffAnalysis" -c "${BASECFG}" "${OPTS[@]}" -l "${SKIMLIST}" -o "${OUT}.root" -t "${OUT}_summary.txt" > "${OUT}.log" 2>&1
      MSG="${TAG}: $(grep 'selected events' "${OUT}.log")"
      if [[ -n "${EVALLIST}" ]]; then
        for P in A B; do
          "${ANA}/bin/TrkEffAnalysis" -c "${BASECFG}" "${OPTS[@]}" -l "${OUTDIR}/split/list${P}" -o "${OUTDIR}/split/TrkEff_${TAG}_${P}.root" \
              -t "${OUTDIR}/split/TrkEff_${TAG}_${P}_summary.txt" > "${OUTDIR}/split/TrkEff_${TAG}_${P}.log" 2>&1
        done
        hadd -f "${OUTDIR}/split/TrkEff_${TAG}_A+B.root" "${OUTDIR}/split/TrkEff_${TAG}_A.root" "${OUTDIR}/split/TrkEff_${TAG}_B.root" > /dev/null
        ( cd "${ANA}/macros" && root -l -b -q "ValidateTrkEff.C+(\"${EVALLIST}\",\"${OUT}.root\",\"${OUTDIR}/validation_${TAG}.log\",\"${OUT}_summary.txt\",\"${OUTDIR}/split/TrkEff_${TAG}_A+B.root\")" \
            > "${OUTDIR}/validation_${TAG}_root.out" 2>&1 )
        MSG="${MSG} | $(tail -1 "${OUTDIR}/validation_${TAG}.log")"
      fi
      echo "${MSG}"
    done
  done
done
mkdir -p "${OUTDIR}/plots"
( cd "${ANA}/macros" && root -l -b -q "CompareTrkEffVariants.C+(\"${OUTDIR}\",\"${OUTDIR}/plots\",\"${LABEL:-}\",\"TrkEffVariants.pdf\")" \
    > "${OUTDIR}/plots/Compare_root.log" 2>&1 )
tail -1 "${OUTDIR}/plots/Compare_root.log"
if [[ -n "${EVALLIST}" && -n "${STEP2}" ]]; then
  ( cd "${ANA}/macros" && root -l -b -q "ValidateTrkEffVariants.C+(\"${OUTDIR}\",\"${STEP2}\",\"${OUTDIR}/validation_variants.log\",\"${OUTDIR}/plots\")" \
      > "${OUTDIR}/validation_variants_root.out" 2>&1 )
  echo "cross-variant checks: $(tail -1 "${OUTDIR}/validation_variants.log")"
fi
