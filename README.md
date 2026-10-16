# PolySeedMCEfficiency

Vertex-finding efficiency and single-particle tracking efficiency in sPHENIX p+p minimum-bias simulation with TPC polyseed tracking.

The module covers the whole chain:

1. simulate p+p collisions at √s = 200 GeV;
2. reconstruct them with the polyseed chain, in a configuration validated for MC;
3. reduce the evaluator output to compact skims;
4. measure:
   - **ε_vtx**, the probability that a collision gets a reconstructed vertex;
   - **ε_trk(y, pT)**, the probability that a generated primary π⁺, π⁻, K⁺, K⁻, p or p̄ is reconstructed and selected.

   ε_trk comes per species and charge, as maps for efficiency corrections in two-particle correlation analyses.

Every result is checked by an **independent recomputation** from the original evaluator files.

---

## Contents

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Repository layout](#3-repository-layout)
4. [Quick test (2 jobs × 100 events)](#4-quick-test-2-jobs--100-events)
5. [Full workflow](#5-full-workflow)
6. [Simulation and reconstruction configuration](#6-simulation-and-reconstruction-configuration)
7. [Definitions](#7-definitions)
8. [Outputs](#8-outputs)
9. [Validation](#9-validation)
10. [Reference results (1M events)](#10-reference-results-1m-events)
11. [Known limitations](#11-known-limitations)

---

## 1. Overview

```
 sim/submit.sh ─► Condor jobs (sim/run_job.sh ─► sim/Fun4All_MBonly_polyseed.C)
                  PYTHIA8 ─► GEANT4 ─► polyseed reconstruction ─► DST + SvtxEvaluator ntuples
                                                                          │
              ┌───────────────────────────────────────────────────────────┤
              ▼                                                           ▼
 analysis/run_analysis.sh                                     analysis/run_skim.sh (EvalSkim)
   VtxEffAnalysis: vertex efficiency                            compact skim, ~2 kB/event
   (+ QA, efficiency breakdown, MVTX audit)                               │
                                                                          ▼
                                       analysis/run_trkeff.sh, run_trkeff_variants.sh,
                                       run_trkeff_plots.sh, run_trkeff_final.sh  (TrkEffAnalysis)
                                       single-particle efficiency, maps, Corral-format histograms
```

The three programs are plain C++ with ROOT (`analysis/src/`). Plotting and validation are ROOT macros (`analysis/macros/`) that the run scripts compile with ACLiC.

Two analyses share the evaluator files:

| analysis | input | program |
|---|---|---|
| vertex efficiency | evaluator ntuples | `VtxEffAnalysis` |
| tracking efficiency | skims | `TrkEffAnalysis` |

`EvalSkim` keeps the information both need, so the slow evaluator files are read only once for the tracking efficiency.

---

## 2. Requirements

- **An sPHENIX environment.** The scripts source `/opt/sphenix/core/bin/sphenix_setup.sh -n <release>` themselves.
- **A release containing coresoftware PR #4473** (the polyseed MC timing fix), e.g. `new.1`.
  - Always give an explicit build: never `new`.
  - Every job log records the build's commits and whether the fix is present.
- **HTCondor** for productions.
  - The Condor log goes to `/tmp` (site policy).
  - Job output goes to the production directory.
- **ROOT** (from the release) for the analysis programs and macros.
- **A writable copy of this repository.** ROOT compiles the macros next to their sources.

No path outside the repository is hard-coded. The external inputs are:
- the release;
- the CDB on `/cvmfs`;
- optionally, evaluator files of other productions (read only).

---

## 3. Repository layout

```
sim/
  Fun4All_MBonly_polyseed.C     simulation + reconstruction macro (G4, polyseed, vertexing, evaluator, DST)
  pythia8_MB_Detroit_Tony.cfg   PYTHIA8 settings (SoftQCD:inelastic, Detroit tune)
  run_job.sh                    one job: environment, provenance log, macro
  submit.sh                     Condor submission (NJOBS x NEVENTS)
  check_jobs.sh                 per-job status and configuration checks
  CHANGES.md                    every change relative to the original production macro, with reasons
analysis/
  build.sh                      builds bin/VtxEffAnalysis, bin/EvalSkim, bin/TrkEffAnalysis
  src/VtxEffAnalysis.cc         vertex efficiency, multiplicity response, track/vertex QA
  src/EvalSkim.cc               evaluator -> skim (EventTree, TruthTree, RecoTree)
  src/TrkEffAnalysis.cc         single-particle efficiency (all cuts from a config file)
  src/EvalCommon.h              shared evaluator helpers
  config/                       vtxeff_{TPC,MVTX3}_z{all,8,10}.cfg, skim.cfg,
                                trkeff_baseline.cfg, trkeff_final_mvtx3.cfg
  run_analysis.sh               vertex efficiency of one sample (QA, 6 configurations, breakdown, MVTX audit)
  run_skim.sh                   skim a list of evaluator files
  run_trkeff.sh                 tracking efficiency of one configuration (+ independent validation)
  run_trkeff_variants.sh        12 selection variants (MVTX 0/2/3 x true/reco z x inclusive/trackable)
  run_trkeff_plots.sh           plot set of one TrkEffAnalysis output (+ plot validation)
  run_trkeff_final.sh           FINAL selection, MVTX >= 2 reference, one-cut steps, checks, maps
  macros/                       plotting, validation and map macros (see §8, §9)
  TRKEFF.md                     detailed documentation of the tracking-efficiency part
compare/
  compare.sh, CompareVtxEff.C, CompareSteps.C   side-by-side comparison of analysed samples
```

---

## 4. Quick test (2 jobs × 100 events)

The test runs the full chain on 2 jobs of 100 events in about 40 minutes (most of it is the simulation; the jobs run in parallel). Keep outputs outside the repository.

**Why two jobs and 200 events:**
- The variant validation splits the file list in two and checks that merging the halves with `hadd` reproduces the full run, so it needs at least 2 files.
- Every species needs at least one generated particle, otherwise its efficiency is undefined (0/0) and the plot validation reports it. p̄ is the rarest: about 3 per 100 events in the acceptance.

```bash
M=/path/to/PolySeedMCEfficiency
T=/path/to/test10                          # any writable directory outside the repository
export LABEL="test 200 ev"                  # sample label printed on the plots

$M/analysis/build.sh                                         # once
$M/sim/submit.sh $T/S200ev 2 100 new.1                       # 2 jobs x 100 events (JOBID 0,1; seeds 100000,100001)
# ... wait for the job (condor_q) ...
$M/sim/check_jobs.sh $T/S200ev                               # must report 2 OK, 0 failed

ls $T/S200ev/evaluator/outputeval_*.root > $T/eval.list
$M/analysis/run_skim.sh  $T/eval.list $T/skim                # skim + structural checks
$M/analysis/run_analysis.sh $T/S200ev "$LABEL"            # vertex efficiency, 6 validations
$M/analysis/run_trkeff.sh $T/skim/skim_files.list $T/trkeff \
      $M/analysis/config/trkeff_baseline.cfg $T/eval.list    # tracking efficiency + independent validation
$M/analysis/run_trkeff_variants.sh $T/skim/skim_files.list $T/trkeff/variants \
      $M/analysis/config/trkeff_baseline.cfg $T/eval.list $T/trkeff/TrkEff_baseline.root
$M/analysis/run_trkeff_plots.sh $T/trkeff/TrkEff_baseline.root $T/trkeff/plots "$LABEL" $T/skim/skim_files.list
$M/analysis/run_trkeff_final.sh $T/skim/skim_files.list $T/trkeff/final_mvtx3 $T/trkeff/variants "$LABEL" $T/eval.list
```

**Pass criterion:** every validation line printed by the scripts ends in `VALIDATION PASSED`, with one possible exception: **check C of the final checks** (`final checks:` line of `run_trkeff_final.sh`).
- Check C is a physics consistency check, not a code check. It compares the MVTX ≥ 2 → MVTX ≥ 3 survival of pions between true-z and reco-z event selection; they must agree within 0.005.
- The true difference is about −0.002. On a 200-event sample the statistical scatter is about ±0.004, so check C can fail by fluctuation alone.
- It passes reliably from about 1k events on. K1 (1k), 10k and 1M all pass.
- All other checks must pass on any sample size.

Fixed seeds make the test reproducible: the same release and macro always give the same events.

---

## 5. Full workflow

### 5.1 Build

```bash
analysis/build.sh
```

Builds three executables in `analysis/bin/`. Needs only ROOT.

### 5.2 Production

```bash
sim/submit.sh OUTDIR NJOBS NEVENTS RELEASE [--dry-run]
sim/check_jobs.sh OUTDIR
```

- **Seeds:** JOBID runs from 0 to NJOBS−1, with seed 100000 + JOBID. The same JOBID always gives the same events, so two configurations can be compared event by event.
- **Speed:** about 10–20 s per event; 1000 events per job take 2–5 h. Memory ≤ 4 GB (6 GB requested).
- **Size:** DST and evaluator files take about 0.8 MB/event each.
- **Random seeds** (independent large samples, as in an official production):
  - write the job file with `--dry-run`;
  - add the line `environment = "RANDOM_SEED=1"` to `OUTDIR/condor.job`;
  - submit with `condor_submit`.

  `PHRandomSeed` then draws every seed from `std::random_device` and prints it in the job log (`PHRandomSeed::GetSeed() seed:`). `ntp_info.seed` is NaN in this mode.
- **Use a different JOBID range per sample** (edit `Queue` in `condor.job`) so file names never collide.

`check_jobs.sh` checks every job log for:
- exit code 0;
- CDB tag MDC2;
- zero tracking alignment;
- a release with the polyseed fix;
- MC TPC timing (t0 = 1.5 time bins, ADC clock 53.3262 ns);
- no data reconstruction mode.

It also checks that the macro uses the release beam-spot defaults.

Interactive single job: `sim/run_job.sh NEVENTS JOBID OUTDIR new.1`. The log is `OUTDIR/logs/job_NNNNNN.log`.

### 5.3 Vertex efficiency

```bash
analysis/run_analysis.sh SAMPLEDIR [LABEL] [EVALLIST]
```

- **SAMPLEDIR** is a production directory. Results go to `SAMPLEDIR/analysis/`.
- With **EVALLIST**, any list of evaluator files is analysed instead, e.g. somebody else's production (read only).
- Calling `analysis/bin/VtxEffAnalysis -c CONFIG -l LIST -o OUT.root` directly also works:
  - relative names in `LIST` are resolved against `input_dir` of the config, by default the evaluator directory of the 1M reference production;
  - the program never writes inside `input_dir`.

**Steps:**
1. evaluator QA (`PlotEvalQA`);
2. `VtxEffAnalysis` for three z windows (all, |z_true| < 10 cm, |z_true| < 8 cm) and two good-track definitions (TPC-only, and ≥ 3 MVTX clusters);
3. after each configuration, `ValidateVtxEff` and the plots (`PlotVtxEff`);
4. `VtxEffBreakdown` per z window: efficiency per multiplicity group, why collisions have no vertex, why tracks are missing;
5. the MVTX audit (`ResidPhi`, `MvtxAudit`, `LostTracks`);
6. `SUMMARY.txt`.

Each step reads the full evaluator set. For large samples, the steps can run as parallel jobs: the steps are independent, except `VtxEffBreakdown`, which needs the `VtxEff_TPC.root` of its z window.

### 5.4 Skim

```bash
analysis/run_skim.sh EVALLIST OUTDIR [MAXEVENTS]
```

Writes one `skim_NNNNNN.root` per job, plus `skim.log` (must report 0 structural errors) and `skim_files.list`.

| tree | content |
|---|---|
| `EventTree` | true and reconstructed vertex, multiplicities |
| `TruthTree` | charged G4 primaries (and V0s) with true pT, η, y (PDG mass) and truth hit counts |
| `RecoTree` | every track with cluster counts, DCA, matched truth particle and charge check |

**Independent check:** `ValidateSkim.C` re-reads the evaluator files with its own code and compares job by job.

```bash
cd analysis/macros
root -l -b -q 'ValidateSkim.C+("eval.list","SKIMDIR","report.log",NEXP,-1,-1,"")'
```

### 5.5 Tracking efficiency

```bash
analysis/run_trkeff.sh          SKIMLIST OUTDIR [CONFIG] [EVALLIST]
analysis/run_trkeff_variants.sh SKIMLIST OUTDIR [BASECFG] [EVALLIST [BASELINE_ROOT]]
analysis/run_trkeff_plots.sh    TRKEFF_ROOT OUTDIR LABEL [SKIMLIST] [PDFNAME]
analysis/run_trkeff_final.sh    SKIMLIST OUTDIR VARIANTSDIR LABEL [EVALLIST]     # TRUTH_ETA_CUT=1 for the truth-eta mode
```

| script | what it produces |
|---|---|
| `run_trkeff.sh` | one configuration. Baseline: `config/trkeff_baseline.cfg` (\|z_true\| < 8 cm, TPC ≥ 18, no MVTX requirement, truth \|η\| < 1.1) |
| `run_trkeff_variants.sh` | 12 variants of the baseline: MVTX ≥ 0/2/3 × true/reco vertex z × all/trackable particles, plus comparison plots |
| `run_trkeff_plots.sh` | 12-page plot set: maps, ε vs pT and y, correct charge, counts, low-pT diagnostics |
| `run_trkeff_final.sh` | **FINAL** selection (§7) with an MVTX ≥ 2 reference, the one-cut steps from the baseline to FINAL, consistency checks and the per-species maps in Corral format |

Run the final chain twice: once as it is (no truth-η cut, the FINAL definition) and once with `TRUTH_ETA_CUT=1`. Then compare the two modes with:
- `macros/ValidateTruthEtaModes.C` (particle-level check);
- `macros/CompareTruthEtaModes.C` (PDF).

**With EVALLIST**, every output is validated against the evaluator files (§9). For large samples, run the validation on a subset of jobs: analyse the subset's skims and give the matching evaluator list.

All keys, definitions and outputs: [`analysis/TRKEFF.md`](analysis/TRKEFF.md).

### 5.6 Comparing samples

```bash
compare/compare.sh OUTDIR "SAMPLEDIR1,SAMPLEDIR2,..." "label1,label2,..."
```

Writes ε_vtx and the QA of up to 8 analysed samples side by side (`CompareVtxEff.pdf`, `qa/CompareSteps.pdf`).

`macros/CompareTrkEffSamples.C` compares the tracking efficiency of two independent samples bin by bin (pulls).

---

## 6. Simulation and reconstruction configuration

| item | setting | why |
|---|---|---|
| release | explicit build with coresoftware **PR #4473**, e.g. `new.1` | correct non-distorted MC TPC timing (ADC clock 53.326 ns, t0 = 1.5 time bins). Without it the TPC drift is stretched by 6.7 % and the vertex efficiency is about halved |
| CDB tag | **`MDC2`** (official simulation tag), timestamp 6, run 29 | the data tag loads 2025 data alignment, which the simulated detector does not have: MVTX clusters 0.3–0.6 mm off, ε_vtx 0.65 instead of 0.73 |
| tracking alignment | all zero (MDC2 payload) | the simulation has no residual misalignment |
| geometry | survey geometry in G4 and reconstruction | the simulated detector is built with it |
| polyseed | one shared `PHGarfield` with `SetIsNonDistortedMC(true)` (387.5 V/cm, no space charge); the four field maps are given explicitly | MDC2 has no polyseed payloads |
| generator | PYTHIA8 `SoftQCD:inelastic`, Detroit tune; one collision per event; no trigger emulation | |
| decays | PYTHIA decays only particles with cτ < 1 nm; K⁰_S, Λ, Σ±, Ξ, Ω, π±, K± decay in GEANT4 | |
| beam (simulation) | release `Input::ApplysPHENIXBeamParameter`: vertex (−0.07, 0.141) cm, σ_z ≈ 16 cm | official p+p beam parameters |
| beam spot (reconstruction) | release defaults (`Tracking_Reco_SiliconSeed_run2pp`, `Tracking_Reco_Vertex_run2pp`); no beam numbers in the macro | official reconstruction; costs −0.005 in ε_vtx (\|z\| < 8 cm) compared with the true beam spot |
| vertex finder | release `Tracking_Reco_Vertex_run2pp`: `PHSimpleVertexFinder` on fitted tracks, ≥ 3 MVTX clusters per track, pT ≥ 0.1, pair DCA < 0.05 cm, χ²/ndf ≤ 300 | the same function as in the data production |

Every change relative to the original production macro, with its reason: [`sim/CHANGES.md`](sim/CHANGES.md).

---

## 7. Definitions

### Vertex efficiency

ε_vtx = (collisions with ≥ 1 vertex in `SvtxVertexMap`) / (generated collisions)

- Every generated inelastic collision is in the denominator, also those without a charged particle in the acceptance.
- Windows on the true z: all, |z_true| < 10 cm, |z_true| < 8 cm (the analysis cut). The MVTX acceptance ends at |z| ≈ 20 cm.
- Uncertainties are Clopper–Pearson 68.3 % intervals.
- A **good track** has pT > 0.1 GeV/c, |η| < 1.1 and ≥ 18 TPC clusters; for the MVTX3 definition also ≥ 3 MVTX clusters.
- Collisions without a vertex are classified by the first failing condition:

  | reason | condition |
  |---|---|
  | R1 | \|z_true\| > 20 cm |
  | R2 | no good track |
  | R3 | no good track with ≥ 3 MVTX clusters |
  | R4 | exactly one such track |
  | R5 | ≥ 2 such tracks, but no vertex (the finder) |

### Single-particle tracking efficiency

ε_h(y, pT) = Reconstructed_h(y, pT) / Generated_h(y, pT), for h = π⁺, π⁻, K⁺, K⁻, p, p̄.

- **Generated:** G4 primaries of species h in the selected events, binned in **true** y (PDG mass) and **true** pT, with pT > 0.1 GeV/c. Weak-decay daughters, conversions and secondaries are not primaries.
- **Reconstructed:** generated particles with ≥ 1 selected track matched to them.
  - The matched particle is the one that contributed most clusters to the track.
  - Each particle is counted once.
  - ε_Q counts only tracks with the correct charge.
- **FINAL selection** (`config/trkeff_final_mvtx3.cfg`), matched to the data analysis:
  - **event:** a vertex from the final finder with |z_reco| < 8 cm;
  - **track:** TPC ≥ 20, MVTX ≥ 3, INTT ≥ 1, |DCA_xy| and |DCA_z| ≤ 1.5 cm to the reconstructed vertex (straight line from the track state);
  - **denominator:** no truth-η cut;
  - **analysis windows:**

    | species | \|y\| | pT (GeV/c) |
    |---|---|---|
    | π | < 1.0 | ≥ 0.1 |
    | K | < 0.7 | 0.1–0.6 |
    | p, p̄ | < 0.6 | 0.4–1.0 |

  The vertex efficiency and the dE/dx PID efficiency are separate quantities and are not included.

Raw, unweighted counts are stored everywhere, so samples can be merged exactly with `hadd`. Bins without generated particles are undefined: they get no point and are hatched in the maps, never set to 0.

---

## 8. Outputs

**Production:**
- `DST/MBonly_polyseed_DST_NNNNNN.root`;
- `evaluator/outputeval_NNNNNN.root`;
- `logs/job_NNNNNN.log`: release, commits, CDB tag, every setting, seeds.

**Vertex efficiency** (`SAMPLEDIR/analysis/`):

| path | content |
|---|---|
| `SUMMARY.txt` | ε_vtx for all z and \|z_true\| < 8 cm; MVTX truth–reco distance |
| `z8/breakdown/VtxEffBreakdown_summary.txt`, `VtxEffBreakdown.pdf` | **main result**: ε_vtx per multiplicity group and why collisions have no vertex |
| `z{all,8,10}/VtxEff_{MVTX3,TPC}.root` | histograms and the per-event tree `evt`, with `validation_*.log`, `plots_*/` and `tables_*.txt` next to them |
| `qa/EvalQA.pdf` | vertex residuals, track quantities, truth matching |
| `audit/` | MVTX truth–reco cluster distance vs φ, MVTX association, lost tracks |

**Tracking efficiency:**

| path | content |
|---|---|
| `TrkEff_<config>.root` | per species: 2D raw counts (`hDen_<sp>_ypt`, `hNum_`, `hNumQ_`), pT and y projections, `TEfficiency` objects |
| `_summary.txt` | integrated numbers |
| `final_mvtx3/maps/TrkEffFinalMaps.{root,pdf}` | **correction maps.** `heffpty_<sp>` is a TH2D (x = y, y = pT) holding ε, or −1 where nothing was generated; errors in `_errlo`, `_errhi` |

---

## 9. Validation

Each analysis has an independent check, written with separate code. Every report ends in `VALIDATION PASSED` or `VALIDATION FAILED`.

| macro | what it checks |
|---|---|
| `ValidateVtxEff.C` | ~1000 internal checks per configuration: event keys, totals, response matrices, numerator ≤ denominator, tree vs histograms, no NaN, unique objects |
| `ValidateSkim.C` | the skim against the evaluator, job by job: event counts, vertices, keys, PDG and charge, rapidity, reco rows, matching, truth retention, NaN; random-seed mode |
| `ValidateTrkEff.C` | Generated, Reconstructed and correct-charge counts recomputed from the evaluator, bin by bin, with its own rapidity formula and its own INTT and DCA code; also NQ ≤ N ≤ D, projections, Clopper–Pearson intervals and split input + `hadd` = full run |
| `ValidateTrkEffVariants.C` | the 12 variants are consistent with each other and with the baseline |
| `ValidateFinalTrkEff.C` | the one-cut steps (Generated unchanged, Reconstructed never increases), the validation table, the vertex-conditioning check |
| `ValidateTruthEtaModes.C` | the truth-η switch removes exactly the particles with \|η_true\| ≥ 1.1, cell by cell |
| `ValidateTrkEffPlots.C` | every plotted point equals the stored counts |
| `VertexProvenance.C` | the skim's event vertex equals the DST `SvtxVertexMap`, event by event. It needs the sPHENIX libraries |

Edge cases that only large samples show are handled explicitly and printed as information:
- **random-seed mode:** the seed is NaN;
- **float rounding at bin edges:** y is rounded to float, as in the skim;
- **vertex ties:** two vertices with the same number of tracks; the first in evaluator order is used (0.04 % of vertexed events);
- **ROOT autosave cycles** of large trees;
- **check C (vertex conditioning)** uses an absolute tolerance of 0.005, because the true effect (−0.002) is significant with 10⁶ events.

---

## 10. Reference results (1M events)

The reference sample has 1000 jobs × 1000 events, release `new.1`, MDC2 and random seeds. Every validation above passed.

| quantity | value |
|---|---|
| ε_vtx, all collisions | 0.470 (469,882 / 1,000,000) |
| ε_vtx, \|z_true\| < 8 cm | **0.730 ± 0.001** (280,080 / 383,541) |
| no vertex at \|z_true\| < 8 cm | 56 % no good track, 21 % no MVTX-matched track, 19 % one MVTX-matched track, 4 % the finder |

**Single-particle efficiency** (statistical errors 0.0005 for π, 0.002 for K, 0.004 for p):

| selection | π⁺ | π⁻ | K⁺ | K⁻ | p | p̄ |
|---|---|---|---|---|---|---|
| baseline (TPC ≥ 18, no MVTX, true z) | 0.868 | 0.873 | 0.772 | 0.758 | 0.828 | 0.800 |
| **FINAL** (analysis windows) | **0.679** | **0.748** | **0.470** | **0.543** | **0.733** | **0.775** |
| FINAL cuts with MVTX ≥ 2 (reference) | 0.718 | 0.783 | 0.494 | 0.568 | 0.771 | 0.812 |

**Main observations:**

- **The vertex efficiency is limited by multiplicity, not by the finder.** Collisions with ≥ 2 good MVTX-matched tracks get a vertex 98.6 % of the time.
- **An MVTX requirement makes every positive species less efficient than its antiparticle.** For pions the difference is 0.069 ± 0.001.
  - **Cause:** the default silicon–TPC matching window, −0.25 < Δφ < +0.05 rad for both charges, acts on a charge-odd Δφ residual from the TPC-seed curvature.
  - On top of this there is a smaller geometric loss of particles crossing only two MVTX layers.
  - The data production uses the same default windows, so the effect is expected in data.
- **Low-pT edge:** below 0.15 GeV/c the efficiency is limited by the polyseed converter acceptance (poly-tracks need > 20 clusters and fitted pT > 0.1 GeV/c).
- **Without an MVTX requirement**, the charge differences are small but resolved at 10⁶ events: π 0.005, K 0.014, p 0.027.

---

## 11. Known limitations

- **The MC is an idealised detector:** no residual misalignment, no TPC distortions or space charge, no TPC/TPOT dead maps. MC efficiencies are therefore a best case relative to data.
- **No trigger emulation and no pileup.** ε_vtx is per generated inelastic collision.
- **MDC2 has no polyseed payloads.** They are set explicitly in the macro until the experts add them to the MC tag.
- **The polyseed MC chain is still developing upstream.** After a release change, re-validate with `check_jobs.sh` and the reference numbers above.
- **The DCA** is a straight-line approximation from the track state to the reconstructed vertex, not the helix DCA.
