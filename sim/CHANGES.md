# Fun4All_MBonly_polyseed.C — changes relative to M. Peters' macro

Reference: `/sphenix/user/mjpeters/analysis/LightFlavorRatios/geometric_acceptance/simulation/Fun4All_HFG_MBonly_polyseed.C`
(version of 2026-10-12 17:30; read only). Target release: **first release containing coresoftware PR #4473 (3e4855d97) and macros
PR #1393** — the rebuilt new.1 of 2026-10-14 (coresoftware bb9d98e, macros f18832e).
Release names cycle (new.1 … new.20); always identify a build by its commits.

## Kept exactly as Peters
- PYTHIA8, Detroit tune: `pythia8_MB_Detroit_Tony.cfg` copied byte-identical
  (PDF:pSet 17, MPI bProfile 2 / ecmRef 200 / pT0Ref 1.40 / ecmPow 0.135 /
  coreRadius 0.56 / coreFraction 0.78, CR range 5.4, SoftQCD:inelastic,
  decays limited to a 1e-6 mm cylinder → K0S/Λ decay in GEANT4).
- p+p beam: `Input::ApplysPHENIXBeamParameter` (vertex (−0.07, 0.141) cm, σz 16 cm, 1.5 mrad).
- Silicon seeding and final vertex: release `Tracking_Reco_SiliconSeed_run2pp()`, `Tracking_Reco_Vertex_run2pp()` (no beam numbers in the macro; since 2026-10-15).
- Detectors: MBDFAKE, PIPE + absorber, MVTX, INTT, TPC, MICROMEGAS; `TRACKING::streaming_mode = true`.
- Run number 29; MVTX misalignment ON in G4 and Acts (open decision).
- Silicon seeding chain, silicon-only vertex finder (`SiliconSvtxVertexMap`, same cuts).
- Polyseed modules and order: module tracks → assembled tracks → crossing finder
  (triggered mode, silicon vertex map) → poly clusterizer → poly tracks → poly vertexer
  → seed converter → cluster converter; then matching and track fit (`_run2pp`).
- `Global_Reco()`, `Tracking_Eval()` with its default settings (baseline evaluator).

## Changed
| # | Peters | Ours | Why |
|---|---|---|---|
| 1 | Local `"G4_Input.C"` | Release `<G4_Input.C>` | His copy is identical to new.20 |
| 2 | Random seed | Fixed seed `100000 + jobId` (`RANDOMSEED`) **by default**; since 2026-10-15 the macro argument `randomSeed = true` (job script: `RANDOM_SEED=1`) restores Peters' behaviour: no `RANDOMSEED` flag, every seed drawn by `PHRandomSeed` and printed in the job log (`ntp_info.seed` = 0) | Reproducible by default (K1, 10k); random seeds for the 1M production on request |
| 3 | CDB tag `newcdbtag` (data tag) | **`MDC2`** (official MC tag; changed 2026-10-14) | Experts: no data cdbtag for simulations. newcdbtag's TRACKINGALIGNMENT = 2025 data alignment, displaces reco silicon in MC by 0.3-0.6 mm (`../audit_mvtx/`); MDC2's is all zero. Test `../mdc2_test/` (K1, same seeds): eps_vtx(abs(z)<8) 0.654 → 0.735, identical to newcdbtag + zero alignment (`../zeroalign_test/`). Other effective differences (no visible effect): gas falls back to `files/gasfiles` (same v_drift 7.551 cm/us at 387.5 V/cm), INTT digitisation mask `INTT_Hotmap_Sim`, MVTX dead pixels 30 instead of 244, `is_data_reco = false` (only the release beam spot: seeding at (0,0) for MC, #6). The 2026-10-13 MDC2 collapse (T500_new20) was on the pre-fix new.20. Previous version: `../macros_record/2026-10-14_before_MDC2/` |
| 4 | `sampa_tzero_bias = −65 ns` | `0` | Reco-only bias; poly-cluster converter applies it asymmetrically (≈1.5 cm on side 0) |
| 5 | `TPC_Clustering()` | removed | Converter deletes standard TPC clusters anyway; official chain omits it |
| 6 | `Tracking_Reco_SiliconSeed_run2pp()` (release) | same release function (**reverted 2026-10-15**; 2026-10-14 version set the simulated beam spot (−0.07, 0.141) by hand) | No beam numbers of our own. With MDC2 the release uses seeding beam spot (0,0) for MC; Peters' `newcdbtag` made the release skip that line, so his seeding used the module default (−0.072, 0.141). Test `../beamspot_release_test/`: track content changed in 1 of 1000 events |
| 7 | Modules: `setT0(0)`, `setCMVoltageDefault(387.5)`, `setIsNonDistortedMC(true)`, `setUseSurveyGeometry(false)` | none on modules; one shared `PHGarfield` with `SetIsNonDistortedMC(true)`, `SetUseSurveyGeometry(true)` | PR #4473 / PR #1393 API; gives t0 = 1.5, ADC clock from geometry (53.326 ns), 387.5 V/cm; `setT0(0)` would override t0 |
| 8 | `DSTClusterPruning` | removed | Not in the official chain |
| 9 | Second loose `PHSimpleVertexFinder` (production) | one finder only, release run2pp cuts | Clean, documented vertexing |
| 10 | `Tracking_Reco_Vertex_run2pp()` (release; box ±0.2 cm around (0,0)) | same release function (**reverted 2026-10-15**; 2026-10-14 version centred the box on (−0.07, 0.141)) | Release default. Test `../beamspot_release_test/` (K1, same seeds): eps_vtx(abs(z)<8) 0.735 → 0.730 (291 → 289 of 396; the lost vertices are 2-track vertices with reco y ≥ 0.18 cm, at the box edge y = 0.2), tracking and MVTX association unchanged |
| 11 | Extra `build_truthreco_tables()` | removed | Already inside `Tracking_Eval()` |
| 12 | DST strips G4 hits, hit sets, truth containers, Acts trajectories | strips only G4HIT_PIPE, G4HIT_SVTXSUPPORT | Evaluator can be re-run on the DST for the fix studies |
| 13 | KFParticle V0 reconstruction (`HF_selections.C`), DecayFinder includes | removed | Not needed for this study |
| 14 | Hard-coded output in Peters' tg01 area | `outDir` argument | Write only into our area |
| 15 | — | `SaveGitTags()`; log records release, commits, LD_LIBRARY_PATH | Provenance |
| 16 | Polyseed field maps from the CDB (data tag) | The same four files given explicitly to `PHGarfield` (`SetElectricFieldMap3D`, `SetFrameElectricFieldMap3D`) | Not available under MDC2 (first test aborted); scaled to zero in non-distorted MC |

## Switch `manualClockFix` (releases without the fix only, e.g. `run_job.sh ... new.20 1`)
Applies the PR #4473 TPC timing fix by hand to the crossing finder and clusterizer:
`setTpcAdcClock(53.326184)` (simulated clock; default 56.881 = data clock), `setT0(1.5)`
(instead of Peters' 0), crossing finder `setStartZ(±102.605)` (clusterizer default already).
Ignored on releases with the fix, which apply it themselves. Removed switches (2026-10-13):
`petersSampaBias` (−65 ns gave only +3 points; 0 is consistent once the clock is right)
and `petersCdbTag` (the tag is fixed: newcdbtag until 2026-10-14, MDC2 since).

## Release-dependent block (releases without PR #4473, e.g. new.20 of 2026-10-13)
new.20 has no shared-PHGarfield design: the clusterizer builds a private PHGarfield
(distortions ON by default). When running on new.20 the macro therefore configures the
crossing finder and clusterizer **exactly as Peters** (`setT0(0)`, `setCMVoltageDefault(387.5)`,
`setIsNonDistortedMC(true)`, clusterizer `setUseSurveyGeometry(false)`) plus the four field maps.
The ADC clock stays at 56.881 ns: a new.20 run is the **pre-fix reference**. On releases with the fix
the block is skipped. The log line `polyseed configuration = ...` states which mode ran.
The mode is chosen by API, not by release name: a release has the fix if
`$OFFLINE_MAIN/include/phgarfield/PHGarfield.h` contains `GetNonDistortedMCT0Tbins`
(`MBPOLY::releaseHasPolyseedFix()`; `run_job.sh` writes the same check to the job log).

## Open points checked in the smoke test (10 events) [historical]
- MDC2 tag: all TPC/polyseed CDB payloads resolve (Peters used a data tag).
- No `TpcConditionsReco` is registered (Peters' MC did not use it; the data macros do).
  If PHGarfield needs it in MC, the smoke test will abort and it is added.
- Log must show from both polyseed modules: `non-distorted MC: t0=1.5 tbins, ADC clock=53.3262 ns`,
  and exactly one `PHSimpleVertexFinder` registration besides `SiliconVertexFinder`.

Previous version (hard-coded beam spot, 2026-10-14): `../macros_record/2026-10-15_before_release_beamspot/`.
