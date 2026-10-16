# Single-particle tracking efficiency (`TrkEffAnalysis`)

Computes the tracking efficiency of identified charged hadrons
(π⁺, π⁻, K⁺, K⁻, p, p̄) versus **true rapidity y** and **true transverse
momentum pT**. It runs on the compact skim files written by `EvalSkim` (see the
README), not on the large evaluator files.

The pipeline:

```
SvtxEvaluator files ──EvalSkim──> skim_NNNNNN.root ──TrkEffAnalysis──> TrkEff_<config>.root (+ summary.txt)
 (outputeval_*.root)               EventTree                            raw counts, projections,
                                    TruthTree                            TEfficiency objects
                                    RecoTree
```

## Build and run

Requirements: ROOT ≥ 6 (`root-config` in `PATH`). On the sPHENIX nodes, the
scripts source `sphenix_setup.sh -n new.1` if ROOT is not already set up.

```bash
analysis/build.sh        # builds bin/VtxEffAnalysis, bin/EvalSkim, bin/TrkEffAnalysis

# 1. skim the evaluator files (one skim file per job)
analysis/run_skim.sh  EVALLIST  SKIMDIR

# 2. tracking efficiency (+ optional independent validation against the evaluator)
analysis/run_trkeff.sh  SKIMDIR/skim_files.list  OUTDIR  [CONFIG]  [EVALLIST]
```

The direct call is:
```bash
analysis/bin/TrkEffAnalysis -c analysis/config/trkeff_baseline.cfg -l skim_files.list \
                            -o TrkEff_baseline.root -t TrkEff_baseline_summary.txt
```

Outputs of separate runs (for example one per production chunk) can be merged
with `hadd`. All histograms are raw counts, and `TEfficiency` objects store
passed and total counts. The merged file is therefore identical to a single
run over all inputs; this is verified by the validation.

## Definitions

### Event selection
| `use_reco_zvtx` | selected collisions |
|---|---|
| `0` (baseline) | \|z_true\| < `zvtx_max` — every simulated collision, with or without a reconstructed vertex |
| `1` | a reconstructed vertex exists **and** \|z_reco\| < `zvtx_max` (z_reco = vertex with the most tracks) |

All event multiplicities are combined. The vertex-finding efficiency
(`VtxEffAnalysis`) does not use this setting; it keeps its own denominator.

### Denominator D(y, pT): selected true particles
- A **G4 primary**, i.e. a particle handed to GEANT4 by PYTHIA. PYTHIA decays
  only states that decay within 1 nm, so weak-decay daughters (K0S, Λ, …),
  photon conversions and hadronic-interaction products are G4 secondaries and
  never enter.
- PDG code of a configured species (`species`).
- In a selected event.
- pT_true > `truth_pt_min`, and |η_true| < `truth_abseta_max`. Set
  `truth_abseta_max = 0` to disable the η cut; the y binning then defines the
  acceptance.
- Binned in **true y**, computed with the PDG mass (y = asinh(pz / mT)), and
  **true pT**.

### Numerators
- **N(y, pT)**: particles of D with at least one **selected reconstructed
  track** that is **matched** to the particle.
- **NQ(y, pT)**: as N, but at least one of those tracks has the **correct
  reconstructed charge**.

Selected track (reconstructed quantities; no DCA cut):
pT > `trk_pt_min`, |η| < `trk_abseta_max`, TPC clusters ≥ `trk_ntpc_min`,
MVTX clusters ≥ `trk_nmvtx_min`.

Efficiencies: ε = N / D (any charge) and ε_Q = NQ / D (correct charge).

### Truth matching
- Each reconstructed track is matched to its **dominant truth particle**: the
  G4 particle that contributed the most clusters (SvtxEvaluator / truth-reco
  tables; there is no minimum purity).
- A true particle is identified by the event-local key **(job, event,
  gtrackID)**, because G4 track ids restart in every event.
- **Each true particle enters the numerator at most once**, however many
  selected tracks are matched to it. Extra tracks are counted as information
  ("duplicate tracks" in `hTrkEffStats`).
- Only tracks matched to the particle **itself** count. For example, a π whose
  track is matched to its decay muon is *not* found (baseline definition; no
  credit through descendants).

### Statistics
`TEfficiency` with Clopper–Pearson 68.3% intervals (the same convention as
ε_vtx). A bin with **D = 0 has an undefined efficiency**: `TEfficiency` reports
0 there, so always check the total (`GetTotalHistogram()`). The text summary
writes `undefined` for such bins.

## Configuration (`analysis/config/trkeff_baseline.cfg`)
| key | baseline | meaning |
|---|---|---|
| `use_reco_zvtx` | 0 | event selection on true (0) or reconstructed (1) z |
| `zvtx_max` | 8.0 cm | \|z\| cut |
| `truth_pt_min` | 0.1 GeV/c | true pT cut |
| `truth_abseta_max` | 1.1 | true \|η\| cut; ≤ 0 disables |
| `species` | 211, −211, 321, −321, 2212, −2212 | PDG codes (only these six are supported) |
| `trk_pt_min` | 0.1 GeV/c | reconstructed pT |
| `trk_abseta_max` | 1.1 | reconstructed \|η\| |
| `trk_ntpc_min` | 18 | TPC clusters |
| `trk_nmvtx_min` | 0 | MVTX clusters (0 = no MVTX requirement) |
| `y_edges` | −1.2 … 1.2 in steps of 0.2 | true-y bin edges |
| `pt_edges` | 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0 GeV/c | true-pT bin edges |

Unknown keys stop the program, so a typo cannot silently fall back to a
default. All values are stored in the output (`TNamed config`).

**Skim coverage check.** At run time the program reads the retention window
stored in every skim file. It refuses to run if the requested truth acceptance
is not fully inside that window, because the denominator would be incomplete.
The window is pT ≥ 0 and (|η| < 1.5 or |y| < 1.5). With the η cut disabled,
the in-range bins are complete, but under/overflow is not (see the note in the
summary).

## Output (`TrkEff_<config>.root`)
`<sp>` ∈ {`pip`, `pim`, `kp`, `km`, `p`, `pbar`}.

| object | type | content |
|---|---|---|
| `hDen_<sp>_ypt`, `hNum_<sp>_ypt`, `hNumQ_<sp>_ypt` | TH2D | raw counts; x = y_true, y = pT_true [GeV/c] |
| `hDen_<sp>_pt`, `hNum_<sp>_pt`, `hNumQ_<sp>_pt` | TH1D | pT projections (sum over the in-range y bins) |
| `hDen_<sp>_y`, `hNum_<sp>_y`, `hNumQ_<sp>_y` | TH1D | y projections (sum over the in-range pT bins) |
| `eff_<sp>_{ypt,pt,y}` | TEfficiency | N / D (any charge), Clopper–Pearson 68.3% |
| `effQ_<sp>_{ypt,pt,y}` | TEfficiency | NQ / D (correct charge) |
| `hTrkEffStats` | TH1D | counters: files, events, selected events, selected / found particles, selected tracks (matched / unmatched / on denominator particles / duplicates), particles outside the binning |
| `config`, `definitions`, `program`, `skim_files` | TNamed | provenance |

In total: 18 TH2D, 36 TH1D projections + 1 stats histogram, 36 TEfficiency,
4 TNamed. The text summary has the integrated counts per species and ε(pT),
ε_Q(pT) per species.

## Validation (`analysis/macros/ValidateTrkEff.C`)
The macro recomputes D, N and NQ **directly from the evaluator files** with its
own code: event selection from `ntp_gpoint` / `ntp_vertex`, true y from the
PDG mass with an independent formula, matching from `ntp_track.gtrackID`, and
charge from `TDatabasePDG`. It then checks:

- the selected events;
- D, N and NQ bin by bin, including under/overflow;
- no duplicate numerator entries;
- NQ ≤ N ≤ D in every bin;
- raw counts are unweighted;
- the projections equal the 2D histograms summed over the in-range bins;
- the `TEfficiency` objects equal the histograms, with the Clopper–Pearson
  option;
- the species charges;
- every empty bin is reported as `undefined`;
- optionally, split input lists merged with `hadd` give a result identical to
  the full run.

The true y is computed in double and **rounded to float** before binning, because the skim stores y as a float. Without the rounding, a particle exactly on a bin edge is binned differently. The 1M sample had one: a π⁺ with y = 0.99999998 in double and 1.0 as a float (validator fix 2026-10-16).

## Known limitations
- **Threshold migration.** The lowest bin (0.10 ≤ pT_true < 0.15 GeV/c) has a
  low efficiency (about 0.35 for pions in K1). Near the 0.1 GeV/c cut, the
  reconstructed pT and the TPC cluster count (low-pT tracks curl) decide whether
  a track passes. This is real, not a bug, but bins at the cut need care.
- **Matching has no purity requirement.** The dominant truth particle is used
  as is. Purity variables are in the RecoTree for later studies.
- **No credit through descendants** (π → μ, interactions). Losses from decays
  and material are part of the inefficiency, by definition.
- **Only G4 primaries** in the denominator. The secondary / fake contamination
  of the reconstructed sample is not part of ε; it needs a separate purity
  correction.
- **Statistics.** 1000 events give about 800 π±, 80 K± and 35 p / p̄ inside
  \|z\| < 8 cm: most (y, pT) bins are empty or have large uncertainties.
  Species-resolved 2D maps need production-scale samples.
- With the truth η cut disabled, the denominator under/overflow (y outside the
  binning) is incomplete; use the in-range bins only.

## Plots (`PlotTrkEff.C`)

```bash
analysis/run_trkeff_plots.sh  TrkEff_baseline.root  PLOTDIR  "sample label"  [skim_files.list]  [TrackEfficiency.pdf]
```
The script runs `macros/PlotTrkEff.C` and then `macros/ValidateTrkEffPlots.C`
(`PLOTDIR/validation_plots.log` must end with `VALIDATION PASSED`). The input
file is opened read-only and nothing is recomputed: every efficiency shown
comes from the stored `TEfficiency` objects.

**Inputs:**
- from the `TrkEffAnalysis` output: `eff_<sp>_{ypt,pt,y}`, `effQ_<sp>_{ypt,pt,y}`,
  `hDen/hNum/hNumQ_<sp>_{ypt,pt,y}`, `hTrkEffStats`, `config`;
- optionally, the skim list (only for the low-pT page).

**Outputs in `PLOTDIR`:**

| file | content |
|---|---|
| `<pdf>` | multipage PDF (12 pages with the skim list, 11 without) |
| `figures/pageNN_*.{pdf,png}` | every page as a single figure |
| `figures/<plot>_<species>.{pdf,png}` | every panel as a stand-alone figure, for slides |
| `TrackEfficiency_plots.root` | every plotted object (graphs `g_*`, maps `map_*`, low-pT fractions) |
| `TrackEfficiency_numbers.txt` | integrated D, N, NQ, ε, ε_Q with Clopper–Pearson intervals, empty-bin counts |
| `LowPtDiagnostics.txt` | low-pT loss categories (with a skim list) |

**Pages:**

| page | content |
|---|---|
| 1 | summary: sample, selections, definitions; table of D, N, NQ, integrated ε and ε_Q, empty bins |
| 2–4 | ε(y, pT) maps: π⁺/π⁻, K⁺/K⁻, p/p̄ |
| 5 | ε vs pT: π±, K±, p/p̄ (positive filled, negative open markers), all six species |
| 6 | ε vs y: same layout |
| 7 | any-charge vs correct-charge ε vs pT, per species |
| 8 | raw counts vs pT: generated D, reconstructed N, correct charge NQ |
| 9–11 | ε_Q(y, pT) maps (correct charge) |
| 12 | low-pT diagnostics (with a skim list), see below |

**Statistical treatment:**
- **1D:** `TEfficiency::CreateGraph`, Clopper–Pearson 68.3% asymmetric
  errors. Bins with D = 0 have **no point**: never a 0, never interpolated.
- **2D:** the colour is passed/total in bins with D > 0; a measured
  efficiency of 0 is painted with the darkest colour. Bins with D = 0 are
  covered by a **hatched grey box** ("undefined"). Every defined bin shows
  `passed/total`, so bins with one or two particles are visible as such.
- **Integrated values:** sums over the in-range bins, with the Clopper–Pearson
  68.3% interval of the sums. These are identical to the Step 2 numbers.
- No smoothing, fitting, rebinning or interpolation.

**Low-pT diagnostics (`LowPtDiagnostics.C`, page 12).** Every denominator
particle (same selections, read from the `config` of the `TrkEffAnalysis`
output) gets exactly one category from stored facts:
- found;
- a matched track exists but its best match (most clusters) fails the pT,
  |η|, n_TPC or n_MVTX cut;
- no matched track, split by the number of TPC layers with truth clusters
  (≥ N_TPC, 1 … N_TPC − 1, 0).

The categories record which reconstruction stage or cut a particle fails, not
why: no physics mechanism is assigned. The "found" counts are checked against
the Step 2 numerator in every bin.

**Validation (`ValidateTrkEffPlots.C`):**
- every plotted point equals passed/total of the Step 2 histograms;
- the y errors equal Clopper–Pearson recomputed from the counts, and the x
  errors equal the bin half-width;
- the number of points equals the number of defined bins;
- every map bin with D > 0 equals passed/total, and every map bin with D = 0
  is masked;
- the integrated numbers equal the Step 2 report, and the CP intervals are
  recomputed;
- the stored projections equal the 2D sums;
- the low-pT "found" and total counts equal the Step 2 numerator and
  denominator.

**Known limitations of the plots:**
- With 1000 events most kaon and (anti)proton bins are empty or hold one or
  two particles. The maps show this as it is.
- The positive/negative comparison is an overlay; no ratio is drawn, because a
  ratio of two low-statistics efficiencies needs its own interval treatment.
- The ε_Q(y) projection is stored by `TrkEffAnalysis` but not plotted.
- The header text is scaled down in narrow pads, so the sample label stays
  inside the pad; keep labels short.

## Selection variants (Step 4A)

The same validated calculation runs once per variant. Only the configuration
changes (`-s key=value` overrides of the baseline config), so there is no
second implementation:

```bash
LABEL="sample label" analysis/run_trkeff_variants.sh  SKIMLIST  OUTDIR  [BASECFG]  [EVALLIST  [STEP2_BASELINE_ROOT]]
analysis/bin/TrkEffAnalysis -c analysis/config/trkeff_baseline.cfg -s trk_nmvtx_min=3 -s use_reco_zvtx=1 ...   # single variant
```

| key | values | meaning |
|---|---|---|
| `trk_nmvtx_min` | 0, 2, 3 | MVTX clusters on the track (TPC ≥ 18, pT and \|η\| cuts as baseline) |
| `use_reco_zvtx` | 0 (`ztrue`), 1 (`zreco`) | event selection on \|z_true\| < 8 cm, or on a reconstructed vertex with \|z_reco\| < 8 cm |
| `truth_min_tpc_layers` | 0 (`incl`), 18 (`track`) | denominator A: all accepted G4 primaries; denominator B: only those with truth clusters in ≥ 18 TPC layers |

**Trackable efficiency (B) is conditional.** It counts only particles that
actually reached at least 18 TPC layers, so losses before the TPC (decays,
interactions, stopping) are excluded by construction. B is never combined with
A, and only A corrects to the primary-particle level. The vertex-finding
efficiency (`VtxEffAnalysis`) is independent of all three switches.

**Outputs:**
- `OUTDIR/TrkEff_mvtx<M>_<ztrue|zreco>_<incl|track>.root` (+ `_summary.txt`,
  `.log`), with identical binning in every variant. Each file stores its
  selection as `TNamed selection`, for example
  "reco vertex, |z_reco| < 8 cm; TPC >= 18, MVTX >= 3; trackable (>= 18 TPC layers)".
- `OUTDIR/plots/` (`CompareTrkEffVariants.C`):
  - page 1: integrated N/D = ε for all 12 variants × 6 species;
  - pages 2–4: pions, kaons, (anti)protons, each with ε vs pT for three
    comparisons: MVTX ≥ 0/2/3; inclusive vs trackable; true-z vs reco-z.
    Every curve shows passed/total and its empty bins.
  - also single figures, `TrkEffVariants_plots.root` and
    `TrkEffVariants_numbers.txt`.

**Validation:**
- `ValidateTrkEff.C` per variant: independent recomputation from the evaluator,
  including the trackable cut on `ntp_gtrack.gnltpc`; the split input lists
  merged with `hadd` equal the full run.
- `ValidateTrkEffVariants.C` across variants:
  - the baseline variant equals Step 2 exactly;
  - NumQ ≤ Num ≤ Den;
  - the denominator is independent of the MVTX cut, and Num(MVTX3) ≤ Num(MVTX2) ≤ Num(MVTX0);
  - trackable ≤ inclusive;
  - selected events depend only on the vertex mode;
  - the labels are correct;
  - the plotted graphs and the numbers table equal the variant files.

## Final selection (v2, Step 5B/6)

`TrkEffAnalysis` v2 adds three optional numerator keys. They are off by default, and with them off every output is identical to v1 (checked during development: all 15 validated baseline and variant configurations, every object bin by bin).

| key | default | meaning |
|---|---|---|
| `trk_nintt_min` | 0 | INTT clusters on the track |
| `trk_dca_xy_max` | 0 (off) | \|`linXYReco`\| ≤ value [cm] |
| `trk_dca_z_max` | 0 (off) | \|`linZReco`\| ≤ value [cm] |

- `linXYReco` / `linZReco` (EvalSkim) is the straight-line DCA to the event vertex: the reconstructed vertex with the most tracks.
- DCA cuts therefore require `use_reco_zvtx = 1`; otherwise the program stops.
- A non-finite DCA fails the cut and is counted (`hTrkEffDcaStats`, summary).

**FINAL configuration** (`config/trkeff_final_mvtx3.cfg`):

ε_trk = P(selected track | generated primary, final reconstructed vertex, |z_reco| < 8 cm)

| part | requirement |
|---|---|
| event | final vertex (`SvtxVertexMap`, Nmvtx ≥ 3 on the silicon seed), \|z_reco\| < 8 cm |
| track | TPC ≥ 20, MVTX ≥ 3, INTT ≥ 1, \|DCA_xy\| and \|DCA_z\| ≤ 1.5 cm |
| denominator | no truth-η cut |
| binning | fine common grid (Δy = 0.1, 17 pT bins) |

The vertex-finding efficiency and the PID efficiency are **not** included; they are separate quantities.

```bash
analysis/run_trkeff_final.sh SKIMLIST OUTDIR VARIANTSDIR "label" [EVALLIST]
```

The script writes:
- FINAL;
- the MVTX ≥ 2 reference (systematic, not a correction map);
- the one-cut steps;
- a same-denominator table input;
- the independent evaluator validation (`ValidateTrkEff.C`, extended with its own INTT and straight-line DCA code);
- `ValidateFinalTrkEff.C`: the steps, the true-z / reco-z / FINAL table, and the conditioning check;
- `MakeFinalTrkEffMaps.C`: maps in the species windows, Corral hook format, `maps/TrkEffFinalMaps.pdf`.

**Species windows** (true y, true pT):

| species | \|y\| | pT [GeV/c] |
|---|---|---|
| π± | < 1.0 | ≥ 0.1 |
| K± | < 0.7 | 0.1–0.6, plus a 0.6–1.0 QA bin |
| p, p̄ | < 0.6 | 0.4–1.0 |

**Hook format:** `heffpty_<sp>` is a TH2D (x = y, y = pT) whose content is ε, or −1 where no particle was generated. The errors are in `_errlo` / `_errhi`.

**Vertex provenance:** `macros/VertexProvenance.C` checks, per event, that the skim's event vertex equals the DST `SvtxVertexMap`.

**Ties:** when several vertices share the maximal number of tracks, any tied vertex is accepted and the tied events are counted. The evaluator row order (used by EvalSkim, first tied row) and the DST map key order can differ. This affects 0.04 % of vertexed events (1M subset).

**Check C of `ValidateFinalTrkEff.C`:** uses an absolute tolerance |diff| < 0.005. At 1M the true conditioning effect (−0.002) is a 5σ difference, so a pure significance test no longer applies.

## Truth-eta switch (v3, Step 6A)

| key | values | meaning |
|---|---|---|
| `use_truth_eta_cut` | `0` / `1`; absent = legacy | `1`: generated primaries restricted to \|η_true\| < `truth_eta_cut_max` (denominator and numerator); `0`: no truth-η cut |
| `truth_eta_cut_max` | 1.1 | threshold used when the switch is `1` |
| `truth_abseta_max` | legacy | used only when `use_truth_eta_cut` is absent (v1/v2 configurations: > 0 = cut on). Must be ≤ 0 when the switch is set; otherwise the program stops |

- **Where the cut acts:** in the common truth selection, before either histogram is filled, so the numerator is always a subset of the denominator. The reconstructed \|η\| < 1.1 cut is unchanged.
- **Final configuration:** `use_truth_eta_cut = 0`.
- **Truth-η mode:**

  ```
  TRUTH_ETA_CUT=1 analysis/run_trkeff_final.sh SKIMLIST OUTDIR_trutheta VARIANTSDIR "label" [EVALLIST]
  ```

  Use a separate output directory per mode.

**Meaning of the two definitions**
- **Switch `0` (inclusive):** every generated primary in the (y, pT) bin is counted. A particle that leaves the \|η\| < 1.1 acceptance is an inefficiency. This corrects a (y, pT)-binned yield to all particles in that bin.
- **Switch `1`:** the efficiency is conditional on the particle being generated inside \|η\| < 1.1. The geometric acceptance must then be corrected separately.

The two definitions differ only where a (y, pT) bin is partly outside \|η\| < 1.1: \|y\| > 0.8 for pions, and from \|y\| ≈ 0.3–0.5 at low pT for K and p.

**Tools**
- `macros/TrkEffConfigUtil.h`: the one rule for the effective truth-η threshold, shared by the macros.
- `macros/ValidateTruthEtaModes.C`: cross-check of the two modes, down to particle level.
- `macros/CompareTruthEtaModes.C`: comparison PDF.
