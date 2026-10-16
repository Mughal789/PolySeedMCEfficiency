// =============================================================================
// VtxEffAnalysis.cc
//
// Vertex-finding efficiency study for sPHENIX p+p 200 GeV.
//
// Input : SvtxEvaluator ntuples (outputeval_NNNNNN.root) of Michael Peters'
//         pileup-free PYTHIA8 minimum-bias sample with polyseed tracking.
//         One collision per event, 1000 events per file.
// Output: ONE ROOT file containing
//           - additive histograms (raw counts only, so outputs of several jobs
//             can later be merged with `hadd` without changing the physics);
//           - optionally a slim per-collision TTree "evt".
//         All ratios (efficiencies, P(n|N), fractions) are computed AFTER
//         merging, by macros/PlotVtxEff.C. Ratios cannot be added, counts can.
//
// Physics quantities, per simulated collision (= per evaluator "event"):
//
//   n        Reconstructed good-track multiplicity (baseline definition):
//            number of ntp_track rows with
//              pt > 0.1 GeV/c, |eta| < 1.1, ntpc >= 18
//            and, optionally, nmaps >= good_nmvtx_min (MVTX clusters; config
//            key, default 0 = no MVTX requirement = baseline definition).
//            ALL reconstructed tracks are considered, whether or not they are
//            attached to a vertex and whether or not a vertex exists. Truth
//            information is NOT used for n, so n is what real data would see.
//
//   N        Generator-level primary charged multiplicity:
//            number of ntp_gtrack rows with
//              gprimary == 1, charge != 0, gpt > 0.1 GeV/c, |geta| < 1.1.
//            "Primary" means a GEANT4 primary particle, i.e. a particle handed
//            to GEANT4 by PYTHIA. Michael's card sets
//            ParticleDecays:limitCylinder = on (1e-6 mm), so PYTHIA decays only
//            extremely short-lived states. Weak decays (K0S, Lambda, Sigma,
//            Xi, ...) happen in GEANT4, and their daughters are G4 secondaries.
//            Photon conversions and hadronic-interaction products are also G4
//            secondaries. All of these are therefore excluded from N
//            automatically. Charged Sigma/Xi/Omega themselves ARE counted
//            (ALICE-like c*tau > 1 cm convention).
//            NOTE: the G4-primary definition is not identical to a PYTHIA
//            "final state" definition; see the summary for the small
//            differences (e.g. pi0 Dalitz e+-, D-meson daughters).
//
//   hasVertex  true if the event has >= 1 ntp_vertex row with a finite
//            position. ntp_vertex is filled from SvtxVertexMap tracking
//            vertices; MBD-only/fake GlobalVertexMap entries never appear
//            there. Events without a vertex still produce one all-NaN row,
//            which must be ignored.
//            We deliberately do NOT use GlobalVertexMap, ntp_track.vertexID
//            (can be a stale 0) or ntp_gpoint.vx (an unmatched reco vertex
//            is missing there). See summary_evaluator.txt.
//
//   epsilon_vtx(n) = N_events(hasVertex, n) / N_events(n)
//   epsilon_vtx(N) = N_events(hasVertex, N) / N_events(N)
//   P(n|N)         = N_events(N, n) / N_events(N)
//
// Truth matching of reconstructed tracks (QA ONLY - never changes n):
//   gtrackID  id of the G4 particle contributing most clusters (NaN or 0 if
//             none); positive for G4 primaries, negative for secondaries.
//   gprimary  1 = matched particle is a G4 primary, 0 = G4 secondary.
//   A "duplicate" is a good track whose matched G4 particle already has
//   another good track (earlier in the same event) matched to it.
//
// Event key: the evaluator "event" counter restarts at 0 in every file, so
// the unique key is (job number from the file name, event).
//
// Build : scripts/build.sh
// Run   : bin/VtxEffAnalysis -c config/vtxeff.cfg [-l list] [-o out.root]
//         (see README.md)
// =============================================================================

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TTree.h>

#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// Configuration (see config/vtxeff.cfg for documentation of every key)
// -----------------------------------------------------------------------------
struct Config
{
  std::string input_dir;
  std::string input_list;
  std::string output_file;
  bool write_event_tree = true;

  double good_pt_min = 0.1;
  double good_abseta_max = 1.1;
  int good_ntpc_min = 18;
  // Optional MVTX requirement on good tracks (ntp_track.nmaps; NaN counts as 0).
  // 0 = no requirement (baseline).
  int good_nmvtx_min = 0;

  double true_pt_min = 0.1;
  double true_abseta_max = 1.1;

  int mult_max = 100;
  double ztrue_max = 40.0;
  // Additional efficiency histograms are filled for |z_true| < zcut_abs.
  // The vertex finder's acceptance depends strongly on z (silicon acceptance),
  // so epsilon(n) inside a central z window separates the multiplicity
  // dependence from the geometric z acceptance. Baseline (all z) is unchanged.
  double zcut_abs = 10.0;
  // Process only events 0 .. max_events_per_file-1 of every file (0 = all).
  // For quick tests only; events beyond the limit are ignored consistently in
  // every tree, so denominators stay correct for the processed events.
  long max_events_per_file = 0;
};

// Trim leading/trailing white space.
static std::string Trim(const std::string &s)
{
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

// Read "key = value" lines; '#' starts a comment. Unknown keys are fatal.
static Config ReadConfig(const std::string &path)
{
  Config c;
  std::ifstream in(path);
  if (!in)
  {
    std::cerr << "FATAL: cannot open config file " << path << std::endl;
    std::exit(1);
  }
  std::string line;
  int lineno = 0;
  while (std::getline(in, line))
  {
    ++lineno;
    const auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = Trim(line);
    if (line.empty()) continue;
    const auto eq = line.find('=');
    if (eq == std::string::npos)
    {
      std::cerr << "FATAL: " << path << ":" << lineno << " has no '='" << std::endl;
      std::exit(1);
    }
    const std::string key = Trim(line.substr(0, eq));
    const std::string val = Trim(line.substr(eq + 1));
    if (key == "input_dir") c.input_dir = val;
    else if (key == "input_list") c.input_list = val;
    else if (key == "output_file") c.output_file = val;
    else if (key == "write_event_tree") c.write_event_tree = (std::stoi(val) != 0);
    else if (key == "good_pt_min") c.good_pt_min = std::stod(val);
    else if (key == "good_abseta_max") c.good_abseta_max = std::stod(val);
    else if (key == "good_ntpc_min") c.good_ntpc_min = std::stoi(val);
    else if (key == "good_nmvtx_min") c.good_nmvtx_min = std::stoi(val);
    else if (key == "true_pt_min") c.true_pt_min = std::stod(val);
    else if (key == "true_abseta_max") c.true_abseta_max = std::stod(val);
    else if (key == "mult_max") c.mult_max = std::stoi(val);
    else if (key == "ztrue_max") c.ztrue_max = std::stod(val);
    else if (key == "zcut_abs") c.zcut_abs = std::stod(val);
    else if (key == "max_events_per_file") c.max_events_per_file = std::stol(val);
    else
    {
      std::cerr << "FATAL: unknown config key '" << key << "' (" << path << ":" << lineno << ")" << std::endl;
      std::exit(1);
    }
  }
  return c;
}

// -----------------------------------------------------------------------------
// PDG code -> electric charge (in units of e)
// -----------------------------------------------------------------------------
// The evaluator stores the PDG code as a Float_t (gflavor). A float holds
// integers exactly only up to 2^24 = 16,777,216. Ordinary hadron/lepton codes
// are far below that, but nuclear codes 10LZZZAAAI (~1e9) are NOT exact:
// e.g. 1000010020 (deuteron) can come back off by up to +-64. Therefore:
//   - round the float to the nearest integer;
//   - for nuclei (|pdg| >= 1e9), take Z from digits ZZZ after rounding to the
//     nearest multiple of 10^4, which is immune to the +-64 float error:
//       Z = ((|pdg| + 5000) / 10000) % 1000
//   - otherwise ask ROOT's TDatabasePDG (charge is stored in units of |e|/3).
// Unknown codes are reported (counter) and treated as neutral.
class ChargeLookup
{
 public:
  // Returns charge in units of e; sets known=false if the code is unknown.
  int Charge(long long pdg, bool &known)
  {
    known = true;
    auto it = cache_.find(pdg);
    if (it != cache_.end())
    {
      known = it->second.second;
      return it->second.first;
    }
    int q = 0;
    const long long apdg = std::llabs(pdg);
    if (apdg >= 1000000000LL)
    {
      // Nucleus: 10LZZZAAAI. Sign of pdg gives matter/antimatter.
      const long long Z = ((apdg + 5000LL) / 10000LL) % 1000LL;
      q = static_cast<int>(pdg > 0 ? Z : -Z);
    }
    else
    {
      TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(static_cast<int>(pdg));
      if (p)
      {
        q = static_cast<int>(std::lround(p->Charge() / 3.0));
      }
      else
      {
        known = false;  // e.g. PDG 0 (G4 geantino/optical photon) or exotic code
      }
    }
    cache_[pdg] = {q, known};
    return q;
  }

 private:
  std::unordered_map<long long, std::pair<int, bool>> cache_;
};

// -----------------------------------------------------------------------------
// Per-collision record, filled while reading one file
// -----------------------------------------------------------------------------
struct EventRecord
{
  // bookkeeping / verification
  int nInfoRows = 0;     // rows in ntp_info for this event (must be 1)
  int nGpointRows = 0;   // rows in ntp_gpoint for this event (must be 1)
  int ntrkInfo = -1;     // ntp_info.ntrk = size of SvtxTrackMap (cross-check)

  // true collision vertex (ntp_gpoint)
  float xTrue = NAN, yTrue = NAN, zTrue = NAN;

  // reconstructed tracking vertices (ntp_vertex rows with finite position)
  int nVtx = 0;
  float xReco = NAN, yReco = NAN, zReco = NAN;  // vertex with most tracks
  int vtxNtracks = -1;

  // reconstructed tracks
  int nAll = 0;          // all ntp_track rows of the event
  int nGood = 0;         // n (good-track multiplicity)
  int nGoodPrim = 0;     // good, matched to a G4 primary (incl. duplicates)
  int nGoodSec = 0;      // good, matched to a G4 secondary (incl. duplicates)
  int nGoodUnm = 0;      // good, no truth match
  int nGoodDup = 0;      // good, matched particle already has a good track
  int nGoodDupPrim = 0;  // ... of which the particle is a primary
  int nGoodMVTX = 0;     // good, >= 1 MVTX cluster
  int nGoodNoMVTX = 0;   // good, 0 MVTX clusters
  int nGoodTPCOnly = 0;  // good, 0 MVTX and 0 INTT clusters
  int nTrkWithVtxId = 0; // tracks with vertexID >= 0 (diagnostic only)
  std::set<long long> goodTruthIds;  // matched G4 ids of good tracks (duplicate finding)

  // Silicon content of GOOD tracks (diagnostics; n itself is unchanged).
  // nGoodMVTX above is n_MVTX1 (>= 1 MVTX cluster).
  int nMVTX2 = 0;        // good tracks with >= 2 MVTX clusters
  int nMVTX3 = 0;        // good tracks with >= 3 MVTX clusters

  // Per-good-track truth information, kept in ntp_track order so that the
  // uniqueness classification can be done after all tracks of the event are
  // read (step 6). QA only.
  struct GoodTrack
  {
    long long id = 0;      // matched G4 track id (valid only if matched)
    bool matched = false;  // gtrackID finite and != 0
    bool prim = false;     // matched particle is a G4 primary (gprimary == 1)
    bool primInAcc = false;// ... and charged with gpt > true_pt_min, |geta| < true_abseta_max
    float nfromtruth = 0;  // clusters contributed by the matched particle
    float nhits = 0;       // all clusters on the track
    float quality = 0;     // track-fit quality (chi2/ndf)
    int ntpc = 0, nmvtx = 0, nintt = 0;
  };
  std::vector<GoodTrack> goodTracks;

  // Results of the uniqueness classification (filled in step 6).
  int nUniquePrim = 0;     // first good track of a G4 primary
  int nUniqueSec = 0;      // first good track of a G4 secondary
  int nAddPrim = 0;        // additional good tracks of an already matched primary
  int nAddSec = 0;         // additional good tracks of an already matched secondary
  int nUniquePrimAcc = 0;  // unique, primary, inside the N acceptance
  int nPrimOutAcc = 0;     // good tracks matched to a primary OUTSIDE the N acceptance (unique)

  // generator-level primaries (ntp_gtrack, gprimary==1)
  int N = 0;             // charged, in kinematic window (true multiplicity N)
  int NchAll = 0;        // charged, any pT / eta
  int Npip = 0, Npim = 0, NKp = 0, NKm = 0, Np = 0, Npbar = 0;  // in window
  int NprimAll = 0;      // all G4 primaries, any charge (cross-check vs gntrkall)
  int gntrkallInfo = -1; // ntp_info.gntrkall

  bool HasVertex() const { return nVtx > 0; }
};

// -----------------------------------------------------------------------------
// Global bookkeeping for the validation log
// -----------------------------------------------------------------------------
struct FileSummary
{
  std::string name;
  int job = -1;
  long nEvents = 0, nVertexed = 0, nVertexless = 0;
  long nTracks = 0, nGood = 0, nGoodInVertexless = 0;
  long nTrkInVertexless = 0;
  long nGtrackRows = 0, nPrimaries = 0;
  long nUnknownPdg = 0;
  long nVtxRowsNaN = 0, nVtxRowsFinite = 0;
  long nTrackNaN = 0;            // tracks with NaN pt/eta/ntpc (cannot be good)
  long nTrkWithVtxId = 0;        // diagnostic: tracks with vertexID>=0
  int errors = 0;                // structural check failures
  double seconds = 0;
};

// -----------------------------------------------------------------------------
// Histograms. All are RAW COUNTS (additive under hadd).
// Naming: h<quantity>_<sample>; sample = all | vtx | novtx.
// -----------------------------------------------------------------------------
struct Histos
{
  // event statistics (labelled bins)
  TH1D *hStats = nullptr;

  // multiplicity, fine (unit) bins
  TH1D *hn_all = nullptr, *hn_vtx = nullptr;     // reco n
  TH1D *hN_all = nullptr, *hN_vtx = nullptr;     // true N
  // multiplicity, wide variable bins (for high multiplicity)
  TH1D *hnW_all = nullptr, *hnW_vtx = nullptr;
  TH1D *hNW_all = nullptr, *hNW_vtx = nullptr;
  // same, restricted to |z_true| < zcut_abs (suffix _zc)
  TH1D *hn_all_zc = nullptr, *hn_vtx_zc = nullptr, *hN_all_zc = nullptr, *hN_vtx_zc = nullptr;
  TH1D *hnW_all_zc = nullptr, *hnW_vtx_zc = nullptr, *hNW_all_zc = nullptr, *hNW_vtx_zc = nullptr;
  // alternative QA multiplicities (all events)
  TH1D *hnAll_all = nullptr;                     // all reco tracks per event
  TH1D *hnGoodPrim_all = nullptr, *hnGoodSec_all = nullptr;
  TH1D *hnGoodMVTX_all = nullptr, *hnGoodNoMVTX_all = nullptr;
  TH1D *hnAll_vtx = nullptr, *hnAll_novtx = nullptr;
  TH1D *hn_novtx = nullptr;

  // vertex
  TH1D *hZtrue_all = nullptr, *hZtrue_vtx = nullptr, *hZtrue_novtx = nullptr;
  TH1D *hNvtx = nullptr;
  TH1D *hDz = nullptr, *hDx = nullptr, *hDy = nullptr;  // reco - true (cm)
  TH1D *hVtxNtracks = nullptr;
  TH2D *hDz_vs_n = nullptr;

  // response (x = true N, y = reco n)
  TH2D *hResp_all = nullptr, *hResp_vtx = nullptr;
  TH1D *hnMinusN_all = nullptr, *hnMinusN_vtx = nullptr;
  TH1D *hnOverN_all = nullptr;  // only N > 0

  // species (generator-level, in kinematic window): per-event multiplicities
  TH1D *hNpip = nullptr, *hNpim = nullptr, *hNKp = nullptr, *hNKm = nullptr, *hNp = nullptr, *hNpbar = nullptr;

  // per-track QA, split by event vertex status
  // index [0] = tracks in vertexed events, [1] = tracks in vertexless events
  TH1D *hPt_all[2]{}, *hPt_good[2]{};
  TH1D *hEta_all[2]{}, *hEta_good[2]{};
  TH1D *hNtpc_all[2]{}, *hNtpc_good[2]{};
  TH1D *hNmvtx_all[2]{}, *hNmvtx_good[2]{};
  TH1D *hNintt_all[2]{}, *hNintt_good[2]{};
  // composition of GOOD tracks, per track: bins
  //   1 primary(unique)  2 secondary(unique)  3 duplicate  4 unmatched
  TH1D *hComp_good[2]{};
  // also with/without MVTX and TPC-only (bins: 1 with MVTX, 2 no MVTX but INTT, 3 TPC-only)
  TH1D *hSilicon_good[2]{};
  // composition vs reco n (x = n, y = category 1..4), weight = #tracks
  TH2D *hComp_vs_n = nullptr;
  TH2D *hSilicon_vs_n = nullptr;
  // truth-matched multiplicities vs n
  TH2D *hGoodPrim_vs_n = nullptr;  // x = n, y = n_good_primary

  // ---- silicon-content diagnostics (n unchanged) ----
  // event counts vs n_MVTXk = #good tracks with >= k MVTX clusters
  TH1D *hnM1_all = nullptr, *hnM1_vtx = nullptr;
  TH1D *hnM2_all = nullptr, *hnM2_vtx = nullptr;
  TH1D *hnM3_all = nullptr, *hnM3_vtx = nullptr;
  // 2D numerator/denominator: x = n, y = n_MVTX2
  TH2D *hEff2D_all = nullptr, *hEff2D_vtx = nullptr;
  // epsilon(n) in classes of n_MVTX1: m0 (= 0), m1 (>= 1), m2 (>= 2)
  TH1D *hn_all_m0 = nullptr, *hn_vtx_m0 = nullptr;
  TH1D *hn_all_m1 = nullptr, *hn_vtx_m1 = nullptr;
  TH1D *hn_all_m2 = nullptr, *hn_vtx_m2 = nullptr;

  // ---- truth-matching diagnostics (QA only) ----
  // A (origin) x B (uniqueness), per good track:
  //   x: 1 G4 primary, 2 G4 secondary, 3 unmatched
  //   y: 1 unique truth match, 2 additional reco match, 3 unmatched
  TH2D *hOriginUnique = nullptr;
  // exclusive track classes with priority (see README): per good track
  //   1 unmatched, 2 additional match, 3 secondary, 4 primary outside N acceptance,
  //   5 primary inside N acceptance
  TH1D *hTrkClass_all = nullptr, *hTrkClass_N0 = nullptr;
  // per-event multiplicities of the categories
  TH1D *hnUniquePrim_all = nullptr, *hnUniqueSec_all = nullptr, *hnAdd_all = nullptr, *hnUnm_all = nullptr;
  // duplicate diagnostics: [0] particle with exactly one good track,
  // [1] first good track of a multiply matched particle, [2] additional tracks
  TH1D *hPurity[3]{}, *hNtpcDup[3]{}, *hQualDup[3]{}, *hMvtxDup[3]{};
  TH1D *hMatchesPerParticle = nullptr;      // good tracks per matched G4 particle
  TH1D *hAddMinusFirst_nfromtruth = nullptr; // nfromtruth(additional) - nfromtruth(first)
  // diagnostic response matrices (x = true N)
  TH2D *hResp_uprim = nullptr;     // y = unique primary-matched good tracks
  TH2D *hResp_uprimAcc = nullptr;  // y = unique primary-matched, primary inside N acceptance

  // ---- vertex-residual diagnostics (vertexed events) ----
  TH2D *hDz_vs_z = nullptr, *hDx_vs_z = nullptr, *hDy_vs_z = nullptr;  // x = true z
  TH1D *hXtrue_vtx = nullptr, *hYtrue_vtx = nullptr, *hXreco = nullptr, *hYreco = nullptr;
};

// Fine bins centred on integers 0..max.
static TH1D *MultHist(const char *name, const char *title, int max)
{
  auto *h = new TH1D(name, title, max + 1, -0.5, max + 0.5);
  h->Sumw2();
  return h;
}

// Wide bins: unit width at low multiplicity (where the efficiency changes most), then
// widening. Edges are half-integers so each bin contains whole multiplicities.
static std::vector<double> WideEdges(int max)
{
  std::vector<double> e;
  const int lows[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70};
  for (int v : lows)
    if (v <= max) e.push_back(v - 0.5);
  e.push_back(max + 0.5);
  return e;
}

static Histos BookHistos(const Config &cfg)
{
  Histos H;
  const int M = cfg.mult_max;
  const auto we = WideEdges(M);

  H.hStats = new TH1D("hStats", "Event statistics;;events", 30, 0.5, 30.5);
  const char *labels[] = {"files", "events", "vertexed", "vertexless", "zero tracks", "zero good tracks",
                          "n=0", "n=1", "n=2", "n=3", "n=4", "n=5",
                          "N=0", "N=0 & n>0", "n>N", "n<N", "n=N",
                          "has secondary", "has duplicate", "has unmatched",
                          "nVtx>1", "tracks", "good tracks", "G4 primaries",
                          "dup G4 primary rows", "uniquePrimAcc>N", "n_MVTX1=0", "n_MVTX1>=1", "n_MVTX1>=2",
                          "has additional match"};
  for (int i = 0; i < 30; ++i) H.hStats->GetXaxis()->SetBinLabel(i + 1, labels[i]);

  H.hn_all = MultHist("hn_all", "All collisions;reco good-track multiplicity n;collisions", M);
  H.hn_vtx = MultHist("hn_vtx", "Vertexed collisions;reco good-track multiplicity n;collisions", M);
  H.hn_novtx = MultHist("hn_novtx", "Vertexless collisions;reco good-track multiplicity n;collisions", M);
  H.hN_all = MultHist("hN_all", "All collisions;true primary charged multiplicity N;collisions", M);
  H.hN_vtx = MultHist("hN_vtx", "Vertexed collisions;true primary charged multiplicity N;collisions", M);
  H.hnW_all = new TH1D("hnW_all", "All collisions (wide bins);n;collisions", we.size() - 1, we.data());
  H.hnW_vtx = new TH1D("hnW_vtx", "Vertexed collisions (wide bins);n;collisions", we.size() - 1, we.data());
  H.hNW_all = new TH1D("hNW_all", "All collisions (wide bins);N;collisions", we.size() - 1, we.data());
  H.hNW_vtx = new TH1D("hNW_vtx", "Vertexed collisions (wide bins);N;collisions", we.size() - 1, we.data());

  {
    const char *zt = Form("|z_{true}| < %g cm", cfg.zcut_abs);
    H.hn_all_zc = MultHist("hn_all_zc", Form("All collisions, %s;n;collisions", zt), M);
    H.hn_vtx_zc = MultHist("hn_vtx_zc", Form("Vertexed collisions, %s;n;collisions", zt), M);
    H.hN_all_zc = MultHist("hN_all_zc", Form("All collisions, %s;N;collisions", zt), M);
    H.hN_vtx_zc = MultHist("hN_vtx_zc", Form("Vertexed collisions, %s;N;collisions", zt), M);
    H.hnW_all_zc = new TH1D("hnW_all_zc", Form("All, %s (wide);n;collisions", zt), we.size() - 1, we.data());
    H.hnW_vtx_zc = new TH1D("hnW_vtx_zc", Form("Vertexed, %s (wide);n;collisions", zt), we.size() - 1, we.data());
    H.hNW_all_zc = new TH1D("hNW_all_zc", Form("All, %s (wide);N;collisions", zt), we.size() - 1, we.data());
    H.hNW_vtx_zc = new TH1D("hNW_vtx_zc", Form("Vertexed, %s (wide);N;collisions", zt), we.size() - 1, we.data());
  }

  H.hnAll_all = MultHist("hnAll_all", "All collisions;all reco tracks per collision;collisions", 3 * M);
  H.hnAll_vtx = MultHist("hnAll_vtx", "Vertexed;all reco tracks per collision;collisions", 3 * M);
  H.hnAll_novtx = MultHist("hnAll_novtx", "Vertexless;all reco tracks per collision;collisions", 3 * M);
  H.hnGoodPrim_all = MultHist("hnGoodPrim_all", "Good tracks matched to G4 primaries;n_{good,prim};collisions", M);
  H.hnGoodSec_all = MultHist("hnGoodSec_all", "Good tracks matched to G4 secondaries;n_{good,sec};collisions", M);
  H.hnGoodMVTX_all = MultHist("hnGoodMVTX_all", "Good tracks with MVTX clusters;n_{good,MVTX};collisions", M);
  H.hnGoodNoMVTX_all = MultHist("hnGoodNoMVTX_all", "Good tracks without MVTX clusters;n_{good,noMVTX};collisions", M);

  const double zm = cfg.ztrue_max;
  H.hZtrue_all = new TH1D("hZtrue_all", "All collisions;true z_{vtx} (cm);collisions", 40, -zm, zm);
  H.hZtrue_vtx = new TH1D("hZtrue_vtx", "Vertexed;true z_{vtx} (cm);collisions", 40, -zm, zm);
  H.hZtrue_novtx = new TH1D("hZtrue_novtx", "Vertexless;true z_{vtx} (cm);collisions", 40, -zm, zm);
  H.hNvtx = new TH1D("hNvtx", "Reconstructed tracking vertices per collision;N_{vtx};collisions", 6, -0.5, 5.5);
  H.hDz = new TH1D("hDz", "Vertexed collisions;z_{reco} - z_{true} (cm);collisions", 200, -0.5, 0.5);
  H.hDx = new TH1D("hDx", "Vertexed collisions;x_{reco} - x_{true} (cm);collisions", 200, -0.3, 0.3);
  H.hDy = new TH1D("hDy", "Vertexed collisions;y_{reco} - y_{true} (cm);collisions", 200, -0.3, 0.3);
  H.hVtxNtracks = new TH1D("hVtxNtracks", "Vertexed collisions;tracks attached to vertex;collisions", 61, -0.5, 60.5);
  H.hDz_vs_n = new TH2D("hDz_vs_n", "Vertexed collisions;n;z_{reco} - z_{true} (cm)", M + 1, -0.5, M + 0.5, 100, -0.5, 0.5);

  H.hResp_all = new TH2D("hResp_all", "All collisions;true N;reco n", M + 1, -0.5, M + 0.5, M + 1, -0.5, M + 0.5);
  H.hResp_vtx = new TH2D("hResp_vtx", "Vertexed collisions;true N;reco n", M + 1, -0.5, M + 0.5, M + 1, -0.5, M + 0.5);
  H.hnMinusN_all = new TH1D("hnMinusN_all", "All collisions;n - N;collisions", 2 * M + 1, -M - 0.5, M + 0.5);
  H.hnMinusN_vtx = new TH1D("hnMinusN_vtx", "Vertexed collisions;n - N;collisions", 2 * M + 1, -M - 0.5, M + 0.5);
  H.hnOverN_all = new TH1D("hnOverN_all", "Collisions with N > 0;n / N;collisions", 100, 0., 10.);

  H.hNpip = MultHist("hNpip", "#pi^{+} primaries;multiplicity;collisions", 50);
  H.hNpim = MultHist("hNpim", "#pi^{-} primaries;multiplicity;collisions", 50);
  H.hNKp = MultHist("hNKp", "K^{+} primaries;multiplicity;collisions", 50);
  H.hNKm = MultHist("hNKm", "K^{-} primaries;multiplicity;collisions", 50);
  H.hNp = MultHist("hNp", "p primaries;multiplicity;collisions", 50);
  H.hNpbar = MultHist("hNpbar", "#bar{p} primaries;multiplicity;collisions", 50);

  const char *tag[2] = {"vtxEvt", "novtxEvt"};
  const char *ttl[2] = {"tracks in vertexed collisions", "tracks in vertexless collisions"};
  for (int k = 0; k < 2; ++k)
  {
    H.hPt_all[k] = new TH1D(Form("hPt_all_%s", tag[k]), Form("All %s;p_{T} (GeV/c);tracks", ttl[k]), 100, 0., 5.);
    H.hPt_good[k] = new TH1D(Form("hPt_good_%s", tag[k]), Form("Good %s;p_{T} (GeV/c);tracks", ttl[k]), 100, 0., 5.);
    H.hEta_all[k] = new TH1D(Form("hEta_all_%s", tag[k]), Form("All %s;#eta;tracks", ttl[k]), 60, -1.5, 1.5);
    H.hEta_good[k] = new TH1D(Form("hEta_good_%s", tag[k]), Form("Good %s;#eta;tracks", ttl[k]), 60, -1.5, 1.5);
    H.hNtpc_all[k] = new TH1D(Form("hNtpc_all_%s", tag[k]), Form("All %s;TPC clusters per track;tracks", ttl[k]), 61, -0.5, 60.5);
    H.hNtpc_good[k] = new TH1D(Form("hNtpc_good_%s", tag[k]), Form("Good %s;TPC clusters per track;tracks", ttl[k]), 61, -0.5, 60.5);
    H.hNmvtx_all[k] = new TH1D(Form("hNmvtx_all_%s", tag[k]), Form("All %s;MVTX clusters per track;tracks", ttl[k]), 8, -0.5, 7.5);
    H.hNmvtx_good[k] = new TH1D(Form("hNmvtx_good_%s", tag[k]), Form("Good %s;MVTX clusters per track;tracks", ttl[k]), 8, -0.5, 7.5);
    H.hNintt_all[k] = new TH1D(Form("hNintt_all_%s", tag[k]), Form("All %s;INTT clusters per track;tracks", ttl[k]), 8, -0.5, 7.5);
    H.hNintt_good[k] = new TH1D(Form("hNintt_good_%s", tag[k]), Form("Good %s;INTT clusters per track;tracks", ttl[k]), 8, -0.5, 7.5);
    H.hComp_good[k] = new TH1D(Form("hComp_good_%s", tag[k]), Form("Good %s;truth category;tracks", ttl[k]), 4, 0.5, 4.5);
    H.hSilicon_good[k] = new TH1D(Form("hSilicon_good_%s", tag[k]), Form("Good %s;silicon content;tracks", ttl[k]), 3, 0.5, 3.5);
  }
  H.hComp_vs_n = new TH2D("hComp_vs_n", "Good tracks;reco n;truth category", M + 1, -0.5, M + 0.5, 4, 0.5, 4.5);
  H.hSilicon_vs_n = new TH2D("hSilicon_vs_n", "Good tracks;reco n;silicon content", M + 1, -0.5, M + 0.5, 3, 0.5, 3.5);
  H.hGoodPrim_vs_n = new TH2D("hGoodPrim_vs_n", "All collisions;reco n;good tracks matched to G4 primaries",
                              M + 1, -0.5, M + 0.5, M + 1, -0.5, M + 0.5);

  // silicon diagnostics
  H.hnM1_all = MultHist("hnM1_all", "All;n_{MVTX1};collisions", 30);
  H.hnM1_vtx = MultHist("hnM1_vtx", "Vertexed;n_{MVTX1};collisions", 30);
  H.hnM2_all = MultHist("hnM2_all", "All;n_{MVTX2};collisions", 30);
  H.hnM2_vtx = MultHist("hnM2_vtx", "Vertexed;n_{MVTX2};collisions", 30);
  H.hnM3_all = MultHist("hnM3_all", "All;n_{MVTX3};collisions", 30);
  H.hnM3_vtx = MultHist("hnM3_vtx", "Vertexed;n_{MVTX3};collisions", 30);
  H.hEff2D_all = new TH2D("hEff2D_all", "All collisions;n;n_{MVTX2}", 41, -0.5, 40.5, 11, -0.5, 10.5);
  H.hEff2D_vtx = new TH2D("hEff2D_vtx", "Vertexed collisions;n;n_{MVTX2}", 41, -0.5, 40.5, 11, -0.5, 10.5);
  H.hn_all_m0 = MultHist("hn_all_m0", "All, n_{MVTX1} = 0;n;collisions", M);
  H.hn_vtx_m0 = MultHist("hn_vtx_m0", "Vertexed, n_{MVTX1} = 0;n;collisions", M);
  H.hn_all_m1 = MultHist("hn_all_m1", "All, n_{MVTX1} #geq 1;n;collisions", M);
  H.hn_vtx_m1 = MultHist("hn_vtx_m1", "Vertexed, n_{MVTX1} #geq 1;n;collisions", M);
  H.hn_all_m2 = MultHist("hn_all_m2", "All, n_{MVTX1} #geq 2;n;collisions", M);
  H.hn_vtx_m2 = MultHist("hn_vtx_m2", "Vertexed, n_{MVTX1} #geq 2;n;collisions", M);

  // truth diagnostics
  H.hOriginUnique = new TH2D("hOriginUnique", "Good tracks;truth origin;reconstruction uniqueness", 3, 0.5, 3.5, 3, 0.5, 3.5);
  {
    const char *ox[3] = {"G4 primary", "G4 secondary", "unmatched"};
    const char *uy[3] = {"unique match", "additional match", "unmatched"};
    for (int b = 0; b < 3; ++b)
    {
      H.hOriginUnique->GetXaxis()->SetBinLabel(b + 1, ox[b]);
      H.hOriginUnique->GetYaxis()->SetBinLabel(b + 1, uy[b]);
    }
  }
  H.hTrkClass_all = new TH1D("hTrkClass_all", "Good tracks, all collisions;track class;tracks", 5, 0.5, 5.5);
  H.hTrkClass_N0 = new TH1D("hTrkClass_N0", "Good tracks in collisions with N = 0;track class;tracks", 5, 0.5, 5.5);
  {
    const char *cl[5] = {"unmatched", "additional match", "secondary", "primary, outside N acc.", "primary, inside N acc."};
    for (int b = 0; b < 5; ++b)
    {
      H.hTrkClass_all->GetXaxis()->SetBinLabel(b + 1, cl[b]);
      H.hTrkClass_N0->GetXaxis()->SetBinLabel(b + 1, cl[b]);
    }
  }
  H.hnUniquePrim_all = MultHist("hnUniquePrim_all", "All;unique primary-matched good tracks;collisions", M);
  H.hnUniqueSec_all = MultHist("hnUniqueSec_all", "All;unique secondary-matched good tracks;collisions", M);
  H.hnAdd_all = MultHist("hnAdd_all", "All;additional reco matches;collisions", M);
  H.hnUnm_all = MultHist("hnUnm_all", "All;unmatched good tracks;collisions", M);
  {
    const char *dn[3] = {"single", "multiFirst", "multiAdd"};
    const char *dt[3] = {"particle with one good track", "first track of multiply matched particle",
                         "additional track of multiply matched particle"};
    for (int i = 0; i < 3; ++i)
    {
      H.hPurity[i] = new TH1D(Form("hPurity_%s", dn[i]), Form("%s;nfromtruth / nhits;tracks", dt[i]), 50, 0., 1.0001);
      H.hNtpcDup[i] = new TH1D(Form("hNtpcDup_%s", dn[i]), Form("%s;TPC clusters;tracks", dt[i]), 61, -0.5, 60.5);
      H.hQualDup[i] = new TH1D(Form("hQualDup_%s", dn[i]), Form("%s;track quality (#chi^{2}/ndf);tracks", dt[i]), 100, 0., 20.);
      H.hMvtxDup[i] = new TH1D(Form("hMvtxDup_%s", dn[i]), Form("%s;MVTX clusters;tracks", dt[i]), 8, -0.5, 7.5);
    }
  }
  H.hMatchesPerParticle = new TH1D("hMatchesPerParticle", "Matched G4 particles;good tracks per G4 particle;particles", 10, 0.5, 10.5);
  H.hAddMinusFirst_nfromtruth = new TH1D("hAddMinusFirst_nfromtruth", "Additional tracks;nfromtruth(additional) - nfromtruth(first);tracks", 121, -60.5, 60.5);
  H.hResp_uprim = new TH2D("hResp_uprim", "All collisions;true N;unique primary-matched good tracks", M + 1, -0.5, M + 0.5, M + 1, -0.5, M + 0.5);
  H.hResp_uprimAcc = new TH2D("hResp_uprimAcc", "All collisions;true N;unique primary-matched good tracks (primary in N acceptance)", M + 1, -0.5, M + 0.5, M + 1, -0.5, M + 0.5);

  // residual diagnostics
  H.hDz_vs_z = new TH2D("hDz_vs_z", "Vertexed;true z (cm);z_{reco} - z_{true} (cm)", 20, -20., 20., 200, -0.5, 0.5);
  H.hDx_vs_z = new TH2D("hDx_vs_z", "Vertexed;true z (cm);x_{reco} - x_{true} (cm)", 20, -20., 20., 200, -0.3, 0.3);
  H.hDy_vs_z = new TH2D("hDy_vs_z", "Vertexed;true z (cm);y_{reco} - y_{true} (cm)", 20, -20., 20., 200, -0.3, 0.3);
  H.hXtrue_vtx = new TH1D("hXtrue_vtx", "Vertexed;x (cm);collisions", 120, -0.3, 0.3);
  H.hYtrue_vtx = new TH1D("hYtrue_vtx", "Vertexed;y (cm);collisions", 120, -0.1, 0.5);
  H.hXreco = new TH1D("hXreco", "Vertexed;x (cm);collisions", 120, -0.3, 0.3);
  H.hYreco = new TH1D("hYreco", "Vertexed;y (cm);collisions", 120, -0.1, 0.5);

  const char *comp[4] = {"primary", "secondary", "duplicate", "unmatched"};
  const char *sil[3] = {"MVTX", "INTT only", "TPC only"};
  for (int k = 0; k < 2; ++k)
  {
    for (int b = 0; b < 4; ++b) H.hComp_good[k]->GetXaxis()->SetBinLabel(b + 1, comp[b]);
    for (int b = 0; b < 3; ++b) H.hSilicon_good[k]->GetXaxis()->SetBinLabel(b + 1, sil[b]);
  }
  for (int b = 0; b < 4; ++b) H.hComp_vs_n->GetYaxis()->SetBinLabel(b + 1, comp[b]);
  for (int b = 0; b < 3; ++b) H.hSilicon_vs_n->GetYaxis()->SetBinLabel(b + 1, sil[b]);
  return H;
}

// Clamp a multiplicity into the histogram range so that overflow events are
// not lost from totals (they land in the last bin; their number is reported).
static int Clamp(int v, int max, long &nOverflow)
{
  if (v > max)
  {
    ++nOverflow;
    return max;
  }
  return v;
}

// -----------------------------------------------------------------------------
// Helpers to read TNtuple branches safely
// -----------------------------------------------------------------------------
// Enable only the listed branches and bind them to float variables. Exit with
// an error if a branch is missing (better than silently reading zeros).
static void Bind(TTree *t, const std::vector<std::pair<const char *, float *>> &br)
{
  t->SetBranchStatus("*", 0);
  for (const auto &b : br)
  {
    if (!t->GetBranch(b.first))
    {
      std::cerr << "FATAL: branch " << b.first << " missing in tree " << t->GetName() << std::endl;
      std::exit(2);
    }
    t->SetBranchStatus(b.first, 1);
    t->SetBranchAddress(b.first, b.second);
  }
}

static bool Finite(float x) { return std::isfinite(x); }

// Extract the job number from ".../outputeval_NNNNNN.root". -1 if not found.
static int JobFromName(const std::string &path)
{
  static const std::regex re("outputeval_([0-9]+)\\.root$");
  std::smatch m;
  if (std::regex_search(path, m, re)) return std::stoi(m[1]);
  return -1;
}

// -----------------------------------------------------------------------------
// Process one evaluator file
// -----------------------------------------------------------------------------
static FileSummary ProcessFile(const std::string &path, const Config &cfg, Histos &H, TTree *evtTree,
                               ChargeLookup &chg, std::map<std::string, long> &warn)
{
  FileSummary S;
  S.name = path;
  S.job = JobFromName(path);
  TStopwatch sw;
  sw.Start();

  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    S.errors++;
    return S;
  }
  // TFile::Get returns the highest key cycle, i.e. the complete tree (some
  // ntuples have an older autosave cycle as well).
  auto *tInfo = dynamic_cast<TTree *>(f->Get("ntp_info"));
  auto *tGpoint = dynamic_cast<TTree *>(f->Get("ntp_gpoint"));
  auto *tVertex = dynamic_cast<TTree *>(f->Get("ntp_vertex"));
  auto *tTrack = dynamic_cast<TTree *>(f->Get("ntp_track"));
  auto *tGtrack = dynamic_cast<TTree *>(f->Get("ntp_gtrack"));
  if (!tInfo || !tGpoint || !tVertex || !tTrack || !tGtrack)
  {
    std::cerr << "ERROR: missing evaluator tree in " << path << std::endl;
    S.errors++;
    return S;
  }

  // ---- 1. ntp_info: defines the list of simulated collisions ---------------
  // One row per event. The number of rows sets the denominator, so every
  // collision counts, including those without tracks or vertices.
  const long nEvtFile = tInfo->GetEntries();
  // Number of collisions actually used (all, unless a test limit is set).
  const long nEvt = (cfg.max_events_per_file > 0 && cfg.max_events_per_file < nEvtFile) ? cfg.max_events_per_file : nEvtFile;
  std::vector<EventRecord> ev(nEvt);
  {
    float event, ntrk, gntrkall;
    Bind(tInfo, {{"event", &event}, {"ntrk", &ntrk}, {"gntrkall", &gntrkall}});
    for (long i = 0; i < nEvtFile; ++i)
    {
      tInfo->GetEntry(i);
      const long e = std::lround(event);
      if (e < 0 || e >= nEvtFile)
      {
        std::cerr << "ERROR: ntp_info event " << event << " outside [0," << nEvtFile << ")" << std::endl;
        S.errors++;
        continue;
      }
      if (e >= nEvt) continue;  // beyond the test limit
      ev[e].nInfoRows++;
      ev[e].ntrkInfo = std::lround(ntrk);
      ev[e].gntrkallInfo = std::lround(gntrkall);
    }
  }

  // Helper: map an event value from another tree onto the ntp_info list.
  auto lookup = [&](float event, const char *tree) -> EventRecord * {
    const long e = std::lround(event);
    if (e >= nEvt && e < nEvtFile) return nullptr;  // beyond the test limit: ignore silently
    if (e < 0 || e >= nEvt || ev[e].nInfoRows == 0)
    {
      warn[std::string("event not in ntp_info (") + tree + ")"]++;
      S.errors++;
      return nullptr;
    }
    return &ev[e];
  };

  // ---- 2. ntp_gpoint: true collision vertex (one row per collision) --------
  {
    float event, gvx, gvy, gvz;
    Bind(tGpoint, {{"event", &event}, {"gvx", &gvx}, {"gvy", &gvy}, {"gvz", &gvz}});
    for (long i = 0; i < tGpoint->GetEntries(); ++i)
    {
      tGpoint->GetEntry(i);
      EventRecord *r = lookup(event, "ntp_gpoint");
      if (!r) continue;
      r->nGpointRows++;
      r->xTrue = gvx;
      r->yTrue = gvy;
      r->zTrue = gvz;
    }
  }

  // ---- 3. ntp_vertex: reconstructed tracking vertices -----------------------
  // A row with finite (vx,vy,vz) is a real SvtxVertexMap vertex. Events with no
  // vertex get one all-NaN placeholder row, which is skipped here.
  {
    float event, vx, vy, vz, ntracks;
    Bind(tVertex, {{"event", &event}, {"vx", &vx}, {"vy", &vy}, {"vz", &vz}, {"ntracks", &ntracks}});
    for (long i = 0; i < tVertex->GetEntries(); ++i)
    {
      tVertex->GetEntry(i);
      EventRecord *r = lookup(event, "ntp_vertex");
      if (!r) continue;
      if (!(Finite(vx) && Finite(vy) && Finite(vz)))
      {
        S.nVtxRowsNaN++;
        continue;
      }
      S.nVtxRowsFinite++;
      r->nVtx++;
      // If ever more than one vertex: keep the one with most tracks for the
      // residual plots (never happened in the audited files).
      const int nt = Finite(ntracks) ? static_cast<int>(std::lround(ntracks)) : 0;
      if (r->nVtx == 1 || nt > r->vtxNtracks)
      {
        r->xReco = vx;
        r->yReco = vy;
        r->zReco = vz;
        r->vtxNtracks = nt;
      }
    }
  }

  // ---- 4. ntp_track: ALL reconstructed tracks --------------------------------
  // No requirement on vertexID, silicon hits, DCA or truth: n is built from
  // every reconstructed track passing the good-track kinematic + TPC cuts.
  {
    float event, pt, eta, ntpc, nmaps, nintt, gtrackID, gprimary, vertexID;
    float gpt, geta, gflavor, nfromtruth, nhits, quality;  // truth-QA only
    Bind(tTrack, {{"event", &event}, {"pt", &pt}, {"eta", &eta}, {"ntpc", &ntpc}, {"nmaps", &nmaps},
                  {"nintt", &nintt}, {"gtrackID", &gtrackID}, {"gprimary", &gprimary}, {"vertexID", &vertexID},
                  {"gpt", &gpt}, {"geta", &geta}, {"gflavor", &gflavor}, {"nfromtruth", &nfromtruth},
                  {"nhits", &nhits}, {"quality", &quality}});
    for (long i = 0; i < tTrack->GetEntries(); ++i)
    {
      tTrack->GetEntry(i);
      EventRecord *r = lookup(event, "ntp_track");
      if (!r) continue;
      r->nAll++;
      S.nTracks++;
      if (Finite(vertexID) && vertexID >= 0) r->nTrkWithVtxId++;

      if (!(Finite(pt) && Finite(eta) && Finite(ntpc)))
      {
        S.nTrackNaN++;  // cannot satisfy the cuts; counted but never "good"
        continue;
      }
      const int itpc = static_cast<int>(std::lround(ntpc));
      const int imvtx = Finite(nmaps) ? static_cast<int>(std::lround(nmaps)) : 0;
      const int iintt = Finite(nintt) ? static_cast<int>(std::lround(nintt)) : 0;

      // Per-track QA needs the event's vertex status, known from step 3.
      const int k = r->HasVertex() ? 0 : 1;
      H.hPt_all[k]->Fill(pt);
      H.hEta_all[k]->Fill(eta);
      H.hNtpc_all[k]->Fill(itpc);
      H.hNmvtx_all[k]->Fill(imvtx);
      H.hNintt_all[k]->Fill(iintt);

      // ---- good-track definition (baseline n) ----
      const bool good = (pt > cfg.good_pt_min) && (std::fabs(eta) < cfg.good_abseta_max) && (itpc >= cfg.good_ntpc_min) &&
                        (imvtx >= cfg.good_nmvtx_min);
      if (!good) continue;
      r->nGood++;
      H.hPt_good[k]->Fill(pt);
      H.hEta_good[k]->Fill(eta);
      H.hNtpc_good[k]->Fill(itpc);
      H.hNmvtx_good[k]->Fill(imvtx);
      H.hNintt_good[k]->Fill(iintt);

      // silicon content (QA)
      if (imvtx > 0) r->nGoodMVTX++;
      else r->nGoodNoMVTX++;
      if (imvtx >= 2) r->nMVTX2++;
      if (imvtx >= 3) r->nMVTX3++;
      if (imvtx == 0 && iintt == 0) r->nGoodTPCOnly++;

      // keep the per-track truth information for the classification in step 6
      {
        EventRecord::GoodTrack g;
        g.matched = Finite(gtrackID) && std::lround(gtrackID) != 0;
        if (g.matched)
        {
          g.id = std::llround(gtrackID);
          g.prim = Finite(gprimary) && std::lround(gprimary) == 1;
          bool known = true;
          const int qTrue = Finite(gflavor) ? chg.Charge(std::llround(gflavor), known) : 0;
          g.primInAcc = g.prim && qTrue != 0 && Finite(gpt) && Finite(geta) && gpt > cfg.true_pt_min &&
                        std::fabs(geta) < cfg.true_abseta_max;
        }
        g.nfromtruth = Finite(nfromtruth) ? nfromtruth : 0.f;
        g.nhits = Finite(nhits) ? nhits : 0.f;
        g.quality = Finite(quality) ? quality : -1.f;
        g.ntpc = itpc;
        g.nmvtx = imvtx;
        g.nintt = iintt;
        r->goodTracks.push_back(g);
      }

      // truth composition (QA only)
      const bool matched = Finite(gtrackID) && std::lround(gtrackID) != 0;
      if (!matched)
      {
        r->nGoodUnm++;
      }
      else
      {
        const bool isPrim = Finite(gprimary) && std::lround(gprimary) == 1;
        if (isPrim) r->nGoodPrim++;
        else r->nGoodSec++;
        const long long id = std::llround(gtrackID);
        if (!r->goodTruthIds.insert(id).second)
        {
          r->nGoodDup++;  // this G4 particle already has a good track
          if (isPrim) r->nGoodDupPrim++;
        }
      }
    }
  }

  // ---- 5. ntp_gtrack: generator-level (G4) primaries -------------------------
  // ~950k rows per file (all G4 particles). Only 5 small branches are read,
  // and only primaries (gprimary==1, ~30-40 per event) are processed further.
  {
    float event, gprimary, gflavor, gpt, geta, gtrackIDg;
    Bind(tGtrack, {{"event", &event}, {"gprimary", &gprimary}, {"gflavor", &gflavor}, {"gpt", &gpt}, {"geta", &geta},
                   {"gtrackID", &gtrackIDg}});
    std::map<long, std::set<long long>> primIdsPerEvent;  // detects duplicated primary rows
    const long nG = tGtrack->GetEntries();
    S.nGtrackRows = nG;
    for (long i = 0; i < nG; ++i)
    {
      tGtrack->GetEntry(i);
      if (!(Finite(gprimary) && std::lround(gprimary) == 1)) continue;
      EventRecord *r = lookup(event, "ntp_gtrack");
      if (!r) continue;
      if (Finite(gtrackIDg) && !primIdsPerEvent[std::lround(event)].insert(std::llround(gtrackIDg)).second)
      {
        warn["duplicated G4-primary row in ntp_gtrack (same event, same gtrackID) - skipped"]++;
        H.hStats->Fill("dup G4 primary rows", 1);
        continue;
      }
      r->NprimAll++;
      S.nPrimaries++;
      const long long pdg = std::llround(gflavor);
      bool known = true;
      const int q = chg.Charge(pdg, known);
      if (!known)
      {
        S.nUnknownPdg++;
        warn["unknown PDG code " + std::to_string(pdg) + " (treated as neutral)"]++;
      }
      if (q == 0) continue;
      r->NchAll++;
      if (!(Finite(gpt) && Finite(geta))) continue;
      if (gpt > cfg.true_pt_min && std::fabs(geta) < cfg.true_abseta_max)
      {
        r->N++;
        switch (pdg)
        {
        case 211: r->Npip++; break;
        case -211: r->Npim++; break;
        case 321: r->NKp++; break;
        case -321: r->NKm++; break;
        case 2212: r->Np++; break;
        case -2212: r->Npbar++; break;
        default: break;
        }
      }
    }
  }

  // ---- 6. Per-event checks and histogram filling ----------------------------
  S.nEvents = nEvt;
  H.hStats->Fill("files", 1);
  // per-event tree variables
  static int t_job, t_event, t_n, t_nAll, t_nGoodPrim, t_nGoodSec, t_nGoodUnm, t_nGoodDup, t_nGoodMVTX,
      t_nGoodNoMVTX, t_nGoodTPCOnly, t_N, t_NchAll, t_Npip, t_Npim, t_NKp, t_NKm, t_Np, t_Npbar, t_hasVtx, t_nVtx,
      t_vtxNtracks, t_nMVTX2, t_nMVTX3, t_nINTTonly, t_nUniquePrim, t_nUniqueSec, t_nAddPrim, t_nAddSec,
      t_nUniquePrimAcc, t_nPrimOutAcc;
  static float t_zTrue, t_xTrue, t_yTrue, t_zReco, t_xReco, t_yReco;
  if (evtTree && evtTree->GetNbranches() == 0)
  {
    evtTree->Branch("job", &t_job, "job/I");
    evtTree->Branch("event", &t_event, "event/I");
    evtTree->Branch("n", &t_n, "n/I");
    evtTree->Branch("nAll", &t_nAll, "nAll/I");
    evtTree->Branch("nGoodPrim", &t_nGoodPrim, "nGoodPrim/I");
    evtTree->Branch("nGoodSec", &t_nGoodSec, "nGoodSec/I");
    evtTree->Branch("nGoodUnm", &t_nGoodUnm, "nGoodUnm/I");
    evtTree->Branch("nGoodDup", &t_nGoodDup, "nGoodDup/I");
    evtTree->Branch("nGoodMVTX", &t_nGoodMVTX, "nGoodMVTX/I");
    evtTree->Branch("nGoodNoMVTX", &t_nGoodNoMVTX, "nGoodNoMVTX/I");
    evtTree->Branch("nGoodTPCOnly", &t_nGoodTPCOnly, "nGoodTPCOnly/I");
    evtTree->Branch("N", &t_N, "N/I");
    evtTree->Branch("NchAll", &t_NchAll, "NchAll/I");
    evtTree->Branch("Npip", &t_Npip, "Npip/I");
    evtTree->Branch("Npim", &t_Npim, "Npim/I");
    evtTree->Branch("NKp", &t_NKp, "NKp/I");
    evtTree->Branch("NKm", &t_NKm, "NKm/I");
    evtTree->Branch("Np", &t_Np, "Np/I");
    evtTree->Branch("Npbar", &t_Npbar, "Npbar/I");
    evtTree->Branch("hasVtx", &t_hasVtx, "hasVtx/I");
    evtTree->Branch("nVtx", &t_nVtx, "nVtx/I");
    evtTree->Branch("vtxNtracks", &t_vtxNtracks, "vtxNtracks/I");
    evtTree->Branch("zTrue", &t_zTrue, "zTrue/F");
    evtTree->Branch("xTrue", &t_xTrue, "xTrue/F");
    evtTree->Branch("yTrue", &t_yTrue, "yTrue/F");
    evtTree->Branch("zReco", &t_zReco, "zReco/F");
    evtTree->Branch("xReco", &t_xReco, "xReco/F");
    evtTree->Branch("yReco", &t_yReco, "yReco/F");
    // n_MVTX1 = nGoodMVTX (branch above)
    evtTree->Branch("nMVTX2", &t_nMVTX2, "nMVTX2/I");
    evtTree->Branch("nMVTX3", &t_nMVTX3, "nMVTX3/I");
    evtTree->Branch("nINTTonly", &t_nINTTonly, "nINTTonly/I");
    evtTree->Branch("nUniquePrim", &t_nUniquePrim, "nUniquePrim/I");
    evtTree->Branch("nUniqueSec", &t_nUniqueSec, "nUniqueSec/I");
    evtTree->Branch("nAddPrim", &t_nAddPrim, "nAddPrim/I");
    evtTree->Branch("nAddSec", &t_nAddSec, "nAddSec/I");
    evtTree->Branch("nUniquePrimAcc", &t_nUniquePrimAcc, "nUniquePrimAcc/I");
    evtTree->Branch("nPrimOutAcc", &t_nPrimOutAcc, "nPrimOutAcc/I");
  }

  static long nOverflow = 0;
  const int M = cfg.mult_max;
  for (long e = 0; e < nEvt; ++e)
  {
    EventRecord &r = ev[e];
    // structural checks: exactly one info row and one gpoint row per collision
    if (r.nInfoRows != 1)
    {
      warn["ntp_info rows != 1 for an event"]++;
      S.errors++;
    }
    if (r.nGpointRows != 1)
    {
      warn["ntp_gpoint rows != 1 for an event"]++;
      S.errors++;
    }
    if (r.ntrkInfo != r.nAll)
    {
      warn["ntp_info.ntrk != number of ntp_track rows"]++;
      S.errors++;
    }
    if (r.gntrkallInfo != r.NprimAll) warn["ntp_info.gntrkall != number of G4 primaries in ntp_gtrack (info)"]++;
    if (!Finite(r.zTrue)) warn["non-finite true vertex z"]++;

    // ---- truth classification of the good tracks (QA only; n unchanged) ----
    // Uniqueness: the FIRST good track (in ntp_track order, i.e. SvtxTrackMap
    // key order, which is unrelated to track quality) matched to a given G4
    // particle is the "unique match"; every further good track matched to the
    // same particle in the same event is an "additional reco match".
    // Unmatched tracks never enter the duplicate test.
    // Exclusive class (priority): 1 unmatched > 2 additional match >
    //   3 secondary > 4 primary outside N acceptance > 5 primary inside N acceptance.
    std::map<long long, std::vector<int>> byParticle;  // G4 id -> indices of its good tracks
    std::vector<int> trkClass(r.goodTracks.size(), 0);
    for (size_t i = 0; i < r.goodTracks.size(); ++i)
    {
      const auto &g = r.goodTracks[i];
      int origin = 3, uniq = 3;  // unmatched
      if (g.matched)
      {
        origin = g.prim ? 1 : 2;
        auto &v = byParticle[g.id];
        uniq = v.empty() ? 1 : 2;
        v.push_back(static_cast<int>(i));
        if (uniq == 1)
        {
          if (g.prim) r.nUniquePrim++;
          else r.nUniqueSec++;
          if (g.primInAcc) r.nUniquePrimAcc++;
          if (g.prim && !g.primInAcc) r.nPrimOutAcc++;
        }
        else
        {
          if (g.prim) r.nAddPrim++;
          else r.nAddSec++;
        }
      }
      H.hOriginUnique->Fill(origin, uniq);
      int cls;
      if (!g.matched) cls = 1;
      else if (uniq == 2) cls = 2;
      else if (!g.prim) cls = 3;
      else if (!g.primInAcc) cls = 4;
      else cls = 5;
      trkClass[i] = cls;
      H.hTrkClass_all->Fill(cls);
      if (r.N == 0) H.hTrkClass_N0->Fill(cls);
    }
    // consistency with the running counters of step 4
    if (r.nUniquePrim + r.nAddPrim != r.nGoodPrim || r.nUniqueSec + r.nAddSec != r.nGoodSec ||
        r.nAddPrim + r.nAddSec != r.nGoodDup)
      warn["uniqueness classification inconsistent with step-4 counters"]++;
    // each true primary in the N acceptance can have at most one UNIQUE match
    if (r.nUniquePrimAcc > r.N)
    {
      warn["more unique in-acceptance primary matches than N"]++;
      H.hStats->Fill("uniquePrimAcc>N", 1);
    }
    // duplicate diagnostics per matched particle
    for (const auto &kv : byParticle)
    {
      const auto &idx = kv.second;
      H.hMatchesPerParticle->Fill(std::min<int>(idx.size(), 10));
      for (size_t j = 0; j < idx.size(); ++j)
      {
        const auto &g = r.goodTracks[idx[j]];
        const int cat = (idx.size() == 1) ? 0 : (j == 0 ? 1 : 2);
        if (g.nhits > 0) H.hPurity[cat]->Fill(g.nfromtruth / g.nhits);
        H.hNtpcDup[cat]->Fill(g.ntpc);
        H.hQualDup[cat]->Fill(g.quality);
        H.hMvtxDup[cat]->Fill(g.nmvtx);
        if (cat == 2) H.hAddMinusFirst_nfromtruth->Fill(g.nfromtruth - r.goodTracks[idx[0]].nfromtruth);
      }
    }

    const bool vtx = r.HasVertex();
    S.nVertexed += vtx;
    S.nVertexless += !vtx;
    S.nGood += r.nGood;
    if (!vtx)
    {
      S.nGoodInVertexless += r.nGood;
      S.nTrkInVertexless += r.nAll;
    }
    S.nTrkWithVtxId += r.nTrkWithVtxId;

    const int n = Clamp(r.nGood, M, nOverflow);
    const int N = Clamp(r.N, M, nOverflow);

    // statistics
    H.hStats->Fill("events", 1);
    H.hStats->Fill(vtx ? "vertexed" : "vertexless", 1);
    if (r.nAll == 0) H.hStats->Fill("zero tracks", 1);
    if (r.nGood == 0) H.hStats->Fill("zero good tracks", 1);
    if (r.nGood <= 5) H.hStats->Fill(Form("n=%d", r.nGood), 1);
    if (r.N == 0) H.hStats->Fill("N=0", 1);
    if (r.N == 0 && r.nGood > 0) H.hStats->Fill("N=0 & n>0", 1);
    if (r.nGood > r.N) H.hStats->Fill("n>N", 1);
    if (r.nGood < r.N) H.hStats->Fill("n<N", 1);
    if (r.nGood == r.N) H.hStats->Fill("n=N", 1);
    if (r.nGoodSec > 0) H.hStats->Fill("has secondary", 1);
    if (r.nGoodDup > 0) H.hStats->Fill("has duplicate", 1);
    if (r.nGoodUnm > 0) H.hStats->Fill("has unmatched", 1);
    if (r.nVtx > 1) H.hStats->Fill("nVtx>1", 1);
    H.hStats->Fill("tracks", r.nAll);
    H.hStats->Fill("good tracks", r.nGood);
    H.hStats->Fill("G4 primaries", r.NprimAll);

    // multiplicities: the SAME binning for numerator (vtx) and denominator (all)
    H.hn_all->Fill(n);
    H.hnW_all->Fill(n);
    H.hN_all->Fill(N);
    H.hNW_all->Fill(N);
    H.hResp_all->Fill(N, n);
    H.hnMinusN_all->Fill(r.nGood - r.N);
    if (r.N > 0) H.hnOverN_all->Fill(double(r.nGood) / r.N);
    H.hnAll_all->Fill(std::min(r.nAll, 3 * M));
    (vtx ? H.hnAll_vtx : H.hnAll_novtx)->Fill(std::min(r.nAll, 3 * M));
    H.hnGoodPrim_all->Fill(Clamp(r.nGoodPrim, M, nOverflow));
    H.hnGoodSec_all->Fill(Clamp(r.nGoodSec, M, nOverflow));
    H.hnGoodMVTX_all->Fill(Clamp(r.nGoodMVTX, M, nOverflow));
    H.hnGoodNoMVTX_all->Fill(Clamp(r.nGoodNoMVTX, M, nOverflow));
    H.hGoodPrim_vs_n->Fill(n, Clamp(r.nGoodPrim, M, nOverflow));
    H.hNvtx->Fill(r.nVtx);
    if (Finite(r.zTrue)) H.hZtrue_all->Fill(r.zTrue);
    H.hNpip->Fill(r.Npip);
    H.hNpim->Fill(r.Npim);
    H.hNKp->Fill(r.NKp);
    H.hNKm->Fill(r.NKm);
    H.hNp->Fill(r.Np);
    H.hNpbar->Fill(r.Npbar);

    // ---- silicon-content diagnostics (event counts; n unchanged) ----
    {
      const int m1 = Clamp(r.nGoodMVTX, 30, nOverflow), m2 = Clamp(r.nMVTX2, 30, nOverflow),
                m3 = Clamp(r.nMVTX3, 30, nOverflow);
      H.hnM1_all->Fill(m1);
      H.hnM2_all->Fill(m2);
      H.hnM3_all->Fill(m3);
      H.hEff2D_all->Fill(std::min(n, 40), std::min(r.nMVTX2, 10));
      TH1D *ca = r.nGoodMVTX == 0 ? H.hn_all_m0 : H.hn_all_m1;
      ca->Fill(n);
      if (r.nGoodMVTX >= 2) H.hn_all_m2->Fill(n);
      H.hStats->Fill(r.nGoodMVTX == 0 ? "n_MVTX1=0" : "n_MVTX1>=1", 1);
      if (r.nGoodMVTX >= 2) H.hStats->Fill("n_MVTX1>=2", 1);
      if (vtx)
      {
        H.hnM1_vtx->Fill(m1);
        H.hnM2_vtx->Fill(m2);
        H.hnM3_vtx->Fill(m3);
        H.hEff2D_vtx->Fill(std::min(n, 40), std::min(r.nMVTX2, 10));
        (r.nGoodMVTX == 0 ? H.hn_vtx_m0 : H.hn_vtx_m1)->Fill(n);
        if (r.nGoodMVTX >= 2) H.hn_vtx_m2->Fill(n);
      }
    }
    // ---- truth-category multiplicities and diagnostic responses ----
    H.hnUniquePrim_all->Fill(Clamp(r.nUniquePrim, M, nOverflow));
    H.hnUniqueSec_all->Fill(Clamp(r.nUniqueSec, M, nOverflow));
    H.hnAdd_all->Fill(Clamp(r.nAddPrim + r.nAddSec, M, nOverflow));
    H.hnUnm_all->Fill(Clamp(r.nGoodUnm, M, nOverflow));
    H.hResp_uprim->Fill(N, Clamp(r.nUniquePrim, M, nOverflow));
    H.hResp_uprimAcc->Fill(N, Clamp(r.nUniquePrimAcc, M, nOverflow));
    if (r.nAddPrim + r.nAddSec > 0) H.hStats->Fill("has additional match", 1);

    // central-z subsample (same binning, same definitions)
    const bool zc = Finite(r.zTrue) && std::fabs(r.zTrue) < cfg.zcut_abs;
    if (zc)
    {
      H.hn_all_zc->Fill(n);
      H.hnW_all_zc->Fill(n);
      H.hN_all_zc->Fill(N);
      H.hNW_all_zc->Fill(N);
      if (vtx)
      {
        H.hn_vtx_zc->Fill(n);
        H.hnW_vtx_zc->Fill(n);
        H.hN_vtx_zc->Fill(N);
        H.hNW_vtx_zc->Fill(N);
      }
    }

    if (vtx)
    {
      H.hn_vtx->Fill(n);
      H.hnW_vtx->Fill(n);
      H.hN_vtx->Fill(N);
      H.hNW_vtx->Fill(N);
      H.hResp_vtx->Fill(N, n);
      H.hnMinusN_vtx->Fill(r.nGood - r.N);
      if (Finite(r.zTrue)) H.hZtrue_vtx->Fill(r.zTrue);
      if (Finite(r.zTrue) && Finite(r.zReco))
      {
        H.hDz->Fill(r.zReco - r.zTrue);
        H.hDz_vs_n->Fill(n, r.zReco - r.zTrue);
      }
      if (Finite(r.xTrue) && Finite(r.xReco)) H.hDx->Fill(r.xReco - r.xTrue);
      if (Finite(r.yTrue) && Finite(r.yReco)) H.hDy->Fill(r.yReco - r.yTrue);
      if (Finite(r.zTrue))
      {
        if (Finite(r.zReco)) H.hDz_vs_z->Fill(r.zTrue, r.zReco - r.zTrue);
        if (Finite(r.xReco) && Finite(r.xTrue)) H.hDx_vs_z->Fill(r.zTrue, r.xReco - r.xTrue);
        if (Finite(r.yReco) && Finite(r.yTrue)) H.hDy_vs_z->Fill(r.zTrue, r.yReco - r.yTrue);
      }
      if (Finite(r.xTrue)) H.hXtrue_vtx->Fill(r.xTrue);
      if (Finite(r.yTrue)) H.hYtrue_vtx->Fill(r.yTrue);
      if (Finite(r.xReco)) H.hXreco->Fill(r.xReco);
      if (Finite(r.yReco)) H.hYreco->Fill(r.yReco);
      H.hVtxNtracks->Fill(r.vtxNtracks);
    }
    else
    {
      H.hn_novtx->Fill(n);
      if (Finite(r.zTrue)) H.hZtrue_novtx->Fill(r.zTrue);
    }

    // composition of good tracks: unique primaries, unique secondaries,
    // duplicates (any origin), unmatched. These four add up to nGood.
    const int k = vtx ? 0 : 1;
    const int nDupSec = r.nGoodDup - r.nGoodDupPrim;
    const int cPrim = r.nGoodPrim - r.nGoodDupPrim;
    const int cSec = r.nGoodSec - nDupSec;
    const int cnt[4] = {cPrim, cSec, r.nGoodDup, r.nGoodUnm};
    for (int b = 0; b < 4; ++b)
    {
      if (cnt[b] <= 0) continue;
      H.hComp_good[k]->Fill(b + 1, cnt[b]);
      H.hComp_vs_n->Fill(n, b + 1, cnt[b]);
    }
    if (cPrim + cSec + r.nGoodDup + r.nGoodUnm != r.nGood) warn["composition does not add up to n"]++;
    const int nInttOnly = r.nGoodNoMVTX - r.nGoodTPCOnly;
    const int sil[3] = {r.nGoodMVTX, nInttOnly, r.nGoodTPCOnly};
    for (int b = 0; b < 3; ++b)
    {
      if (sil[b] <= 0) continue;
      H.hSilicon_good[k]->Fill(b + 1, sil[b]);
      H.hSilicon_vs_n->Fill(n, b + 1, sil[b]);
    }

    if (evtTree)
    {
      t_job = S.job;
      t_event = static_cast<int>(e);
      t_n = r.nGood;
      t_nAll = r.nAll;
      t_nGoodPrim = r.nGoodPrim;
      t_nGoodSec = r.nGoodSec;
      t_nGoodUnm = r.nGoodUnm;
      t_nGoodDup = r.nGoodDup;
      t_nGoodMVTX = r.nGoodMVTX;
      t_nGoodNoMVTX = r.nGoodNoMVTX;
      t_nGoodTPCOnly = r.nGoodTPCOnly;
      t_N = r.N;
      t_NchAll = r.NchAll;
      t_Npip = r.Npip;
      t_Npim = r.Npim;
      t_NKp = r.NKp;
      t_NKm = r.NKm;
      t_Np = r.Np;
      t_Npbar = r.Npbar;
      t_hasVtx = vtx;
      t_nVtx = r.nVtx;
      t_vtxNtracks = r.vtxNtracks;
      t_zTrue = r.zTrue;
      t_xTrue = r.xTrue;
      t_yTrue = r.yTrue;
      t_zReco = r.zReco;
      t_xReco = r.xReco;
      t_yReco = r.yReco;
      t_nMVTX2 = r.nMVTX2;
      t_nMVTX3 = r.nMVTX3;
      t_nINTTonly = r.nGoodNoMVTX - r.nGoodTPCOnly;
      t_nUniquePrim = r.nUniquePrim;
      t_nUniqueSec = r.nUniqueSec;
      t_nAddPrim = r.nAddPrim;
      t_nAddSec = r.nAddSec;
      t_nUniquePrimAcc = r.nUniquePrimAcc;
      t_nPrimOutAcc = r.nPrimOutAcc;
      evtTree->Fill();
    }
  }
  if (nOverflow) warn["multiplicity above mult_max (clamped into last bin)"] = nOverflow;

  sw.Stop();
  S.seconds = sw.RealTime();
  return S;
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
static void Usage()
{
  std::cout << "Usage: VtxEffAnalysis -c config.cfg [-l file.list] [-o output.root] [-n maxEventsPerFile] [file1.root ...]\n"
               "  -c  configuration file (required)\n"
               "  -l  overrides input_list from the config\n"
               "  -o  overrides output_file from the config\n"
               "  -n  process only the first N events of each file (test runs; 0 = all)\n"
               "  extra arguments are input files (override the list)\n";
}

int main(int argc, char **argv)
{
  std::string cfgPath, listOverride, outOverride;
  long maxEvOverride = -1;
  std::vector<std::string> files;
  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "-c" && i + 1 < argc) cfgPath = argv[++i];
    else if (a == "-l" && i + 1 < argc) listOverride = argv[++i];
    else if (a == "-o" && i + 1 < argc) outOverride = argv[++i];
    else if (a == "-n" && i + 1 < argc) maxEvOverride = std::stol(argv[++i]);
    else if (a == "-h" || a == "--help")
    {
      Usage();
      return 0;
    }
    else files.push_back(a);
  }
  if (cfgPath.empty())
  {
    Usage();
    return 1;
  }
  Config cfg = ReadConfig(cfgPath);
  if (!listOverride.empty()) cfg.input_list = listOverride;
  if (!outOverride.empty()) cfg.output_file = outOverride;
  if (maxEvOverride >= 0) cfg.max_events_per_file = maxEvOverride;

  // Build the input file list (absolute paths).
  if (files.empty())
  {
    std::ifstream in(cfg.input_list);
    if (!in)
    {
      std::cerr << "FATAL: cannot open input list " << cfg.input_list << std::endl;
      return 1;
    }
    std::string line;
    while (std::getline(in, line))
    {
      const auto hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      line = Trim(line);
      if (line.empty()) continue;
      files.push_back(line[0] == '/' ? line : cfg.input_dir + "/" + line);
    }
  }
  if (files.empty())
  {
    std::cerr << "FATAL: no input files" << std::endl;
    return 1;
  }

  // Refuse to write anywhere inside the input directory (read-only data).
  if (!cfg.input_dir.empty() && cfg.output_file.rfind(cfg.input_dir, 0) == 0)
  {
    std::cerr << "FATAL: output file is inside input_dir - refusing" << std::endl;
    return 1;
  }

  std::cout << "=== VtxEffAnalysis ===\n"
            << "config      : " << cfgPath << "\n"
            << "output      : " << cfg.output_file << "\n"
            << "files       : " << files.size() << "\n"
            << "good track  : pt > " << cfg.good_pt_min << " GeV/c, |eta| < " << cfg.good_abseta_max
            << ", ntpc >= " << cfg.good_ntpc_min
            << (cfg.good_nmvtx_min > 0 ? ", nmvtx >= " + std::to_string(cfg.good_nmvtx_min) : std::string()) << "\n"
            << "true N      : G4 primary, charged, pt > " << cfg.true_pt_min << " GeV/c, |eta| < "
            << cfg.true_abseta_max << "\n"
            << "z window    : additional histograms for |z_true| < " << cfg.zcut_abs << " cm\n"
            << "events/file : " << (cfg.max_events_per_file > 0 ? std::to_string(cfg.max_events_per_file) : std::string("all")) << "\n";

  // Histograms must not be attached to the input files.
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  TFile out(cfg.output_file.c_str(), "RECREATE");
  if (out.IsZombie())
  {
    std::cerr << "FATAL: cannot create " << cfg.output_file << std::endl;
    return 1;
  }
  Histos H = BookHistos(cfg);
  TTree *evtTree = nullptr;
  if (cfg.write_event_tree)
  {
    out.cd();
    evtTree = new TTree("evt", "One row per simulated collision (job, event) - Vertex efficiency study");
  }

  ChargeLookup chg;
  std::map<std::string, long> warn;
  std::set<int> jobsSeen;
  std::vector<FileSummary> sums;
  int totalErrors = 0;
  for (const auto &fpath : files)
  {
    const int job = JobFromName(fpath);
    if (job < 0)
    {
      std::cerr << "ERROR: cannot parse job number from " << fpath << " - skipped" << std::endl;
      totalErrors++;
      continue;
    }
    if (!jobsSeen.insert(job).second)
    {
      // The same job twice would double count its events.
      std::cerr << "ERROR: job " << job << " appears twice in the input - skipped" << std::endl;
      totalErrors++;
      continue;
    }
    FileSummary S = ProcessFile(fpath, cfg, H, evtTree, chg, warn);
    totalErrors += S.errors;
    sums.push_back(S);
    std::cout << Form("[file] job=%06d events=%ld vertexed=%ld vertexless=%ld tracks=%ld good=%ld "
                      "goodInVertexless=%ld tracksInVertexless=%ld vtxRowsFinite=%ld vtxRowsNaN=%ld "
                      "tracksWithVertexID>=0=%ld trackNaN=%ld G4rows=%ld primaries=%ld unknownPDG=%ld errors=%d time=%.1fs",
                      S.job, S.nEvents, S.nVertexed, S.nVertexless, S.nTracks, S.nGood, S.nGoodInVertexless,
                      S.nTrkInVertexless, S.nVtxRowsFinite, S.nVtxRowsNaN, S.nTrkWithVtxId, S.nTrackNaN,
                      S.nGtrackRows, S.nPrimaries, S.nUnknownPdg, S.errors, S.seconds)
              << std::endl;
  }

  // Store the configuration with the histograms for provenance.
  out.cd();
  // " nmvtx>=K" appears only if an MVTX requirement is set, so the baseline
  // string is unchanged.
  const std::string mvtxCut = cfg.good_nmvtx_min > 0 ? Form(" nmvtx>=%d", cfg.good_nmvtx_min) : "";
  TNamed cfgObj("config", Form("good: pt>%g |eta|<%g ntpc>=%d%s ; N: G4 primary charged pt>%g |eta|<%g ; "
                               "zcut |z_true|<%g cm ; files=%zu",
                               cfg.good_pt_min, cfg.good_abseta_max, cfg.good_ntpc_min, mvtxCut.c_str(),
                               cfg.true_pt_min, cfg.true_abseta_max, cfg.zcut_abs, sums.size()));
  cfgObj.Write();
  TH1 *all[] = {H.hStats, H.hn_all, H.hn_vtx, H.hn_novtx, H.hN_all, H.hN_vtx, H.hnW_all, H.hnW_vtx, H.hNW_all, H.hNW_vtx,
                H.hn_all_zc, H.hn_vtx_zc, H.hN_all_zc, H.hN_vtx_zc, H.hnW_all_zc, H.hnW_vtx_zc, H.hNW_all_zc,
                H.hNW_vtx_zc, H.hnAll_all, H.hnAll_vtx, H.hnAll_novtx, H.hnGoodPrim_all, H.hnGoodSec_all, H.hnGoodMVTX_all,
                H.hnGoodNoMVTX_all, H.hZtrue_all, H.hZtrue_vtx, H.hZtrue_novtx, H.hNvtx, H.hDz, H.hDx, H.hDy,
                H.hVtxNtracks, H.hDz_vs_n, H.hResp_all, H.hResp_vtx, H.hnMinusN_all, H.hnMinusN_vtx, H.hnOverN_all,
                H.hNpip, H.hNpim, H.hNKp, H.hNKm, H.hNp, H.hNpbar, H.hComp_vs_n, H.hSilicon_vs_n, H.hGoodPrim_vs_n,
                H.hnM1_all, H.hnM1_vtx, H.hnM2_all, H.hnM2_vtx, H.hnM3_all, H.hnM3_vtx, H.hEff2D_all, H.hEff2D_vtx,
                H.hn_all_m0, H.hn_vtx_m0, H.hn_all_m1, H.hn_vtx_m1, H.hn_all_m2, H.hn_vtx_m2, H.hOriginUnique,
                H.hTrkClass_all, H.hTrkClass_N0, H.hnUniquePrim_all, H.hnUniqueSec_all, H.hnAdd_all, H.hnUnm_all,
                H.hMatchesPerParticle, H.hAddMinusFirst_nfromtruth, H.hResp_uprim, H.hResp_uprimAcc, H.hDz_vs_z,
                H.hDx_vs_z, H.hDy_vs_z, H.hXtrue_vtx, H.hYtrue_vtx, H.hXreco, H.hYreco,
                H.hPurity[0], H.hPurity[1], H.hPurity[2], H.hNtpcDup[0], H.hNtpcDup[1], H.hNtpcDup[2],
                H.hQualDup[0], H.hQualDup[1], H.hQualDup[2], H.hMvtxDup[0], H.hMvtxDup[1], H.hMvtxDup[2]};
  for (TH1 *h : all) h->Write();
  for (int k = 0; k < 2; ++k)
  {
    TH1 *pk[] = {H.hPt_all[k], H.hPt_good[k], H.hEta_all[k], H.hEta_good[k], H.hNtpc_all[k], H.hNtpc_good[k],
                 H.hNmvtx_all[k], H.hNmvtx_good[k], H.hNintt_all[k], H.hNintt_good[k], H.hComp_good[k],
                 H.hSilicon_good[k]};
    for (TH1 *h : pk) h->Write();
  }
  if (evtTree) evtTree->Write();
  out.Close();

  std::cout << "--- warnings / diagnostics ---" << std::endl;
  if (warn.empty()) std::cout << "(none)" << std::endl;
  for (const auto &w : warn) std::cout << "  " << w.second << " x " << w.first << std::endl;
  std::cout << "structural errors: " << totalErrors << std::endl;
  std::cout << "wrote " << cfg.output_file << std::endl;
  return totalErrors == 0 ? 0 : 3;
}
