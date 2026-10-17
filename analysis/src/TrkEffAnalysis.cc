// =============================================================================
// TrkEffAnalysis.cc - single-particle tracking efficiency from EvalSkim files
//
// Reads the EventTree, TruthTree and RecoTree written by EvalSkim (step 1) and
// fills, for each particle species, raw-count histograms of TRUE rapidity y
// and TRUE transverse momentum pT:
//
//   denominator  D(y,pT) = selected true particles
//   numerator    N(y,pT) = selected true particles with >= 1 selected
//                          reconstructed track matched to them
//   numerator    NQ(y,pT)= as N, but at least one of those tracks has the
//                          correct reconstructed charge
//
//   efficiency   eps(y,pT)  = N / D        (any reconstructed charge)
//                epsQ(y,pT) = NQ / D       (correct reconstructed charge)
//
// Definitions (all cut values come from the configuration file):
//
//   Event selection
//     use_reco_zvtx = 0 : |z_true| < zvtx_max            (all collisions count)
//     use_reco_zvtx = 1 : event has a reconstructed vertex and
//                         |z_reco| < zvtx_max            (z_reco = vertex with
//                                                         the most tracks)
//     The vertex-finding efficiency (VtxEffAnalysis) is NOT affected by this
//     setting; it has its own denominator and is computed elsewhere.
//
//   True particle (denominator)
//     G4 primary (particle handed to GEANT4 by PYTHIA; weak-decay daughters,
//     conversion and interaction products are G4 secondaries and never enter),
//     PDG code of one of the configured species, in a selected event,
//     pT_true > truth_pt_min and, if the truth-eta cut is on, |eta_true| <
//     threshold. y and pT are the TRUE values (y with the PDG mass,
//     computed in EvalSkim).
//     Truth-eta cut (v3, Step 6A):
//       use_truth_eta_cut = 1 : |eta_true| < truth_eta_cut_max (default 1.1)
//       use_truth_eta_cut = 0 : no truth-eta cut (the (y, pT) binning defines
//                               the acceptance; final analysis)
//       key absent ("auto")   : legacy rule of v1/v2 configurations - the cut
//                               is on iff truth_abseta_max > 0 (threshold =
//                               truth_abseta_max). Keeps every validated
//                               configuration reproducible.
//     Setting use_truth_eta_cut explicitly together with truth_abseta_max > 0
//     is refused (two controls for one cut). The cut is applied once, in the
//     common truth selection, so the numerator is always a subset of the
//     denominator; the reconstructed-track |eta| cut is independent of it.
//     Two denominator definitions (never combined in one output):
//       A "inclusive" (truth_min_tpc_layers = 0, baseline): all such particles;
//       B "trackable" (truth_min_tpc_layers = N > 0): only particles that left
//         truth clusters in >= N TPC layers (TruthTree.gnltpc). Efficiency B is
//         CONDITIONAL on the particle reaching enough TPC layers; losses before
//         the TPC (decays, interactions, stopping) are excluded from it.
//
//   Reconstructed track (numerator)
//     pT > trk_pt_min, |eta| < trk_abseta_max, ntpc >= trk_ntpc_min,
//     nmaps >= trk_nmvtx_min, nintt >= trk_nintt_min (reconstructed
//     quantities), and matched to a truth particle.
//     Optional DCA cuts (v2; off by default = no DCA cut, as in v1):
//       trk_dca_xy_max > 0 : |linXYReco| <= trk_dca_xy_max
//       trk_dca_z_max  > 0 : |linZReco|  <= trk_dca_z_max
//     linXYReco / linZReco (EvalSkim): straight-line DCA from the track state
//     along the reconstructed momentum, at the transverse closest approach to
//     the event vertex (reconstructed vertex with the most tracks). They are
//     defined only in events with a reconstructed vertex, so DCA cuts require
//     use_reco_zvtx = 1; a non-finite DCA fails the cut and is counted.
//     With trk_nintt_min = 0 and both DCA keys <= 0, the result is identical
//     to v1.
//
//   Truth matching
//     A track is matched to the G4 particle that contributed most of its
//     clusters (evaluator ntp_track.gtrackID). The key of a true particle is
//     (job, event, gtrackID): G4 track ids restart in every event.
//     A true particle enters the numerator at most ONCE, however many selected
//     tracks are matched to it (duplicates are counted separately as
//     information). Only tracks matched to the particle ITSELF count: a pion
//     whose track is matched to its decay muon is not found (baseline).
//
// Output (one ROOT file; all histograms are raw counts, so outputs of several
// runs can be merged with hadd; TEfficiency objects merge correctly too):
//   hDen_<sp>_ypt, hNum_<sp>_ypt, hNumQ_<sp>_ypt   TH2D, x = y_true, y = pT_true [GeV/c]
//   hDen_<sp>_pt,  hNum_<sp>_pt,  hNumQ_<sp>_pt    projections on pT (y inside the y range)
//   hDen_<sp>_y,   hNum_<sp>_y,   hNumQ_<sp>_y     projections on y (pT inside the pT range)
//   eff_<sp>_{ypt,pt,y}, effQ_<sp>_{ypt,pt,y}      TEfficiency (Clopper-Pearson, 68.3%)
//   hTrkEffStats                                   bookkeeping counters
//   config, definitions, skim_files                TNamed (provenance)
//   <sp> = pip, pim, kp, km, p, pbar.
// Bins whose denominator is 0 have an UNDEFINED efficiency (TEfficiency then
// returns 0 - always check the total, e.g. TEfficiency::GetTotalHistogram).
//
// Usage: TrkEffAnalysis -c trkeff.cfg -o out.root [-t summary.txt] [-l skim.list | skim files...]
// Build: analysis/build.sh
// =============================================================================

#include "EffTypes.h"
#include "EvalCommon.h"

#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <unordered_map>
#include <unordered_set>

using namespace evalcommon;

namespace
{
const char *kVersion = "TrkEffAnalysis v4 (2026-10-16)";  // v2: optional INTT and DCA cuts; v3: use_truth_eta_cut switch; v4: species from EffTypes.h

// -----------------------------------------------------------------------------
// Particle species
// -----------------------------------------------------------------------------
struct Species
{
  int pdg;
  const char *tag;    // used in histogram names
  const char *label;  // ROOT LaTeX label for titles
  int charge;         // electric charge in units of e
};

// the identified species of EffTypes.h (kTrack)
const std::vector<Species> &KnownSpecies()
{
  static const std::vector<Species> v = []
  {
    std::vector<Species> out;
    for (const auto *s : efftypes::MeasuredSpecies(false)) out.push_back({s->pdg, s->tag, s->label, s->charge});
    return out;
  }();
  return v;
}

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------
struct TrkEffConfig
{
  // event selection
  bool use_reco_zvtx = false;
  double zvtx_max = 8.0;  // cm
  // true particles
  double truth_pt_min = 0.1;      // GeV/c
  double truth_abseta_max = 1.1;  // EFFECTIVE truth |eta| threshold used by the selection; <= 0: no cut
  std::string truthEtaMode = "auto";  // "auto" (legacy key truth_abseta_max), "0" or "1" (use_truth_eta_cut)
  double truth_eta_cut_max = 1.1;     // threshold of the switch (use_truth_eta_cut = 1)
  bool TruthEtaExplicit() const { return truthEtaMode != "auto"; }
  int truth_min_tpc_layers = 0;   // 0 = inclusive denominator (A); N > 0 = trackable (B)
  // reconstructed tracks
  double trk_pt_min = 0.1;  // GeV/c
  double trk_abseta_max = 1.1;
  int trk_ntpc_min = 18;
  int trk_nmvtx_min = 0;
  int trk_nintt_min = 0;        // 0 = no INTT requirement
  double trk_dca_xy_max = 0;    // cm, |linXYReco|; <= 0 = no cut
  double trk_dca_z_max = 0;     // cm, |linZReco|;  <= 0 = no cut
  bool DcaCut() const { return trk_dca_xy_max > 0 || trk_dca_z_max > 0; }
  bool V2Cuts() const { return trk_nintt_min > 0 || DcaCut(); }
  // binning (true y, true pT)
  std::vector<double> y_edges, pt_edges;
  std::vector<Species> species;
  std::string dump;  // all keys, stored in the output
};

void RequireIncreasing(const std::vector<double> &edges, const char *key)
{
  if (edges.size() < 2)
  {
    std::cerr << "FATAL: " << key << " needs at least two edges" << std::endl;
    std::exit(1);
  }
  for (size_t i = 1; i < edges.size(); ++i)
    if (!(edges[i] > edges[i - 1]))
    {
      std::cerr << "FATAL: " << key << " must be strictly increasing" << std::endl;
      std::exit(1);
    }
}

// overrides: "key=value" strings from the command line (-s), applied after the file
TrkEffConfig ReadTrkEffConfig(const std::string &path, const std::vector<std::string> &overrides)
{
  KeyValueConfig k;
  k.Declare("use_reco_zvtx", "0");
  k.Declare("zvtx_max", "8.0");
  k.Declare("truth_pt_min", "0.1");
  k.Declare("truth_abseta_max", "1.1");    // legacy control (used only when use_truth_eta_cut is absent)
  k.Declare("use_truth_eta_cut", "auto");  // v3 switch: 0 / 1; "auto" = legacy rule
  k.Declare("truth_eta_cut_max", "1.1");   // v3 threshold for use_truth_eta_cut = 1
  k.Declare("truth_min_tpc_layers", "0");
  k.Declare("trk_pt_min", "0.1");
  k.Declare("trk_abseta_max", "1.1");
  k.Declare("trk_ntpc_min", "18");
  k.Declare("trk_nmvtx_min", "0");
  k.Declare("trk_nintt_min", "0");
  k.Declare("trk_dca_xy_max", "0");
  k.Declare("trk_dca_z_max", "0");
  k.Declare("y_edges", "");
  k.Declare("pt_edges", "");
  k.Declare("species", "211, -211, 321, -321, 2212, -2212");
  k.Read(path);
  for (const auto &o : overrides)
  {
    const auto eq = o.find('=');
    if (eq == std::string::npos)
    {
      std::cerr << "FATAL: -s expects key=value, got '" << o << "'" << std::endl;
      std::exit(1);
    }
    k.Set(Trim(o.substr(0, eq)), Trim(o.substr(eq + 1)));  // unknown keys are fatal
  }

  TrkEffConfig c;
  c.use_reco_zvtx = k.Bool("use_reco_zvtx");
  c.zvtx_max = k.Dbl("zvtx_max");
  c.truth_pt_min = k.Dbl("truth_pt_min");
  // truth-eta cut: the switch (v3) or the legacy key - never both
  c.truthEtaMode = Trim(k.Str("use_truth_eta_cut"));
  c.truth_eta_cut_max = k.Dbl("truth_eta_cut_max");
  const double legacyEta = k.Dbl("truth_abseta_max");
  if (c.truthEtaMode == "auto")
    c.truth_abseta_max = legacyEta;  // v1/v2 rule: > 0 = cut on
  else if (c.truthEtaMode == "0" || c.truthEtaMode == "1")
  {
    if (legacyEta > 0)
    {
      std::cerr << "FATAL: use_truth_eta_cut is set explicitly, so truth_abseta_max must be <= 0 (it is " << legacyEta
                << "); use truth_eta_cut_max for the threshold" << std::endl;
      std::exit(1);
    }
    if (c.truthEtaMode == "1" && !(c.truth_eta_cut_max > 0))
    {
      std::cerr << "FATAL: use_truth_eta_cut = 1 needs truth_eta_cut_max > 0" << std::endl;
      std::exit(1);
    }
    c.truth_abseta_max = c.truthEtaMode == "1" ? c.truth_eta_cut_max : 0.0;
  }
  else
  {
    std::cerr << "FATAL: use_truth_eta_cut must be 0 or 1 (or absent), got '" << c.truthEtaMode << "'" << std::endl;
    std::exit(1);
  }
  c.truth_min_tpc_layers = k.Int("truth_min_tpc_layers");
  c.trk_pt_min = k.Dbl("trk_pt_min");
  c.trk_abseta_max = k.Dbl("trk_abseta_max");
  c.trk_ntpc_min = k.Int("trk_ntpc_min");
  c.trk_nmvtx_min = k.Int("trk_nmvtx_min");
  c.trk_nintt_min = k.Int("trk_nintt_min");
  c.trk_dca_xy_max = k.Dbl("trk_dca_xy_max");
  c.trk_dca_z_max = k.Dbl("trk_dca_z_max");
  if (c.trk_nintt_min < 0)
  {
    std::cerr << "FATAL: trk_nintt_min must be >= 0" << std::endl;
    std::exit(1);
  }
  if (c.DcaCut() && !c.use_reco_zvtx)
  {
    std::cerr << "FATAL: DCA cuts (trk_dca_xy_max / trk_dca_z_max > 0) need use_reco_zvtx = 1: the DCA is "
                 "computed to the reconstructed event vertex and is undefined in events without one"
              << std::endl;
    std::exit(1);
  }
  c.y_edges = ParseNumberList(k.Str("y_edges"), "y_edges");
  c.pt_edges = ParseNumberList(k.Str("pt_edges"), "pt_edges");
  RequireIncreasing(c.y_edges, "y_edges");
  RequireIncreasing(c.pt_edges, "pt_edges");
  if (c.zvtx_max <= 0)
  {
    std::cerr << "FATAL: zvtx_max must be > 0" << std::endl;
    std::exit(1);
  }
  for (double pdg : ParseNumberList(k.Str("species"), "species"))
  {
    bool found = false;
    for (const auto &s : KnownSpecies())
      if (s.pdg == static_cast<int>(std::lround(pdg)))
      {
        c.species.push_back(s);
        found = true;
      }
    if (!found)
    {
      std::cerr << "FATAL: species " << pdg << " not supported (supported:";
      for (const auto &s : KnownSpecies()) std::cerr << " " << s.pdg;
      std::cerr << ")" << std::endl;
      std::exit(1);
    }
  }
  c.dump = k.Dump();
  return c;
}

// -----------------------------------------------------------------------------
// Selections. Each takes the values it needs, so the physics definition is in
// one place and can be read without the I/O code around it.
// -----------------------------------------------------------------------------
bool PassEventSelection(const TrkEffConfig &c, float zTrue, short hasVtx, float zReco)
{
  if (c.use_reco_zvtx) return hasVtx != 0 && std::isfinite(zReco) && std::fabs(zReco) < c.zvtx_max;
  return std::isfinite(zTrue) && std::fabs(zTrue) < c.zvtx_max;
}

// Truth acceptance of the denominator: kinematics and, for the trackable
// definition, the number of TPC layers with truth clusters.
bool PassTruthSelection(const TrkEffConfig &c, float pt, float eta, short nTpcLayers)
{
  if (c.truth_min_tpc_layers > 0 && nTpcLayers < c.truth_min_tpc_layers) return false;
  if (!(std::isfinite(pt) && pt > c.truth_pt_min)) return false;
  if (c.truth_abseta_max > 0 && !(std::isfinite(eta) && std::fabs(eta) < c.truth_abseta_max)) return false;
  return true;
}

bool PassTrackSelection(const TrkEffConfig &c, float pt, float eta, short ntpc, short nmaps, short nintt)
{
  return std::isfinite(pt) && std::isfinite(eta) && pt > c.trk_pt_min && std::fabs(eta) < c.trk_abseta_max &&
         ntpc >= c.trk_ntpc_min && nmaps >= c.trk_nmvtx_min && nintt >= c.trk_nintt_min;
}

// DCA cut (inclusive, as in the data analysis); always true when no DCA cut is configured
bool PassDcaSelection(const TrkEffConfig &c, float dcaXY, float dcaZ)
{
  if (c.trk_dca_xy_max > 0 && !(std::isfinite(dcaXY) && std::fabs(dcaXY) <= c.trk_dca_xy_max)) return false;
  if (c.trk_dca_z_max > 0 && !(std::isfinite(dcaZ) && std::fabs(dcaZ) <= c.trk_dca_z_max)) return false;
  return true;
}

// The skim keeps charged G4 primaries only inside its retention window
// (pT >= skim truth_pt_min AND (|eta| < skim truth_abseta_max OR |y| <
// skim truth_absy_max)). The denominator is complete only if every particle
// that this analysis selects lies inside that window. Returns an empty string
// if it does, otherwise the reason.
std::string CheckSkimCoverage(const TrkEffConfig &c, const std::map<std::string, std::string> &skimCfg)
{
  auto get = [&](const char *key) -> double {
    auto it = skimCfg.find(key);
    return it == skimCfg.end() ? NAN : std::stod(it->second);
  };
  const double skimPt = get("truth_pt_min"), skimEta = get("truth_abseta_max"), skimY = get("truth_absy_max");
  if (!std::isfinite(skimPt) || !std::isfinite(skimEta) || !std::isfinite(skimY))
    return "skim file has no readable retention configuration";
  if (c.truth_pt_min < skimPt) return Form("truth_pt_min %g < skim retention %g", c.truth_pt_min, skimPt);
  const bool etaCovered = c.truth_abseta_max > 0 && c.truth_abseta_max <= skimEta;
  const bool yCovered = std::fabs(c.y_edges.front()) <= skimY && std::fabs(c.y_edges.back()) <= skimY;
  if (!etaCovered && !yCovered)
    return Form("truth |eta| cut %g and y range [%g, %g] both exceed the skim retention (|eta| < %g or |y| < %g)",
                c.truth_abseta_max, c.y_edges.front(), c.y_edges.back(), skimEta, skimY);
  return "";
}

// -----------------------------------------------------------------------------
// Histograms of one species
// -----------------------------------------------------------------------------
struct SpeciesHistograms
{
  Species sp;
  TH2D *den = nullptr, *num = nullptr, *numQ = nullptr;

  SpeciesHistograms(const Species &s, const TrkEffConfig &c) : sp(s)
  {
    auto book = [&](const char *kind, const char *what) {
      auto *h = new TH2D(Form("h%s_%s_ypt", kind, s.tag), Form("%s: %s;y_{true};p_{T,true} (GeV/c)", s.label, what),
                         static_cast<int>(c.y_edges.size()) - 1, c.y_edges.data(), static_cast<int>(c.pt_edges.size()) - 1,
                         c.pt_edges.data());
      h->Sumw2();
      return h;
    };
    den = book("Den", "selected true particles");
    num = book("Num", "true particles with a selected matched track");
    numQ = book("NumQ", "true particles with a selected matched track of correct charge");
  }
};

// Per-particle matching result, built from the RecoTree of one file.
struct MatchInfo
{
  int nSelectedTracks = 0;  // selected tracks whose dominant truth particle is this one
  bool correctCharge = false;  // at least one of them has the correct charge
};

struct TruthKey
{
  int job, event, gtrackID;
  bool operator==(const TruthKey &o) const { return job == o.job && event == o.event && gtrackID == o.gtrackID; }
};
struct TruthKeyHash
{
  size_t operator()(const TruthKey &k) const
  {
    return std::hash<long long>()((static_cast<long long>(k.job) << 40) ^ (static_cast<long long>(k.event) << 24) ^
                                  static_cast<long long>(k.gtrackID & 0xffffff));
  }
};

// Bookkeeping counters (also written as hTrkEffStats)
struct Counters
{
  long files = 0, events = 0, eventsSelected = 0;
  long tracksSelected = 0, tracksSelectedMatched = 0, tracksSelectedUnmatched = 0;
  long tracksOnDenominatorParticles = 0;  // selected tracks matched to a denominator particle
  long duplicateTracks = 0;               // selected tracks beyond the first per denominator particle
  long truthSelected = 0, truthFound = 0, truthFoundQ = 0;
  long truthOutsideBinning = 0;  // selected true particles outside the (y, pT) histogram range
  // v2, only with DCA cuts: tracks passing all other cuts, then failing the DCA cut
  long tracksFailDca = 0, tracksDcaNonFinite = 0;
};

// -----------------------------------------------------------------------------
// One skim file: EventTree -> selected events, RecoTree -> matches,
// TruthTree -> fill. Maps are per file (job), so memory does not grow with
// the number of files.
// -----------------------------------------------------------------------------
bool ProcessSkimFile(const std::string &path, const TrkEffConfig &cfg, std::vector<SpeciesHistograms> &H, Counters &C)
{
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return false;
  }
  auto *tEvent = dynamic_cast<TTree *>(f->Get("EventTree"));
  auto *tTruth = dynamic_cast<TTree *>(f->Get("TruthTree"));
  auto *tReco = dynamic_cast<TTree *>(f->Get("RecoTree"));
  auto *skimCfg = dynamic_cast<TNamed *>(f->Get("config"));
  if (!tEvent || !tTruth || !tReco || !skimCfg)
  {
    std::cerr << "ERROR: " << path << " is not an EvalSkim file (missing trees or config)" << std::endl;
    return false;
  }
  const std::string coverage = CheckSkimCoverage(cfg, ParseKeyValueText(skimCfg->GetTitle()));
  if (!coverage.empty())
  {
    std::cerr << "FATAL: incomplete denominator for " << path << ": " << coverage << std::endl;
    std::exit(1);
  }
  C.files++;

  // ---- 1. selected events --------------------------------------------------
  std::unordered_set<long long> selectedEvents;  // key = job * 2^32 + event
  auto eventKey = [](int job, int event) { return (static_cast<long long>(job) << 32) | static_cast<unsigned>(event); };
  {
    int job, event;
    float zTrue, zReco;
    short hasVtx;
    tEvent->SetBranchAddress("job", &job);
    tEvent->SetBranchAddress("event", &event);
    tEvent->SetBranchAddress("zTrue", &zTrue);
    tEvent->SetBranchAddress("zReco", &zReco);
    tEvent->SetBranchAddress("hasVtx", &hasVtx);
    for (long i = 0; i < tEvent->GetEntries(); ++i)
    {
      tEvent->GetEntry(i);
      C.events++;
      if (!PassEventSelection(cfg, zTrue, hasVtx, zReco)) continue;
      selectedEvents.insert(eventKey(job, event));
      C.eventsSelected++;
    }
  }

  // ---- 2. selected reconstructed tracks -> matched truth particles ---------
  std::unordered_map<TruthKey, MatchInfo, TruthKeyHash> matches;
  {
    int job, event, gtrackID;
    float pt, eta, dcaXY = NAN, dcaZ = NAN;
    short ntpc, nmaps, nintt = 0, matched, chargeOK;
    tReco->SetBranchStatus("*", 0);
    for (const char *b : {"job", "event", "gtrackID", "pt", "eta", "ntpc", "nmaps", "matched", "chargeOK"})
      tReco->SetBranchStatus(b, 1);
    tReco->SetBranchAddress("job", &job);
    tReco->SetBranchAddress("event", &event);
    tReco->SetBranchAddress("gtrackID", &gtrackID);
    tReco->SetBranchAddress("pt", &pt);
    tReco->SetBranchAddress("eta", &eta);
    tReco->SetBranchAddress("ntpc", &ntpc);
    tReco->SetBranchAddress("nmaps", &nmaps);
    tReco->SetBranchAddress("matched", &matched);
    tReco->SetBranchAddress("chargeOK", &chargeOK);
    if (cfg.trk_nintt_min > 0)
    {
      tReco->SetBranchStatus("nintt", 1);
      tReco->SetBranchAddress("nintt", &nintt);
    }
    if (cfg.DcaCut())
    {
      for (const char *b : {"linXYReco", "linZReco"}) tReco->SetBranchStatus(b, 1);
      tReco->SetBranchAddress("linXYReco", &dcaXY);
      tReco->SetBranchAddress("linZReco", &dcaZ);
    }
    for (long i = 0; i < tReco->GetEntries(); ++i)
    {
      tReco->GetEntry(i);
      if (!selectedEvents.count(eventKey(job, event))) continue;
      if (!PassTrackSelection(cfg, pt, eta, ntpc, nmaps, nintt)) continue;
      if (!PassDcaSelection(cfg, dcaXY, dcaZ))
      {
        C.tracksFailDca++;
        if (!std::isfinite(dcaXY) || !std::isfinite(dcaZ)) C.tracksDcaNonFinite++;
        continue;
      }
      C.tracksSelected++;
      if (!matched)
      {
        C.tracksSelectedUnmatched++;
        continue;
      }
      C.tracksSelectedMatched++;
      MatchInfo &m = matches[{job, event, gtrackID}];
      m.nSelectedTracks++;
      if (chargeOK == 1) m.correctCharge = true;
    }
  }

  // ---- 3. true particles: denominator and numerators ------------------------
  {
    int job, event, gtrackID, pdg;
    short q, isPrimary;
    float pt, eta, y;
    tTruth->SetBranchStatus("*", 0);
    short gnltpc;
    for (const char *b : {"job", "event", "gtrackID", "pdg", "q", "isPrimary", "pt", "eta", "y", "gnltpc"})
      tTruth->SetBranchStatus(b, 1);
    tTruth->SetBranchAddress("job", &job);
    tTruth->SetBranchAddress("event", &event);
    tTruth->SetBranchAddress("gtrackID", &gtrackID);
    tTruth->SetBranchAddress("pdg", &pdg);
    tTruth->SetBranchAddress("q", &q);
    tTruth->SetBranchAddress("isPrimary", &isPrimary);
    tTruth->SetBranchAddress("gnltpc", &gnltpc);
    tTruth->SetBranchAddress("pt", &pt);
    tTruth->SetBranchAddress("eta", &eta);
    tTruth->SetBranchAddress("y", &y);
    for (long i = 0; i < tTruth->GetEntries(); ++i)
    {
      tTruth->GetEntry(i);
      if (isPrimary != 1) continue;
      SpeciesHistograms *h = nullptr;
      for (auto &s : H)
        if (s.sp.pdg == pdg) h = &s;
      if (!h) continue;
      if (q != h->sp.charge)
      {
        std::cerr << "FATAL: TruthTree charge " << q << " for PDG " << pdg << " in " << path << std::endl;
        std::exit(2);
      }
      if (!selectedEvents.count(eventKey(job, event))) continue;
      if (!PassTruthSelection(cfg, pt, eta, gnltpc)) continue;
      if (!std::isfinite(y))
      {
        std::cerr << "FATAL: non-finite rapidity for PDG " << pdg << " in " << path << std::endl;
        std::exit(2);
      }
      C.truthSelected++;
      if (y < cfg.y_edges.front() || y >= cfg.y_edges.back() || pt < cfg.pt_edges.front() || pt >= cfg.pt_edges.back())
        C.truthOutsideBinning++;
      h->den->Fill(y, pt);
      auto it = matches.find({job, event, gtrackID});
      if (it == matches.end()) continue;
      // counted once, whatever the number of matched selected tracks
      h->num->Fill(y, pt);
      C.truthFound++;
      C.tracksOnDenominatorParticles += it->second.nSelectedTracks;
      C.duplicateTracks += it->second.nSelectedTracks - 1;
      if (it->second.correctCharge)
      {
        h->numQ->Fill(y, pt);
        C.truthFoundQ++;
      }
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
// Output
// -----------------------------------------------------------------------------
// Projection of a (y, pT) histogram onto pT (over the y bins inside the range)
// or onto y (over the pT bins inside the range). Under/overflow of the summed
// axis is excluded, so numerator and denominator projections stay consistent.
TH1D *ProjectPt(const TH2D *h, const char *name)
{
  auto *p = h->ProjectionY(name, 1, h->GetNbinsX(), "e");
  p->SetTitle(Form("%s;p_{T,true} (GeV/c);particles", h->GetTitle()));
  return p;
}
TH1D *ProjectY(const TH2D *h, const char *name)
{
  auto *p = h->ProjectionX(name, 1, h->GetNbinsY(), "e");
  p->SetTitle(Form("%s;y_{true};particles", h->GetTitle()));
  return p;
}

TEfficiency *MakeEfficiency(const TH1 &pass, const TH1 &total, const char *name, const char *title)
{
  auto *e = new TEfficiency(pass, total);
  e->SetName(name);
  e->SetTitle(title);
  e->SetStatisticOption(TEfficiency::kFCP);  // Clopper-Pearson
  e->SetConfidenceLevel(0.683);
  return e;
}

// Text of the v2 DCA cut, e.g. "|DCA_xy| <= 1.5 cm, |DCA_z| <= 1.5 cm (to the reco vertex)"
std::string DcaText(const TrkEffConfig &c)
{
  if (!c.DcaCut()) return "no DCA cut";
  std::string s;
  if (c.trk_dca_xy_max > 0) s += Form("|DCA_xy| <= %g cm", c.trk_dca_xy_max);
  if (c.trk_dca_z_max > 0) s += std::string(s.empty() ? "" : ", ") + Form("|DCA_z| <= %g cm", c.trk_dca_z_max);
  return s + " (to the reco vertex)";
}

// Text summary: integrated counts per species and epsilon(pT); empty
// denominator bins are written as "undefined".
void WriteSummary(const std::string &path, const TrkEffConfig &cfg, const std::vector<SpeciesHistograms> &H,
                  const Counters &C)
{
  FILE *out = path.empty() ? stdout : fopen(path.c_str(), "w");
  if (!out)
  {
    std::cerr << "ERROR: cannot write " << path << std::endl;
    return;
  }
  fprintf(out, "%s\n", kVersion);
  fprintf(out, "event selection: %s < %g cm%s\n", cfg.use_reco_zvtx ? "|z_reco|" : "|z_true|", cfg.zvtx_max,
          cfg.use_reco_zvtx ? " (reconstructed vertex required)" : "");
  fprintf(out, "true particles : G4 primary, pT > %g GeV/c%s\n", cfg.truth_pt_min,
          cfg.truth_abseta_max > 0 ? Form(", |eta| < %g", cfg.truth_abseta_max) : ", no eta cut");
  if (cfg.truth_min_tpc_layers > 0)
    fprintf(out, "denominator    : TRACKABLE - only particles with truth clusters in >= %d TPC layers\n"
                 "                 (efficiency conditional on reaching the TPC)\n",
            cfg.truth_min_tpc_layers);
  if (!cfg.V2Cuts())
    fprintf(out, "tracks         : pT > %g GeV/c, |eta| < %g, ntpc >= %d, nmvtx >= %d, matched (no DCA cut)\n",
            cfg.trk_pt_min, cfg.trk_abseta_max, cfg.trk_ntpc_min, cfg.trk_nmvtx_min);
  else
    fprintf(out, "tracks         : pT > %g GeV/c, |eta| < %g, ntpc >= %d, nmvtx >= %d, nintt >= %d, %s, matched\n",
            cfg.trk_pt_min, cfg.trk_abseta_max, cfg.trk_ntpc_min, cfg.trk_nmvtx_min, cfg.trk_nintt_min,
            DcaText(cfg).c_str());
  fprintf(out, "files %ld, events %ld, selected events %ld\n\n", C.files, C.events, C.eventsSelected);
  fprintf(out, "%-6s %10s %10s %10s %18s %18s %10s\n", "", "D (den)", "N (num)", "NQ", "eps = N/D", "epsQ = NQ/D",
          "out of bins");
  for (const auto &h : H)
  {
    const double d = h.den->GetEntries(), n = h.num->GetEntries(), nq = h.numQ->GetEntries();
    const double inD = h.den->Integral(1, h.den->GetNbinsX(), 1, h.den->GetNbinsY());
    auto ratio = [&](double a) {
      return d > 0 ? std::string(Form("%.3f +- %.3f", a / d, std::sqrt(a / d * (1 - a / d) / d))) : "undefined";
    };
    const std::string eps = ratio(n), epsQ = ratio(nq);
    fprintf(out, "%-6s %10.0f %10.0f %10.0f %18s %18s %10.0f\n", h.sp.tag, d, n, nq, eps.c_str(), epsQ.c_str(),
            d - inD);
  }
  fprintf(out, "(integrated values with a binomial error, for orientation only; use the TEfficiency objects)\n");
  if (cfg.truth_abseta_max <= 0)
    fprintf(out, "NOTE: truth eta cut disabled - the y binning defines the acceptance. D, N, NQ and the\n"
                 "      'out of bins' column include under/overflow, which is INCOMPLETE: the skim keeps\n"
                 "      only particles with |y| or |eta| inside its retention window. Use in-range bins.\n");
  fprintf(out, "\n");
  for (const auto &h : H)
  {
    fprintf(out, "%s: eps(pT) and epsQ(pT), y inside [%g, %g]\n", h.sp.tag, cfg.y_edges.front(), cfg.y_edges.back());
    for (int b = 1; b <= h.den->GetNbinsY(); ++b)
    {
      const double d = h.den->Integral(1, h.den->GetNbinsX(), b, b), n = h.num->Integral(1, h.num->GetNbinsX(), b, b),
                   nq = h.numQ->Integral(1, h.numQ->GetNbinsX(), b, b);
      const std::string eps = d > 0 ? std::string(Form("%.3f", n / d)) : "undefined";
      const std::string epsQ = d > 0 ? std::string(Form("%.3f", nq / d)) : "undefined";
      fprintf(out, "  pT [%5.2f, %5.2f)  D %6.0f  N %6.0f  NQ %6.0f  eps %s  epsQ %s\n",
              h.den->GetYaxis()->GetBinLowEdge(b), h.den->GetYaxis()->GetBinUpEdge(b), d, n, nq, eps.c_str(),
              epsQ.c_str());
    }
  }
  fprintf(out, "\nselected tracks %ld: matched %ld, unmatched %ld; on denominator particles %ld "
               "(of which duplicates %ld = tracks beyond the first per particle)\n",
          C.tracksSelected, C.tracksSelectedMatched, C.tracksSelectedUnmatched, C.tracksOnDenominatorParticles,
          C.duplicateTracks);
  if (cfg.DcaCut())
    fprintf(out, "DCA cut: %ld tracks passing all other cuts failed it (%ld of them with a non-finite DCA)\n",
            C.tracksFailDca, C.tracksDcaNonFinite);
  if (out != stdout) fclose(out);
}

TH1D *StatsHistogram(const Counters &C)
{
  const std::vector<std::pair<const char *, long>> rows = {
      {"files", C.files},
      {"events", C.events},
      {"selected events", C.eventsSelected},
      {"selected true particles", C.truthSelected},
      {"found (any charge)", C.truthFound},
      {"found (correct charge)", C.truthFoundQ},
      {"true particles outside binning", C.truthOutsideBinning},
      {"selected tracks", C.tracksSelected},
      {"selected tracks, matched", C.tracksSelectedMatched},
      {"selected tracks, unmatched", C.tracksSelectedUnmatched},
      {"selected tracks on den. particles", C.tracksOnDenominatorParticles},
      {"duplicate tracks", C.duplicateTracks}};
  auto *h = new TH1D("hTrkEffStats", "TrkEffAnalysis bookkeeping;;count", static_cast<int>(rows.size()), 0.5,
                     rows.size() + 0.5);
  for (size_t i = 0; i < rows.size(); ++i)
  {
    h->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, rows[i].first);
    h->SetBinContent(static_cast<int>(i) + 1, static_cast<double>(rows[i].second));
  }
  return h;
}

const char *kDefinitions =
    "Denominator D(y,pT): G4-primary particles of the species in selected events, pT_true > truth_pt_min, "
    "|eta_true| < truth_abseta_max (if > 0); binned in TRUE y (PDG mass) and TRUE pT. "
    "Trackable denominator (truth_min_tpc_layers = N > 0): only particles with truth clusters in >= N TPC layers; "
    "the efficiency is then CONDITIONAL on reaching the TPC (never combined with the inclusive definition). "
    "Numerator N: particles of D with >= 1 selected track (pT, |eta|, ntpc, nmaps cuts; no DCA) whose dominant truth "
    "particle (most clusters) is this particle; key (job, event, gtrackID); counted once. "
    "NQ: as N with >= 1 such track of correct reconstructed charge. No credit through descendants (pi -> mu). "
    "Event selection: |z_true| < zvtx_max, or (use_reco_zvtx) reconstructed vertex with |z_reco| < zvtx_max. "
    "Histograms are raw counts; efficiencies = TEfficiency (Clopper-Pearson 68.3%); D = 0 bins are undefined.";

// Short human-readable description of the selection, stored as TNamed
// "selection" and printed by the plotting macros on every comparison.
std::string SelectionLabel(const TrkEffConfig &c)
{
  std::string s = c.use_reco_zvtx ? Form("reco vertex, |z_reco| < %g cm", c.zvtx_max) : Form("|z_true| < %g cm", c.zvtx_max);
  s += Form("; TPC >= %d, MVTX >= %d", c.trk_ntpc_min, c.trk_nmvtx_min);
  if (c.trk_nintt_min > 0) s += Form(", INTT >= %d", c.trk_nintt_min);
  if (c.DcaCut()) s += ", " + DcaText(c);
  if (c.TruthEtaExplicit())
    s += c.truth_abseta_max > 0 ? Form("; truth |eta| < %g", c.truth_abseta_max) : "; no truth-eta cut";
  else if (c.V2Cuts() && c.truth_abseta_max <= 0)
    s += "; no truth-eta cut";  // v1/v2 labels unchanged
  s += c.truth_min_tpc_layers > 0 ? Form("; trackable (>= %d TPC layers)", c.truth_min_tpc_layers) : "; inclusive";
  return s;
}

void Usage()
{
  std::cout << "Usage: TrkEffAnalysis -c trkeff.cfg -o out.root [-s key=value ...] [-t summary.txt] [-l skim.list | skim1.root ...]\n"
               "  -c  configuration file (required)\n"
               "  -s  override one configuration key (repeatable), e.g. -s trk_nmvtx_min=3\n"
               "  -o  output ROOT file (required)\n"
               "  -t  text summary (default: printed to stdout)\n"
               "  -l  text file with skim files (one per line); or give the files as arguments\n";
}
}  // namespace

int main(int argc, char **argv)
{
  std::string cfgPath, outPath, summaryPath, listPath;
  std::vector<std::string> files, overrides;
  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "-c" && i + 1 < argc) cfgPath = argv[++i];
    else if (a == "-o" && i + 1 < argc) outPath = argv[++i];
    else if (a == "-t" && i + 1 < argc) summaryPath = argv[++i];
    else if (a == "-l" && i + 1 < argc) listPath = argv[++i];
    else if (a == "-s" && i + 1 < argc) overrides.push_back(argv[++i]);
    else if (a == "-h" || a == "--help")
    {
      Usage();
      return 0;
    }
    else if (!a.empty() && a[0] == '-')
    {
      std::cerr << "FATAL: unknown option " << a << std::endl;
      Usage();
      return 1;
    }
    else files.push_back(a);
  }
  if (cfgPath.empty() || outPath.empty())
  {
    Usage();
    return 1;
  }
  if (!listPath.empty())
    for (const auto &f : ReadFileList(listPath, "")) files.push_back(f);
  if (files.empty())
  {
    std::cerr << "FATAL: no skim files given" << std::endl;
    return 1;
  }

  const TrkEffConfig cfg = ReadTrkEffConfig(cfgPath, overrides);
  TH1::AddDirectory(false);
  std::vector<SpeciesHistograms> H;
  for (const auto &s : cfg.species) H.emplace_back(s, cfg);

  std::cout << "=== " << kVersion << " ===\nconfig: " << cfgPath;
  for (const auto &o : overrides) std::cout << "  [-s " << o << "]";
  std::cout << "\nselection: " << SelectionLabel(cfg) << "\nfiles : " << files.size() << std::endl;
  Counters C;
  int errors = 0;
  std::set<std::string> seen;
  for (const auto &f : files)
  {
    if (!seen.insert(f).second)
    {
      std::cerr << "ERROR: " << f << " given twice - skipped" << std::endl;
      errors++;
      continue;
    }
    if (!ProcessSkimFile(f, cfg, H, C)) errors++;
  }

  TFile out(outPath.c_str(), "RECREATE");
  if (out.IsZombie())
  {
    std::cerr << "FATAL: cannot create " << outPath << std::endl;
    return 1;
  }
  std::string fileList;
  for (const auto &f : files) fileList += f + "\n";
  TNamed("config", cfg.dump.c_str()).Write();
  std::string definitions = kDefinitions;
  if (cfg.V2Cuts())
    definitions += Form(" v2 track cuts in this output: nintt >= %d; %s. The DCA is the straight-line DCA (EvalSkim "
                        "linXYReco/linZReco) to the reconstructed event vertex with the most tracks "
                        "(SvtxVertexMap, final tracking vertex finder).",
                        cfg.trk_nintt_min, DcaText(cfg).c_str());
  if (cfg.TruthEtaExplicit() && cfg.truth_abseta_max > 0)
    definitions += Form(" Truth-eta mode (use_truth_eta_cut = 1): denominator and numerator restricted to generated "
                        "primaries with |eta_true| < %g; the reconstructed |eta| < %g cut is unchanged.",
                        cfg.truth_abseta_max, cfg.trk_abseta_max);
  TNamed("definitions", definitions.c_str()).Write();
  TNamed("selection", SelectionLabel(cfg).c_str()).Write();
  TNamed("program", kVersion).Write();
  TNamed("skim_files", fileList.c_str()).Write();
  StatsHistogram(C)->Write();
  if (cfg.DcaCut())  // separate histogram, so hTrkEffStats keeps the v1 layout
  {
    TH1D hDca("hTrkEffDcaStats", "DCA cut bookkeeping (tracks passing all other cuts);;tracks", 2, 0.5, 2.5);
    hDca.GetXaxis()->SetBinLabel(1, "failed DCA cut");
    hDca.GetXaxis()->SetBinLabel(2, "non-finite DCA");
    hDca.SetBinContent(1, static_cast<double>(C.tracksFailDca));
    hDca.SetBinContent(2, static_cast<double>(C.tracksDcaNonFinite));
    hDca.Write();
  }
  for (auto &h : H)
  {
    h.den->Write();
    h.num->Write();
    h.numQ->Write();
    const char *t = h.sp.tag;
    TH1D *denPt = ProjectPt(h.den, Form("hDen_%s_pt", t)), *numPt = ProjectPt(h.num, Form("hNum_%s_pt", t)),
         *numQPt = ProjectPt(h.numQ, Form("hNumQ_%s_pt", t));
    TH1D *denY = ProjectY(h.den, Form("hDen_%s_y", t)), *numY = ProjectY(h.num, Form("hNum_%s_y", t)),
         *numQY = ProjectY(h.numQ, Form("hNumQ_%s_y", t));
    for (TH1 *p : {static_cast<TH1 *>(denPt), static_cast<TH1 *>(numPt), static_cast<TH1 *>(numQPt),
                   static_cast<TH1 *>(denY), static_cast<TH1 *>(numY), static_cast<TH1 *>(numQY)})
      p->Write();
    const std::string L = h.sp.label;
    MakeEfficiency(*h.num, *h.den, Form("eff_%s_ypt", t), (L + " efficiency (any charge);y_{true};p_{T,true} (GeV/c)").c_str())->Write();
    MakeEfficiency(*h.numQ, *h.den, Form("effQ_%s_ypt", t), (L + " efficiency (correct charge);y_{true};p_{T,true} (GeV/c)").c_str())->Write();
    MakeEfficiency(*numPt, *denPt, Form("eff_%s_pt", t), (L + " efficiency (any charge);p_{T,true} (GeV/c);#varepsilon").c_str())->Write();
    MakeEfficiency(*numQPt, *denPt, Form("effQ_%s_pt", t), (L + " efficiency (correct charge);p_{T,true} (GeV/c);#varepsilon").c_str())->Write();
    MakeEfficiency(*numY, *denY, Form("eff_%s_y", t), (L + " efficiency (any charge);y_{true};#varepsilon").c_str())->Write();
    MakeEfficiency(*numQY, *denY, Form("effQ_%s_y", t), (L + " efficiency (correct charge);y_{true};#varepsilon").c_str())->Write();
  }
  out.Close();
  WriteSummary(summaryPath, cfg, H, C);
  std::cout << "selected events " << C.eventsSelected << " of " << C.events << "; selected true particles "
            << C.truthSelected << ", found " << C.truthFound << " (correct charge " << C.truthFoundQ << ")\n"
            << "wrote " << outPath << (summaryPath.empty() ? "" : " and " + summaryPath) << "\nerrors: " << errors
            << std::endl;
  return errors == 0 ? 0 : 3;
}
