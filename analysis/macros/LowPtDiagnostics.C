// =============================================================================
// LowPtDiagnostics.C - where do true particles get lost at low pT?
//
// For every true particle of the TrkEffAnalysis denominator (same event, truth
// and track selection, read from the "config" stored in the TrkEffAnalysis
// output) this macro records ONE category, using only information stored in
// the EvalSkim trees:
//
//   found          >= 1 selected matched track            (= the efficiency numerator)
//   track, pT cut  matched track(s) exist, none selected; the best one (most
//                  clusters from this particle) fails the reco pT cut first
//   track, eta cut ... fails the reco |eta| cut first
//   track, nTPC    ... fails the TPC-cluster cut first
//   track, nMVTX   ... fails the MVTX-cluster cut first (only if trk_nmvtx_min > 0)
//   no track, >= N TPC layers   no matched track although the particle left truth
//                               clusters in >= trk_ntpc_min TPC layers
//   no track, 1..N-1 TPC layers no matched track; fewer truth TPC layers
//   no track, 0 TPC layers      no matched track; no truth TPC cluster at all
//
// These are bookkeeping facts (which cut or reconstruction stage the particle
// fails), not a physical explanation; the macro deliberately does not assign a
// physics mechanism.
//
// Used by PlotTrkEff.C (page "low-pT diagnostics"); can also be run alone:
//   root -l -b -q 'LowPtDiagnostics.C+("TrkEff_baseline.root","skim_files.list","out.txt")'
// =============================================================================
#ifndef LOWPTDIAGNOSTICS_C
#define LOWPTDIAGNOSTICS_C

#include "TrkEffConfigUtil.h"

#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace lowpt
{
const std::vector<std::string> kCategories = {"found",
                                              "track, fails pT cut",
                                              "track, fails |#eta| cut",
                                              "track, fails n_{TPC} cut",
                                              "track, fails n_{MVTX} cut",
                                              "no track, #geq N_{TPC} truth TPC layers",
                                              "no track, 1..N_{TPC}-1 truth TPC layers",
                                              "no track, 0 truth TPC layers"};
// same categories as plain text (tables)
const std::vector<std::string> kCategoriesText = {"found", "trk:pT", "trk:eta", "trk:nTPC", "trk:nMVTX",
                                                  "none:>=N TPC lay", "none:1..N-1 lay", "none:0 lay"};

// Particle groups (charge-combined, to have meaningful statistics with 1000 events).
struct Group
{
  const char *name, *label;
  int absPdg;
};
const std::vector<Group> kGroups = {{"pions", "#pi^{#pm}", 211}, {"kaons", "K^{#pm}", 321}, {"protons", "p + #bar{p}", 2212}};

// Selections, read from the TrkEffAnalysis output (TNamed "config").
struct Cuts
{
  bool useRecoZ = false;
  double zmax = 8, truthPtMin = 0.1, truthEtaMax = 1.1, trkPtMin = 0.1, trkEtaMax = 1.1;
  int ntpcMin = 18, nmvtxMin = 0;
  std::vector<double> ptEdges;
};

inline std::map<std::string, std::string> ParseConfig(const char *text)
{
  std::map<std::string, std::string> kv;
  std::istringstream in(text);
  std::string l;
  while (std::getline(in, l))
  {
    const auto eq = l.find('=');
    if (eq == std::string::npos) continue;
    auto trim = [](std::string s) {
      s.erase(0, s.find_first_not_of(" \t"));
      s.erase(s.find_last_not_of(" \t") + 1);
      return s;
    };
    kv[trim(l.substr(0, eq))] = trim(l.substr(eq + 1));
  }
  return kv;
}

inline Cuts ReadCuts(TFile &trkEff)
{
  auto *cfg = dynamic_cast<TNamed *>(trkEff.Get("config"));
  if (!cfg) throw std::runtime_error("TrkEffAnalysis output has no 'config'");
  auto kv = ParseConfig(cfg->GetTitle());
  Cuts c;
  c.useRecoZ = std::stoi(kv.at("use_reco_zvtx")) != 0;
  c.zmax = std::stod(kv.at("zvtx_max"));
  c.truthPtMin = std::stod(kv.at("truth_pt_min"));
  c.truthEtaMax = trkeffcfg::EffectiveTruthEtaMax({kv.begin(), kv.end()});
  c.trkPtMin = std::stod(kv.at("trk_pt_min"));
  c.trkEtaMax = std::stod(kv.at("trk_abseta_max"));
  c.ntpcMin = std::stoi(kv.at("trk_ntpc_min"));
  c.nmvtxMin = std::stoi(kv.at("trk_nmvtx_min"));
  std::string e = kv.at("pt_edges");
  for (char &ch : e)
    if (ch == ',') ch = ' ';
  std::istringstream in(e);
  double x;
  while (in >> x) c.ptEdges.push_back(x);
  return c;
}

// Result: counts[group][pT bin][category]
struct Result
{
  Cuts cuts;
  std::vector<std::vector<std::vector<long>>> counts;
  std::vector<TH1D *> recoPtLowestBin;  // per group: reco pT of the best matched track, true pT in the lowest bin
};

struct TrackInfo
{
  float pt, eta;
  short ntpc, nmaps, rank;
};

inline Result Run(const char *trkEffFile, const char *skimList)
{
  TFile f(trkEffFile, "READ");
  Result R;
  R.cuts = ReadCuts(f);
  const Cuts &c = R.cuts;
  const int nPt = static_cast<int>(c.ptEdges.size()) - 1;
  R.counts.assign(kGroups.size(), std::vector<std::vector<long>>(nPt, std::vector<long>(kCategories.size(), 0)));
  for (const auto &g : kGroups)
  {
    auto *h = new TH1D(Form("hLowPtRecoPt_%s", g.name),
                       Form("%s, %.2f #leq p_{T,true} < %.2f GeV/c;p_{T,reco} of the best matched track (GeV/c);particles",
                            g.label, c.ptEdges[0], c.ptEdges[1]),
                       40, 0., 0.4);
    h->SetDirectory(nullptr);
    R.recoPtLowestBin.push_back(h);
  }

  std::ifstream in(skimList);
  std::string path;
  while (std::getline(in, path))
  {
    if (path.empty() || path[0] == '#') continue;
    TFile s(path.c_str(), "READ");
    auto *tE = dynamic_cast<TTree *>(s.Get("EventTree"));
    auto *tT = dynamic_cast<TTree *>(s.Get("TruthTree"));
    auto *tR = dynamic_cast<TTree *>(s.Get("RecoTree"));
    if (!tE || !tT || !tR) throw std::runtime_error("not a skim file: " + path);

    // selected events (same rule as TrkEffAnalysis::PassEventSelection)
    std::map<std::pair<int, int>, bool> selected;
    {
      int job, event;
      float zT, zR;
      short hasVtx;
      tE->SetBranchAddress("job", &job);
      tE->SetBranchAddress("event", &event);
      tE->SetBranchAddress("zTrue", &zT);
      tE->SetBranchAddress("zReco", &zR);
      tE->SetBranchAddress("hasVtx", &hasVtx);
      for (long i = 0; i < tE->GetEntries(); ++i)
      {
        tE->GetEntry(i);
        selected[{job, event}] = c.useRecoZ ? (hasVtx && std::isfinite(zR) && std::fabs(zR) < c.zmax)
                                            : (std::isfinite(zT) && std::fabs(zT) < c.zmax);
      }
    }
    // all matched tracks per (event, gtrackID); job is fixed within a skim file
    std::map<std::pair<int, int>, std::vector<TrackInfo>> tracks;
    {
      int event, gid;
      float pt, eta;
      short ntpc, nmaps, matched, rank;
      tR->SetBranchAddress("event", &event);
      tR->SetBranchAddress("gtrackID", &gid);
      tR->SetBranchAddress("pt", &pt);
      tR->SetBranchAddress("eta", &eta);
      tR->SetBranchAddress("ntpc", &ntpc);
      tR->SetBranchAddress("nmaps", &nmaps);
      tR->SetBranchAddress("matched", &matched);
      tR->SetBranchAddress("rankSameTruth", &rank);
      for (long i = 0; i < tR->GetEntries(); ++i)
      {
        tR->GetEntry(i);
        if (matched) tracks[{event, gid}].push_back({pt, eta, ntpc, nmaps, rank});
      }
    }
    // true particles of the denominator
    {
      int job, event, gid, pdg;
      short isPrimary, gnltpc;
      float pt, eta;
      tT->SetBranchAddress("job", &job);
      tT->SetBranchAddress("event", &event);
      tT->SetBranchAddress("gtrackID", &gid);
      tT->SetBranchAddress("pdg", &pdg);
      tT->SetBranchAddress("isPrimary", &isPrimary);
      tT->SetBranchAddress("gnltpc", &gnltpc);
      tT->SetBranchAddress("pt", &pt);
      tT->SetBranchAddress("eta", &eta);
      for (long i = 0; i < tT->GetEntries(); ++i)
      {
        tT->GetEntry(i);
        if (isPrimary != 1 || !selected[{job, event}]) continue;
        int g = -1;
        for (size_t k = 0; k < kGroups.size(); ++k)
          if (std::abs(pdg) == kGroups[k].absPdg) g = static_cast<int>(k);
        if (g < 0) continue;
        if (!(pt > c.truthPtMin)) continue;
        if (c.truthEtaMax > 0 && !(std::fabs(eta) < c.truthEtaMax)) continue;
        int ib = -1;
        for (int b = 0; b < nPt; ++b)
          if (pt >= c.ptEdges[b] && pt < c.ptEdges[b + 1]) ib = b;
        if (ib < 0) continue;

        int cat;
        auto it = tracks.find({event, gid});
        if (it == tracks.end())
          cat = gnltpc >= c.ntpcMin ? 5 : (gnltpc > 0 ? 6 : 7);
        else
        {
          bool anySelected = false;
          const TrackInfo *best = nullptr;
          for (const auto &t : it->second)
          {
            if (std::isfinite(t.pt) && std::isfinite(t.eta) && t.pt > c.trkPtMin && std::fabs(t.eta) < c.trkEtaMax &&
                t.ntpc >= c.ntpcMin && t.nmaps >= c.nmvtxMin)
              anySelected = true;
            if (t.rank == 1 || !best) best = &t;  // rank 1 = most clusters from this particle
          }
          if (anySelected) cat = 0;
          else if (!(best->pt > c.trkPtMin)) cat = 1;
          else if (!(std::fabs(best->eta) < c.trkEtaMax)) cat = 2;
          else if (best->ntpc < c.ntpcMin) cat = 3;
          else cat = 4;
          if (ib == 0) R.recoPtLowestBin[g]->Fill(best->pt);
        }
        R.counts[g][ib][cat]++;
      }
    }
  }
  return R;
}

// Text table of the categories (fractions per pT bin).
inline void WriteTable(const Result &R, FILE *out)
{
  const auto &e = R.cuts.ptEdges;
  for (size_t g = 0; g < kGroups.size(); ++g)
  {
    fprintf(out, "%s: categories per true-pT bin (counts)\n", kGroups[g].name);
    fprintf(out, "  %-21s", "pT_true (GeV/c)");
    for (const auto &cat : kCategoriesText) fprintf(out, " | %16s", cat.c_str());
    fprintf(out, "\n");
    for (size_t b = 0; b < R.counts[g].size(); ++b)
    {
      long tot = 0;
      for (long n : R.counts[g][b]) tot += n;
      fprintf(out, "  [%4.2f,%4.2f) N=%4ld", e[b], e[b + 1], tot);
      for (long n : R.counts[g][b]) fprintf(out, " | %16ld", n);
      fprintf(out, "\n");
    }
  }
}
}  // namespace lowpt

// Stand-alone entry point.
inline void LowPtDiagnostics(const char *trkEffFile, const char *skimList, const char *outTxt)
{
  const auto R = lowpt::Run(trkEffFile, skimList);
  FILE *out = fopen(outTxt, "w");
  lowpt::WriteTable(R, out);
  fclose(out);
}

#endif
