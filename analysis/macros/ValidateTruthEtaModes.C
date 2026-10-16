// =============================================================================
// ValidateTruthEtaModes.C - consistency of the two truth-eta definitions
// (TrkEffAnalysis v3, use_truth_eta_cut = 0 / 1) on the same sample
//
// Inputs: the two TrkEffAnalysis outputs (same configuration apart from the
// truth-eta switch) and the skim files they were made from.
//  C1 configurations identical apart from the truth-eta keys; the reconstructed
//     |eta| cut (trk_abseta_max) and every track / event cut are the same
//  C2 the effective truth-eta threshold is off in the inclusive output and > 0
//     in the restricted output
//  B1 NQ <= N <= D in every bin of both outputs
//  B2 restricted <= inclusive in every bin, for D, N and NQ
//  P1 particle level, recomputed here from the skim (own event, truth and track
//     selection, cut values read from the config record):
//       D_incl - D_restr == generated primaries with |eta_true| >= threshold
//       N_incl - N_restr == found particles with |eta_true| >= threshold
//       NQ_incl - NQ_restr == found, correct-charge particles with |eta_true| >= threshold
//     in every (species, y, pT) cell incl. under/overflow - i.e. the cut removes
//     exactly the particles outside the truth acceptance from denominator AND
//     numerator, and no accepted numerator particle lies outside it
//  P2 the independent recomputation also reproduces D, N, NQ of both outputs
//
// Usage: root -l -b -q 'ValidateTruthEtaModes.C+("incl.root","restr.root","skim_files.list","log.txt")'
// =============================================================================
#include "TrkEffConfigUtil.h"

#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>

namespace valeta
{
FILE *gLog = nullptr;
int gN = 0, gFail = 0;
void Out(const std::string &s)
{
  printf("%s\n", s.c_str());
  if (gLog) fprintf(gLog, "%s\n", s.c_str());
}
void Check(bool ok, const std::string &what)
{
  ++gN;
  if (!ok) ++gFail;
  Out(std::string(ok ? "PASS  " : "FAIL  ") + what);
}
std::map<std::string, std::string> Keys(const char *text)
{
  std::map<std::string, std::string> m;
  std::istringstream in(text);
  std::string line;
  auto trim = [](std::string x) {
    x.erase(0, x.find_first_not_of(" \t"));
    x.erase(x.find_last_not_of(" \t") + 1);
    return x;
  };
  while (std::getline(in, line))
  {
    const auto eq = line.find('=');
    if (eq != std::string::npos) m[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
  }
  return m;
}
const char *kSp[6] = {"pip", "pim", "kp", "km", "p", "pbar"};
const int kPdg[6] = {211, -211, 321, -321, 2212, -2212};
}  // namespace valeta

void ValidateTruthEtaModes(const char *inclPath, const char *restrPath, const char *skimList, const char *logPath)
{
  using namespace valeta;
  gLog = fopen(logPath, "w");
  Out("=== ValidateTruthEtaModes ===");
  Out(std::string("inclusive (switch off): ") + inclPath);
  Out(std::string("restricted (switch on): ") + restrPath);
  TFile I(inclPath, "READ"), Rf(restrPath, "READ");
  const auto ci = Keys(I.Get<TNamed>("config")->GetTitle()), cr = Keys(Rf.Get<TNamed>("config")->GetTitle());
  Out(std::string("  inclusive : ") + I.Get<TNamed>("selection")->GetTitle());
  Out(std::string("  restricted: ") + Rf.Get<TNamed>("selection")->GetTitle());

  // ---- C: configurations ---------------------------------------------------------------
  const std::set<std::string> etaKeys = {"truth_abseta_max", "use_truth_eta_cut", "truth_eta_cut_max"};
  std::string diff;
  for (const auto &[k, v] : ci)
    if (!etaKeys.count(k) && (!cr.count(k) || cr.at(k) != v)) diff += k + " ";
  Check(diff.empty(), "C1 configurations identical apart from the truth-eta keys (incl. trk_abseta_max = " +
                          ci.at("trk_abseta_max") + ")" + (diff.empty() ? "" : "; differ: " + diff));
  const double etaI = trkeffcfg::EffectiveTruthEtaMax(ci), etaR = trkeffcfg::EffectiveTruthEtaMax(cr);
  Check(etaI <= 0 && etaR > 0, Form("C2 effective truth-eta threshold: inclusive %g (off), restricted %g", etaI, etaR));

  // ---- B: bin relations -----------------------------------------------------------------------
  long badOrder = 0, badSubset = 0, cells = 0;
  for (const char *sp : kSp)
    for (const char *k : {"Den", "Num", "NumQ"})
    {
      const TH2D *a = I.Get<TH2D>(Form("h%s_%s_ypt", k, sp)), *b = Rf.Get<TH2D>(Form("h%s_%s_ypt", k, sp));
      for (int c = 0; c < a->GetNcells(); ++c)
      {
        ++cells;
        badSubset += b->GetBinContent(c) > a->GetBinContent(c);
      }
    }
  for (TFile *F : {&I, &Rf})
    for (const char *sp : kSp)
    {
      const TH2D *d = F->Get<TH2D>(Form("hDen_%s_ypt", sp)), *n = F->Get<TH2D>(Form("hNum_%s_ypt", sp)),
                 *q = F->Get<TH2D>(Form("hNumQ_%s_ypt", sp));
      for (int c = 0; c < d->GetNcells(); ++c)
        badOrder += !(q->GetBinContent(c) <= n->GetBinContent(c) && n->GetBinContent(c) <= d->GetBinContent(c));
    }
  Check(badOrder == 0, Form("B1 NQ <= N <= D in every cell of both outputs (%ld violations)", badOrder));
  Check(badSubset == 0, Form("B2 restricted <= inclusive for D, N, NQ in every cell (%ld cells, %ld violations)", cells, badSubset));

  // ---- P: particle-level recomputation from the skim --------------------------------------------
  // cuts from the config record (same keys TrkEffAnalysis uses)
  const bool useReco = std::stoi(ci.at("use_reco_zvtx")) != 0;
  const double zmax = std::stod(ci.at("zvtx_max")), tptmin = std::stod(ci.at("truth_pt_min"));
  const double kpt = std::stod(ci.at("trk_pt_min")), keta = std::stod(ci.at("trk_abseta_max"));
  const int kntpc = std::stoi(ci.at("trk_ntpc_min")), knmvtx = std::stoi(ci.at("trk_nmvtx_min"));
  const int knintt = ci.count("trk_nintt_min") ? std::stoi(ci.at("trk_nintt_min")) : 0;
  const double kdxy = ci.count("trk_dca_xy_max") ? std::stod(ci.at("trk_dca_xy_max")) : 0;
  const double kdz = ci.count("trk_dca_z_max") ? std::stod(ci.at("trk_dca_z_max")) : 0;
  TChain ev("EventTree"), rt("RecoTree"), tt("TruthTree");
  std::ifstream in(skimList);
  std::string f;
  while (std::getline(in, f))
    if (!f.empty())
    {
      ev.Add(f.c_str());
      rt.Add(f.c_str());
      tt.Add(f.c_str());
    }
  using Key = std::tuple<int, int, int>;
  std::set<std::pair<int, int>> selEv;
  {
    int job, event;
    float zT, zR;
    short hv;
    ev.SetBranchAddress("job", &job);
    ev.SetBranchAddress("event", &event);
    ev.SetBranchAddress("zTrue", &zT);
    ev.SetBranchAddress("zReco", &zR);
    ev.SetBranchAddress("hasVtx", &hv);
    for (Long64_t i = 0; i < ev.GetEntries(); ++i)
    {
      ev.GetEntry(i);
      if (useReco ? (hv && std::isfinite(zR) && std::fabs(zR) < zmax) : (std::isfinite(zT) && std::fabs(zT) < zmax))
        selEv.insert({job, event});
    }
  }
  std::map<Key, bool> found;  // value: any selected track with correct charge
  {
    int job, event, g;
    float pt, eta, lxy = NAN, lz = NAN;
    short ntpc, nmaps, nintt, m, cq;
    rt.SetBranchAddress("job", &job);
    rt.SetBranchAddress("event", &event);
    rt.SetBranchAddress("gtrackID", &g);
    rt.SetBranchAddress("pt", &pt);
    rt.SetBranchAddress("eta", &eta);
    rt.SetBranchAddress("ntpc", &ntpc);
    rt.SetBranchAddress("nmaps", &nmaps);
    rt.SetBranchAddress("nintt", &nintt);
    rt.SetBranchAddress("matched", &m);
    rt.SetBranchAddress("chargeOK", &cq);
    rt.SetBranchAddress("linXYReco", &lxy);
    rt.SetBranchAddress("linZReco", &lz);
    for (Long64_t i = 0; i < rt.GetEntries(); ++i)
    {
      rt.GetEntry(i);
      if (!m || !selEv.count({job, event})) continue;
      if (!(pt > kpt && std::fabs(eta) < keta && ntpc >= kntpc && nmaps >= knmvtx && nintt >= knintt)) continue;
      if (kdxy > 0 && !(std::isfinite(lxy) && std::fabs(lxy) <= kdxy)) continue;
      if (kdz > 0 && !(std::isfinite(lz) && std::fabs(lz) <= kdz)) continue;
      const Key k{job, event, g};
      found[k] = found[k] || cq == 1;
    }
  }
  // per-species histograms (same binning as the outputs) of the recomputation
  std::map<int, int> sIdx;
  for (int s = 0; s < 6; ++s) sIdx[kPdg[s]] = s;
  std::vector<TH2D *> D[2], N[2], Q[2], Dout, Nout, Qout;  // [0] inclusive, [1] restricted; *out = |eta| >= thr
  for (int s = 0; s < 6; ++s)
  {
    const TH2D *ref = I.Get<TH2D>(Form("hDen_%s_ypt", kSp[s]));
    auto mk = [&](const char *n) {
      auto *h = static_cast<TH2D *>(ref->Clone(Form("%s_%s", n, kSp[s])));
      h->Reset();
      h->SetDirectory(nullptr);
      return h;
    };
    for (int m = 0; m < 2; ++m)
    {
      D[m].push_back(mk(Form("D%d", m)));
      N[m].push_back(mk(Form("N%d", m)));
      Q[m].push_back(mk(Form("Q%d", m)));
    }
    Dout.push_back(mk("Dout"));
    Nout.push_back(mk("Nout"));
    Qout.push_back(mk("Qout"));
  }
  long nOutFound = 0;
  {
    int job, event, g, pdg;
    short prim;
    float pt, eta, y;
    tt.SetBranchAddress("job", &job);
    tt.SetBranchAddress("event", &event);
    tt.SetBranchAddress("gtrackID", &g);
    tt.SetBranchAddress("pdg", &pdg);
    tt.SetBranchAddress("isPrimary", &prim);
    tt.SetBranchAddress("pt", &pt);
    tt.SetBranchAddress("eta", &eta);
    tt.SetBranchAddress("y", &y);
    for (Long64_t i = 0; i < tt.GetEntries(); ++i)
    {
      tt.GetEntry(i);
      if (prim != 1 || !sIdx.count(pdg) || !selEv.count({job, event}) || !(pt > tptmin)) continue;
      const int s = sIdx[pdg];
      const bool inAcc = std::isfinite(eta) && std::fabs(eta) < etaR;
      auto it = found.find(Key{job, event, g});
      const bool isFound = it != found.end(), isQ = isFound && it->second;
      for (int m = 0; m < 2; ++m)
      {
        if (m == 1 && !inAcc) continue;
        D[m][s]->Fill(y, pt);
        if (isFound) N[m][s]->Fill(y, pt);
        if (isQ) Q[m][s]->Fill(y, pt);
      }
      if (!inAcc)
      {
        Dout[s]->Fill(y, pt);
        if (isFound)
        {
          Nout[s]->Fill(y, pt);
          ++nOutFound;
        }
        if (isQ) Qout[s]->Fill(y, pt);
      }
    }
  }
  long badRemoved = 0, badRepro = 0, removedD = 0, removedN = 0;
  for (int s = 0; s < 6; ++s)
    for (int c = 0; c < Dout[s]->GetNcells(); ++c)
    {
      const auto cell = [&](TFile &F, const char *k) { return F.Get<TH2D>(Form("h%s_%s_ypt", k, kSp[s]))->GetBinContent(c); };
      const double dD = cell(I, "Den") - cell(Rf, "Den"), dN = cell(I, "Num") - cell(Rf, "Num"), dQ = cell(I, "NumQ") - cell(Rf, "NumQ");
      removedD += std::lround(dD);
      removedN += std::lround(dN);
      badRemoved += dD != Dout[s]->GetBinContent(c) || dN != Nout[s]->GetBinContent(c) || dQ != Qout[s]->GetBinContent(c);
      badRepro += D[0][s]->GetBinContent(c) != cell(I, "Den") || N[0][s]->GetBinContent(c) != cell(I, "Num") ||
                  Q[0][s]->GetBinContent(c) != cell(I, "NumQ") || D[1][s]->GetBinContent(c) != cell(Rf, "Den") ||
                  N[1][s]->GetBinContent(c) != cell(Rf, "Num") || Q[1][s]->GetBinContent(c) != cell(Rf, "NumQ");
    }
  Check(badRemoved == 0, Form("P1 removed by the cut (incl - restr) == particles with |eta_true| >= %g, cell by cell, for D, N and "
                              "NQ (removed: D %ld, N %ld; %ld cells differ)",
                              etaR, removedD, removedN, badRemoved));
  Check(badRepro == 0, Form("P2 independent recomputation from the skim == both outputs (D, N, NQ, all cells; %ld differ)", badRepro));
  Out(Form("  info: %ld found (numerator) particles of the inclusive output have |eta_true| >= %g; none of them is in "
           "the restricted numerator (P1)",
           nOutFound, etaR));
  Out(Form("%d checks, %d failed", gN, gFail));
  Out(gFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  if (gLog) fclose(gLog);
}
