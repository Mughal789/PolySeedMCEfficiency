// =============================================================================
// ValidateTrkEff.C - independent validation of a TrkEffAnalysis output
//
// Recomputes the denominator and both numerators DIRECTLY from the original
// SvtxEvaluator files (ntp_gpoint, ntp_vertex, ntp_gtrack, ntp_track), with its
// own code (no EvalSkim / skim trees involved), using the cut values stored in
// the TrkEffAnalysis output (TNamed "config"), and compares bin by bin.
//
// Usage:
//   root -l -b -q 'ValidateTrkEff.C+("eval.list","TrkEff.root","report.log","summary.txt","merged.root")'
//     eval.list    evaluator files used for the skims
//     TrkEff.root  TrkEffAnalysis output to validate
//     summary.txt  its text summary ("" = skip the "undefined" check)
//     merged.root  optional: hadd of TrkEffAnalysis runs on split input lists;
//                  must be identical to TrkEff.root ("" = skip)
// =============================================================================
#include "../src/EffTypes.h"
#include "TrkEffConfigUtil.h"

#include <TDatabasePDG.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TNamed.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace vt
{
FILE *gOut = nullptr;
std::map<std::string, std::pair<long, long>> gChk;
std::map<std::string, std::string> gFirst;
void Log(const std::string &s)
{
  printf("%s\n", s.c_str());
  if (gOut) fprintf(gOut, "%s\n", s.c_str());
}
void Test(const std::string &n, bool ok, const std::string &m = "")
{
  auto &c = gChk[n];
  c.first++;
  if (!ok)
  {
    c.second++;
    if (!gFirst.count(n)) gFirst[n] = m;
  }
}
std::map<std::string, std::string> KV(const char *text)
{
  std::map<std::string, std::string> kv;
  std::istringstream in(text);
  std::string l;
  while (std::getline(in, l))
  {
    auto eq = l.find('=');
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
std::vector<double> Nums(const std::string &s)
{
  std::vector<double> v;
  std::string t = s;
  for (char &c : t)
    if (c == ',') c = ' ';
  std::istringstream in(t);
  double x;
  while (in >> x) v.push_back(x);
  return v;
}
int JobOf(const std::string &p)
{
  std::smatch m;
  static const std::regex re("_([0-9]+)\\.root$");
  return std::regex_search(p, m, re) ? std::stoi(m[1]) : -1;
}
}  // namespace vt

void ValidateTrkEff(const char *evalList, const char *trkEffFile, const char *report, const char *summary = "",
                    const char *merged = "")
{
  using namespace vt;
  gOut = fopen(report, "w");
  Log("=== ValidateTrkEff ===");
  Log(std::string("TrkEff file   : ") + trkEffFile);
  Log(std::string("evaluator list: ") + evalList);

  TFile *R = TFile::Open(trkEffFile);
  auto cfg = KV(((TNamed *)R->Get("config"))->GetTitle());
  const bool useReco = std::stoi(cfg["use_reco_zvtx"]) != 0;
  const double zmax = std::stod(cfg["zvtx_max"]), tptmin = std::stod(cfg["truth_pt_min"]),
               tetamax = trkeffcfg::EffectiveTruthEtaMax(cfg), kptmin = std::stod(cfg["trk_pt_min"]),
               ketamax = std::stod(cfg["trk_abseta_max"]);
  const int kntpc = std::stoi(cfg["trk_ntpc_min"]), knmvtx = std::stoi(cfg["trk_nmvtx_min"]);
  // v2 track cuts (absent in outputs written before the keys existed -> off)
  const int knintt = cfg.count("trk_nintt_min") ? std::stoi(cfg["trk_nintt_min"]) : 0;
  const double kdcaxy = cfg.count("trk_dca_xy_max") ? std::stod(cfg["trk_dca_xy_max"]) : 0;
  const double kdcaz = cfg.count("trk_dca_z_max") ? std::stod(cfg["trk_dca_z_max"]) : 0;
  const bool dcaCut = kdcaxy > 0 || kdcaz > 0;
  // trackable denominator (absent in outputs written before the key existed -> inclusive)
  const int tminlayers = cfg.count("truth_min_tpc_layers") ? std::stoi(cfg["truth_min_tpc_layers"]) : 0;
  const std::vector<double> ye = Nums(cfg["y_edges"]), pe = Nums(cfg["pt_edges"]);
  Log(Form("config: use_reco_zvtx=%d zvtx_max=%g truth pt>%g |eta|<%g(%s) TPC layers>=%d ; track pt>%g |eta|<%g ntpc>=%d nmvtx>=%d",
           useReco, zmax, tptmin, tetamax, tetamax > 0 ? "on" : "off", tminlayers, kptmin, ketamax, kntpc, knmvtx));
  if (knintt > 0 || dcaCut)
    Log(Form("v2 track cuts: nintt>=%d, |DCA_xy|<=%g, |DCA_z|<=%g cm (<= 0: off); DCA recomputed here as a straight "
             "line from ntp_track pca along (px,py,pz) to the ntp_vertex vertex with the most tracks",
             knintt, kdcaxy, kdcaz));
  Test("C DCA cuts only together with a reconstructed vertex", !dcaCut || useReco);

  // the identified species of EffTypes.h
  std::vector<int> pdgs, charges;
  std::vector<const char *> tags;
  for (const auto *e : efftypes::MeasuredSpecies(false))
  {
    pdgs.push_back(e->pdg);
    tags.push_back(e->tag);
    charges.push_back(e->charge);
  }
  const int NS = static_cast<int>(pdgs.size());
  std::map<int, int> idx;
  for (int s = 0; s < NS; ++s) idx[pdgs[s]] = s;
  // independent reference histograms (same binning)
  std::vector<TH2D *> rD, rN, rQ;
  for (int s = 0; s < NS; ++s)
  {
    rD.push_back(new TH2D(Form("rD_%s", tags[s]), "", ye.size() - 1, ye.data(), pe.size() - 1, pe.data()));
    rN.push_back(new TH2D(Form("rN_%s", tags[s]), "", ye.size() - 1, ye.data(), pe.size() - 1, pe.data()));
    rQ.push_back(new TH2D(Form("rQ_%s", tags[s]), "", ye.size() - 1, ye.data(), pe.size() - 1, pe.data()));
    TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(pdgs[s]);
    Test("F species charge from TDatabasePDG", p && std::lround(p->Charge() / 3.) == charges[s], tags[s]);
  }

  long nSelEv = 0, nUniqueFound = 0;
  std::ifstream in(evalList);
  std::string f;
  while (std::getline(in, f))
  {
    if (f.empty()) continue;
    const int job = JobOf(f);
    TFile *E = TFile::Open(f.c_str());
    // --- event selection from ntp_gpoint / ntp_vertex --------------------------
    std::map<int, float> zT, zR, xR, yR;
    std::map<int, int> ntrR;
    {
      TTree *t = (TTree *)E->Get("ntp_gpoint");
      float ev, z;
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("event", 1);
      t->SetBranchStatus("gvz", 1);
      t->SetBranchAddress("event", &ev);
      t->SetBranchAddress("gvz", &z);
      for (long i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        zT[(int)std::lround(ev)] = z;
      }
      TTree *v = (TTree *)E->Get("ntp_vertex");
      float vz, vx, vy, nt;
      v->SetBranchStatus("*", 0);
      for (const char *b : {"event", "vx", "vy", "vz", "ntracks"}) v->SetBranchStatus(b, 1);
      v->SetBranchAddress("event", &ev);
      v->SetBranchAddress("vx", &vx);
      v->SetBranchAddress("vy", &vy);
      v->SetBranchAddress("vz", &vz);
      v->SetBranchAddress("ntracks", &nt);
      for (long i = 0; i < v->GetEntries(); ++i)
      {
        v->GetEntry(i);
        if (!(std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz))) continue;
        const int e = std::lround(ev), n = std::isfinite(nt) ? std::lround(nt) : 0;
        if (!zR.count(e) || n > ntrR[e])
        {
          zR[e] = vz;
          xR[e] = vx;
          yR[e] = vy;
          ntrR[e] = n;
        }
      }
    }
    std::set<int> sel;
    for (auto &kv : zT)
    {
      const bool ok = useReco ? (zR.count(kv.first) && std::fabs(zR[kv.first]) < zmax) : std::fabs(kv.second) < zmax;
      if (ok) sel.insert(kv.first);
    }
    nSelEv += sel.size();

    // --- selected tracks -> (event, gtrackID) with charge-correct flag --------
    std::map<std::pair<int, int>, int> pdgOf;  // from ntp_track's own gflavor
    std::map<std::pair<int, int>, bool> found;  // value: any track with correct charge
    {
      TTree *t = (TTree *)E->Get("ntp_track");
      float ev, pt, eta, ntpc, nmaps, gid, q, gfl, nintt = 0, pcax = NAN, pcay = NAN, pcaz = NAN, px = NAN, py = NAN, pz = NAN;
      t->SetBranchStatus("*", 0);
      for (const char *b : {"event", "pt", "eta", "ntpc", "nmaps", "gtrackID", "charge", "gflavor"}) t->SetBranchStatus(b, 1);
      if (knintt > 0)
      {
        t->SetBranchStatus("nintt", 1);
        t->SetBranchAddress("nintt", &nintt);
      }
      if (dcaCut)
      {
        for (const char *b : {"pcax", "pcay", "pcaz", "px", "py", "pz"}) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("pcax", &pcax);
        t->SetBranchAddress("pcay", &pcay);
        t->SetBranchAddress("pcaz", &pcaz);
        t->SetBranchAddress("px", &px);
        t->SetBranchAddress("py", &py);
        t->SetBranchAddress("pz", &pz);
      }
      t->SetBranchAddress("event", &ev);
      t->SetBranchAddress("pt", &pt);
      t->SetBranchAddress("eta", &eta);
      t->SetBranchAddress("ntpc", &ntpc);
      t->SetBranchAddress("nmaps", &nmaps);
      t->SetBranchAddress("gtrackID", &gid);
      t->SetBranchAddress("charge", &q);
      t->SetBranchAddress("gflavor", &gfl);
      for (long i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        const int e = std::lround(ev);
        if (!sel.count(e)) continue;
        const int nm = std::isfinite(nmaps) ? std::lround(nmaps) : 0;
        if (!(pt > kptmin && std::fabs(eta) < ketamax && std::lround(ntpc) >= kntpc && nm >= knmvtx)) continue;
        if (knintt > 0 && !(std::isfinite(nintt) && std::lround(nintt) >= knintt)) continue;
        if (dcaCut)
        {
          // straight line r(s) = pca + s p; s at the transverse closest approach to the vertex V
          const double vx0 = xR.count(e) ? xR[e] : NAN, vy0 = yR.count(e) ? yR[e] : NAN, vz0 = zR.count(e) ? zR[e] : NAN;
          const double pt2 = (double)px * px + (double)py * py;
          const double dx = pcax - vx0, dy = pcay - vy0;
          const double sMin = -(dx * px + dy * py) / pt2;
          const double cx = dx + sMin * px, cy = dy + sMin * py, cz = pcaz + sMin * pz - vz0;
          const double dxy = std::sqrt(cx * cx + cy * cy);
          const bool okxy = kdcaxy <= 0 || (std::isfinite(dxy) && dxy <= kdcaxy);
          const bool okz = kdcaz <= 0 || (std::isfinite(cz) && std::fabs(cz) <= kdcaz);
          if (!(okxy && okz)) continue;
        }
        if (!std::isfinite(gid) || std::lround(gid) == 0) continue;
        TParticlePDG *p = TDatabasePDG::Instance()->GetParticle((int)std::lround(gfl));
        const bool qok = p && std::lround(p->Charge() / 3.) == std::lround(q);
        auto k = std::make_pair(e, (int)std::lround(gid));
        found[k] = found[k] || qok;
      }
    }
    // --- true particles from ntp_gtrack ------------------------------------------
    {
      TTree *t = (TTree *)E->Get("ntp_gtrack");
      float ev, gid, gfl, gpr, gpt, geta, gpx, gpy, gpz, gnltpc;
      t->SetBranchStatus("*", 0);
      for (const char *b : {"event", "gtrackID", "gflavor", "gprimary", "gpt", "geta", "gpx", "gpy", "gpz", "gnltpc"}) t->SetBranchStatus(b, 1);
      t->SetBranchAddress("gnltpc", &gnltpc);
      t->SetBranchAddress("event", &ev);
      t->SetBranchAddress("gtrackID", &gid);
      t->SetBranchAddress("gflavor", &gfl);
      t->SetBranchAddress("gprimary", &gpr);
      t->SetBranchAddress("gpt", &gpt);
      t->SetBranchAddress("geta", &geta);
      t->SetBranchAddress("gpx", &gpx);
      t->SetBranchAddress("gpy", &gpy);
      t->SetBranchAddress("gpz", &gpz);
      for (long i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        if (std::lround(gpr) != 1) continue;
        auto it = idx.find((int)std::lround(gfl));
        if (it == idx.end()) continue;
        const int e = std::lround(ev);
        if (!sel.count(e)) continue;
        if (!(gpt > tptmin)) continue;
        if (tetamax > 0 && !(std::fabs(geta) < tetamax)) continue;
        if (tminlayers > 0 && !(gnltpc >= tminlayers)) continue;  // trackable denominator
        const double m = TDatabasePDG::Instance()->GetParticle(pdgs[it->second])->Mass();
        const double E = std::sqrt(gpx * gpx + gpy * gpy + gpz * gpz + m * m);
        // independent formula in double, rounded to float = the precision at which the skim stores y,
        // so a value on a bin edge (e.g. 0.99999998 -> 1.0f) is binned as in the analysis
        const float y = static_cast<float>(0.5 * std::log((E + gpz) / (E - gpz)));
        const int s = it->second;
        rD[s]->Fill(y, gpt);
        auto fit = found.find({e, (int)std::lround(gid)});
        if (fit == found.end()) continue;
        rN[s]->Fill(y, gpt);
        nUniqueFound++;
        if (fit->second) rQ[s]->Fill(y, gpt);
      }
    }
    E->Close();
    (void)job;
  }

  // ---- compare ---------------------------------------------------------------------
  TH1 *stats = (TH1 *)R->Get("hTrkEffStats");
  auto statOf = [&](const char *label) { return stats->GetBinContent(stats->GetXaxis()->FindBin(label)); };
  Test("E selected events == independent evaluator selection", statOf("selected events") == nSelEv,
       Form("%g vs %ld", statOf("selected events"), nSelEv));
  Test("D no duplicate numerator entries: found particles == unique matched keys",
       statOf("found (any charge)") == nUniqueFound, Form("%g vs %ld", statOf("found (any charge)"), nUniqueFound));
  long nEmpty = 0, nBins = 0;
  double nOverflowMissing = 0;
  for (int s = 0; s < NS; ++s)
  {
    TH2D *hD = (TH2D *)R->Get(Form("hDen_%s_ypt", tags[s])), *hN = (TH2D *)R->Get(Form("hNum_%s_ypt", tags[s])),
         *hQ = (TH2D *)R->Get(Form("hNumQ_%s_ypt", tags[s]));
    for (int bx = 0; bx <= hD->GetNbinsX() + 1; ++bx)
      for (int by = 0; by <= hD->GetNbinsY() + 1; ++by)
      {
        const double d = hD->GetBinContent(bx, by), n = hN->GetBinContent(bx, by), q = hQ->GetBinContent(bx, by);
        const bool inRange = bx >= 1 && bx <= hD->GetNbinsX() && by >= 1 && by <= hD->GetNbinsY();
        // Under/overflow of the denominator is complete only if the truth eta cut lies inside the
        // skim retention window (|eta| < 1.5); without it, particles beyond |y| = 1.5 are not in the skim.
        const bool overflowComplete = tetamax > 0 && tetamax <= 1.5;
        if (inRange)
          Test("T denominator bin-by-bin == evaluator (true y, true pT), in-range bins", d == rD[s]->GetBinContent(bx, by),
               Form("%s bin (%d,%d): %g vs %g", tags[s], bx, by, d, rD[s]->GetBinContent(bx, by)));
        else if (overflowComplete)
          Test("T denominator under/overflow == evaluator (complete: eta cut inside skim window)",
               d == rD[s]->GetBinContent(bx, by), Form("%s bin (%d,%d): %g vs %g", tags[s], bx, by, d, rD[s]->GetBinContent(bx, by)));
        else
        {
          Test("T denominator under/overflow <= evaluator (eta cut off: skim keeps |y|<1.5 or |eta|<1.5 only)",
               d <= rD[s]->GetBinContent(bx, by), Form("%s bin (%d,%d): %g vs %g", tags[s], bx, by, d, rD[s]->GetBinContent(bx, by)));
          nOverflowMissing += rD[s]->GetBinContent(bx, by) - d;
        }
        Test("N numerator (any charge) bin-by-bin == evaluator", n == rN[s]->GetBinContent(bx, by),
             Form("%s bin (%d,%d): %g vs %g", tags[s], bx, by, n, rN[s]->GetBinContent(bx, by)));
        Test("Q numerator (correct charge) bin-by-bin == evaluator", q == rQ[s]->GetBinContent(bx, by),
             Form("%s bin (%d,%d): %g vs %g", tags[s], bx, by, q, rQ[s]->GetBinContent(bx, by)));
        Test("B NumQ <= Num <= Den in every bin", q <= n && n <= d, Form("%s (%d,%d)", tags[s], bx, by));
        Test("B raw counts: error^2 == content (unweighted)", std::fabs(hD->GetBinError(bx, by) * hD->GetBinError(bx, by) - d) < 1e-9);
        if (bx >= 1 && bx <= hD->GetNbinsX() && by >= 1 && by <= hD->GetNbinsY())
        {
          nBins++;
          if (d == 0) nEmpty++;
        }
      }
    // projections and TEfficiency consistency
    TH1D *dPt = (TH1D *)R->Get(Form("hDen_%s_pt", tags[s])), *nPt = (TH1D *)R->Get(Form("hNum_%s_pt", tags[s]));
    TH1D *dY = (TH1D *)R->Get(Form("hDen_%s_y", tags[s]));
    for (int b = 1; b <= dPt->GetNbinsX(); ++b)
    {
      Test("P pT projection == sum of in-range y bins", dPt->GetBinContent(b) == hD->Integral(1, hD->GetNbinsX(), b, b) &&
                                                             nPt->GetBinContent(b) == hN->Integral(1, hN->GetNbinsX(), b, b));
    }
    for (int b = 1; b <= dY->GetNbinsX(); ++b)
      Test("P y projection == sum of in-range pT bins", dY->GetBinContent(b) == hD->Integral(b, b, 1, hD->GetNbinsY()));
    for (const char *k : {"ypt", "pt", "y"})
    {
      TEfficiency *e = (TEfficiency *)R->Get(Form("eff_%s_%s", tags[s], k));
      TEfficiency *eq = (TEfficiency *)R->Get(Form("effQ_%s_%s", tags[s], k));
      TH1 *den = (TH1 *)R->Get(Form("hDen_%s_%s", tags[s], k)), *num = (TH1 *)R->Get(Form("hNum_%s_%s", tags[s], k)),
          *numq = (TH1 *)R->Get(Form("hNumQ_%s_%s", tags[s], k));
      bool ok = e && eq;
      for (int g = 0; ok && g < den->GetNcells(); ++g)
        ok = e->GetTotalHistogram()->GetBinContent(g) == den->GetBinContent(g) &&
             e->GetPassedHistogram()->GetBinContent(g) == num->GetBinContent(g) &&
             eq->GetPassedHistogram()->GetBinContent(g) == numq->GetBinContent(g);
      Test("P TEfficiency passed/total == numerator/denominator histograms", ok, Form("%s %s", tags[s], k));
      if (e) Test("P TEfficiency uses Clopper-Pearson 68.3%", e->GetStatisticOption() == TEfficiency::kFCP && std::fabs(e->GetConfidenceLevel() - 0.683) < 1e-9);
    }
  }

  // ---- undefined bins in the text summary --------------------------------------
  if (summary && *summary)
  {
    std::ifstream sf(summary);
    std::string line, sp;
    long undefLines = 0, expected = 0;
    while (std::getline(sf, line))
      if (line.find("eps undefined") != std::string::npos) undefLines++;
    for (int s = 0; s < NS; ++s)
    {
      TH1D *dPt = (TH1D *)R->Get(Form("hDen_%s_pt", tags[s]));
      for (int b = 1; b <= dPt->GetNbinsX(); ++b) expected += dPt->GetBinContent(b) == 0;
    }
    Test("U empty pT bins are reported as 'undefined' in the summary", undefLines == expected,
         Form("%ld lines vs %ld empty bins", undefLines, expected));
  }

  // ---- split vs full ----------------------------------------------------------------
  if (merged && *merged)
  {
    TFile *M = TFile::Open(merged);
    TIter next(R->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
      TObject *a = key->ReadObj();
      TObject *b = M->Get(key->GetName());
      if (a->InheritsFrom(TH1::Class()))
      {
        bool ok = b != nullptr;
        TH1 *ha = (TH1 *)a, *hb = (TH1 *)b;
        for (int g = 0; ok && g < ha->GetNcells(); ++g)
          ok = ha->GetBinContent(g) == hb->GetBinContent(g) && std::fabs(ha->GetBinError(g) - hb->GetBinError(g)) < 1e-9;
        // hTrkEffStats: "files" differs only in the order of summation, contents must still agree
        Test("S split+hadd == full run (all histograms, every bin incl. errors)", ok, key->GetName());
      }
      else if (a->InheritsFrom(TEfficiency::Class()))
      {
        TEfficiency *ea = (TEfficiency *)a, *eb = (TEfficiency *)b;
        bool ok = eb != nullptr;
        for (int g = 0; ok && g < ea->GetTotalHistogram()->GetNcells(); ++g)
          ok = ea->GetTotalHistogram()->GetBinContent(g) == eb->GetTotalHistogram()->GetBinContent(g) &&
               ea->GetPassedHistogram()->GetBinContent(g) == eb->GetPassedHistogram()->GetBinContent(g);
        Test("S split+hadd == full run (TEfficiency passed/total)", ok, key->GetName());
      }
    }
  }

  Log("");
  Log("--- checks (tests / failures) ---");
  long nFail = 0;
  for (auto &kv : gChk)
  {
    const bool ok = kv.second.second == 0;
    if (!ok) nFail++;
    Log(Form("%s  %-84s %7ld tests %5ld failed%s", ok ? "PASS" : "FAIL", kv.first.c_str(), kv.second.first, kv.second.second,
             ok ? "" : Form("  (first: %s)", gFirst[kv.first].c_str())));
  }
  Log("");
  Log(Form("information: selected events %ld; in-range (y,pT) bins %ld, of which empty denominator (undefined efficiency) %ld",
           nSelEv, nBins, nEmpty));
  if (nOverflowMissing > 0)
    Log(Form("  truth eta cut off: %.0f true particles in under/overflow (y beyond the binning) are not in the skim (by design)",
             nOverflowMissing));
  for (int s = 0; s < NS; ++s)
    Log(Form("  %-5s D %6.0f  N %6.0f  NQ %6.0f   (evaluator recomputation)", tags[s], rD[s]->GetEntries(),
             rN[s]->GetEntries(), rQ[s]->GetEntries()));
  Log("");
  Log(Form("SUMMARY: %zu checks, %ld failed", gChk.size(), nFail));
  Log(nFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  fclose(gOut);
  gOut = nullptr;
}
