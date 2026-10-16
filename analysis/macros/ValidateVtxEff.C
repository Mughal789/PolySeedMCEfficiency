// =============================================================================
// ValidateVtxEff.C
//
// Validation of the VtxEffAnalysis output on the two audit files
// (outputeval_002642.root and outputeval_000000.root).
//
//  A. Reproduce the counts established in the earlier audit
//     (summary_evaluator.txt / summary_MBonly_polyseed.txt).
//  B. Event-key uniqueness: (job, event) unique; events 0..999 per job.
//  C. Internal consistency of histograms (totals, response-matrix sums,
//     numerator <= denominator, tree vs histograms).
//  D. No NaN/Inf in any histogram; unique object names in the output file.
//
// Usage (from VertexEfficiencyStudy/):
//   root -l -b -q 'macros/ValidateVtxEff.C+("output/VtxEff_validation.root","logs/validation.log")'
// Prints PASS/FAIL per check; the final line says whether validation passed.
// =============================================================================

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TKey.h>
#include <TNamed.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::ofstream gLog;
int gFail = 0, gPass = 0;

void Out(const std::string &s)
{
  std::cout << s << std::endl;
  gLog << s << std::endl;
}

void Check(bool ok, const std::string &what)
{
  (ok ? gPass : gFail)++;
  Out(std::string(ok ? "  PASS  " : "  FAIL  ") + what);
}

template <class T>
std::string Cmp(const char *what, T got, T exp)
{
  std::ostringstream o;
  o << what << ": got " << got << ", expected " << exp;
  return o.str();
}

// Expected numbers from the audit (per job).
struct Expected
{
  long events, vertexed, vertexless, tracks, good;
  long vtxlessWithTrk, vtxlessWithGood, zeroTrk;  // -1 = not audited
  long goodInVtxless, trkInVtxless;               // -1 = not audited
  std::vector<std::pair<long, long>> effN;        // n = 0..5 : (vertexed, all)
  std::pair<long, long> effGe6;                   // n >= 6
};
}  // namespace

void ValidateVtxEff(const char *inFile = "output/VtxEff_validation.root", const char *logFile = "logs/validation.log")
{
  gLog.open(logFile);
  TFile *f = TFile::Open(inFile, "READ");
  if (!f || f->IsZombie())
  {
    Out("FATAL: cannot open " + std::string(inFile));
    return;
  }
  Out("Validation of " + std::string(inFile));
  Out("Reference numbers: summary_evaluator.txt and summary_MBonly_polyseed.txt (earlier audit)");

  std::map<int, Expected> expected;
  expected[2642] = {1000, 211, 789, 10710, 8442, 665, 642, 124, 5575, 7173,
                    {{0, 147}, {1, 46}, {2, 50}, {5, 57}, {8, 57}, {14, 69}}, {181, 574}};
  expected[0] = {1000, 218, 782, 10545, 8336, 644, 626, 138, -1, -1,
                 {{0, 156}, {0, 41}, {3, 45}, {5, 56}, {13, 67}, {8, 64}}, {189, 571}};

  // ---------------------------------------------------------------------------
  // A + B: per-job counts from the per-collision tree
  // ---------------------------------------------------------------------------
  TTree *t = dynamic_cast<TTree *>(f->Get("evt"));
  if (!t)
  {
    Out("FATAL: per-event tree 'evt' missing (write_event_tree = 0?)");
    return;
  }
  int job, event, n, nAll, hasVtx, N, nGoodPrim, nGoodSec, nGoodUnm, nGoodDup;
  t->SetBranchAddress("job", &job);
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("nAll", &nAll);
  t->SetBranchAddress("hasVtx", &hasVtx);
  t->SetBranchAddress("N", &N);
  t->SetBranchAddress("nGoodPrim", &nGoodPrim);
  t->SetBranchAddress("nGoodSec", &nGoodSec);
  t->SetBranchAddress("nGoodUnm", &nGoodUnm);
  t->SetBranchAddress("nGoodDup", &nGoodDup);

  struct Got
  {
    long events = 0, vertexed = 0, vertexless = 0, tracks = 0, good = 0;
    long vtxlessWithTrk = 0, vtxlessWithGood = 0, zeroTrk = 0, goodInVtxless = 0, trkInVtxless = 0;
    long eff[6][2] = {{0}};
    long ge6[2] = {0, 0};
    double sumN = 0, sumn = 0;
    long N0 = 0, prim = 0, sec = 0, unm = 0, dup = 0;
    std::set<int> events_seen;
    long dupKeys = 0;
  };
  std::map<int, Got> got;
  std::set<long long> keys;
  long dupGlobal = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    Got &g = got[job];
    if (!keys.insert(1000000LL * job + event).second) dupGlobal++;
    if (!g.events_seen.insert(event).second) g.dupKeys++;
    g.events++;
    g.vertexed += hasVtx;
    g.vertexless += !hasVtx;
    g.tracks += nAll;
    g.good += n;
    g.sumN += N;
    g.sumn += n;
    g.N0 += (N == 0);
    g.prim += nGoodPrim;
    g.sec += nGoodSec;
    g.unm += nGoodUnm;
    g.dup += nGoodDup;
    if (nAll == 0) g.zeroTrk++;
    if (!hasVtx)
    {
      if (nAll > 0) g.vtxlessWithTrk++;
      if (n > 0) g.vtxlessWithGood++;
      g.goodInVtxless += n;
      g.trkInVtxless += nAll;
    }
    if (n <= 5)
    {
      g.eff[n][0] += hasVtx;
      g.eff[n][1]++;
    }
    else
    {
      g.ge6[0] += hasVtx;
      g.ge6[1]++;
    }
  }

  // The audit counted good tracks without an MVTX requirement. With
  // good_nmvtx_min > 0 (config string contains " nmvtx>=K") only the counts that
  // do not depend on n can be compared; the n-dependent ones are skipped.
  int nmvtxMin = 0;
  if (auto *cfgObj = dynamic_cast<TNamed *>(f->Get("config")))
  {
    const std::string c = cfgObj->GetTitle();
    const auto p = c.find("nmvtx>=");
    if (p != std::string::npos) nmvtxMin = std::stoi(c.substr(p + 7));
  }
  const bool nAudit = (nmvtxMin == 0);

  Out("\n[A] Reproduction of audited counts");
  if (!nAudit)
    Out(Form("  NOTE  good tracks require nmvtx >= %d: audit checks that depend on n are skipped "
             "(events, vertices and all-track counts are still compared)", nmvtxMin));
  for (auto &[j, e] : expected)
  {
    if (!got.count(j))
    {
      // In a larger production the audit jobs may not be included.
      Out(Form("  SKIP  job %06d not in this output - audit comparison not possible", j));
      continue;
    }
    Got &g = got[j];
    if (g.events < e.events)
    {
      // Test run with an event limit (-n): the audit numbers are for full files.
      Out(Form("  SKIP  job %06d: only %ld of %ld events processed (event limit) - audit comparison not applicable", j,
               g.events, e.events));
      continue;
    }
    Out(Form(" job %06d", j));
    Check(g.events == e.events, Cmp("events", g.events, e.events));
    Check(g.vertexed == e.vertexed, Cmp("events with tracking vertex", g.vertexed, e.vertexed));
    Check(g.vertexless == e.vertexless, Cmp("vertexless events", g.vertexless, e.vertexless));
    Check(g.tracks == e.tracks, Cmp("reconstructed tracks (all ntp_track rows)", g.tracks, e.tracks));
    if (nAudit) Check(g.good == e.good, Cmp("good tracks", g.good, e.good));
    if (e.vtxlessWithTrk >= 0)
      Check(g.vtxlessWithTrk == e.vtxlessWithTrk, Cmp("vertexless events with >=1 track", g.vtxlessWithTrk, e.vtxlessWithTrk));
    if (nAudit && e.vtxlessWithGood >= 0)
      Check(g.vtxlessWithGood == e.vtxlessWithGood,
            Cmp("vertexless events with >=1 good track", g.vtxlessWithGood, e.vtxlessWithGood));
    if (e.zeroTrk >= 0) Check(g.zeroTrk == e.zeroTrk, Cmp("events with zero tracks", g.zeroTrk, e.zeroTrk));
    if (nAudit && e.goodInVtxless >= 0)
      Check(g.goodInVtxless == e.goodInVtxless, Cmp("good tracks in vertexless events", g.goodInVtxless, e.goodInVtxless));
    if (e.trkInVtxless >= 0)
      Check(g.trkInVtxless == e.trkInVtxless, Cmp("tracks in vertexless events (not filtered by vertexID)", g.trkInVtxless, e.trkInVtxless));
    for (int k = 0; nAudit && k <= 5; ++k)
      Check(g.eff[k][0] == e.effN[k].first && g.eff[k][1] == e.effN[k].second,
            Form("eps_vtx(n=%d): got %ld/%ld, expected %ld/%ld", k, g.eff[k][0], g.eff[k][1], e.effN[k].first,
                 e.effN[k].second));
    if (nAudit)
      Check(g.ge6[0] == e.effGe6.first && g.ge6[1] == e.effGe6.second,
          Form("eps_vtx(n>=6): got %ld/%ld, expected %ld/%ld", g.ge6[0], g.ge6[1], e.effGe6.first, e.effGe6.second));
    Out(Form("   info: <n> = %.3f, <N> = %.3f, N=0 fraction = %.3f; good tracks: primary-matched %ld, "
             "secondary-matched %ld, unmatched %ld, duplicates %ld",
             g.sumn / g.events, g.sumN / g.events, double(g.N0) / g.events, g.prim, g.sec, g.unm, g.dup));
    if (j == 2642 && nAudit)
    {
      // From summary_evaluator.txt: <n> = 8.44, <N> = 4.60, 18% N=0,
      // 4,075 primary-matched, 4,367 secondary-matched, 2,644 duplicates.
      Check(std::fabs(g.sumn / g.events - 8.44) < 0.006, "<n> = 8.44 (audit)");
      Check(std::fabs(g.sumN / g.events - 4.60) < 0.006, "<N> = 4.60 (audit)");
      Check(g.prim == 4075 && g.sec == 4367, Form("primary/secondary-matched good tracks %ld/%ld vs audit 4075/4367", g.prim, g.sec));
      Check(g.dup == 2644, Cmp("duplicate good tracks", g.dup, 2644L));
    }
  }

  Out("\n[B] Event-key uniqueness");
  Check(dupGlobal == 0, Form("(job,event) keys unique: %ld duplicates", dupGlobal));
  for (auto &[j, g] : got)
  {
    const bool contiguous = !g.events_seen.empty() && *g.events_seen.begin() == 0 &&
                            *g.events_seen.rbegin() == static_cast<int>(g.events_seen.size()) - 1;
    Check(g.dupKeys == 0 && contiguous, Form("job %06d: events 0..%zu each exactly once", j, g.events_seen.size() - 1));
  }

  // ---------------------------------------------------------------------------
  // C: histogram internal consistency
  // ---------------------------------------------------------------------------
  Out("\n[C] Histogram consistency");
  auto H1 = [&](const char *nm) { return dynamic_cast<TH1D *>(f->Get(nm)); };
  auto H2 = [&](const char *nm) { return dynamic_cast<TH2D *>(f->Get(nm)); };
  TH1D *hStats = H1("hStats");
  auto stat = [&](const char *lab) { return hStats->GetBinContent(hStats->GetXaxis()->FindBin(lab)); };
  const double nEv = stat("events"), nV = stat("vertexed"), nNV = stat("vertexless");
  const long treeEv = t->GetEntries();
  long treeV = 0;
  for (auto &[j, g] : got) treeV += g.vertexed;
  Check(nEv == treeEv, Form("hStats events (%.0f) == tree rows (%ld)", nEv, treeEv));
  Check(nV + nNV == nEv, Form("vertexed + vertexless (%.0f + %.0f) == events", nV, nNV));
  Check(nV == treeV, Form("hStats vertexed (%.0f) == tree vertexed (%ld)", nV, treeV));
  const char *allH[] = {"hn_all", "hN_all", "hnW_all", "hNW_all", "hnAll_all", "hZtrue_all"};
  for (const char *nm : allH)
  {
    TH1D *h = H1(nm);
    const double s = h->Integral(0, h->GetNbinsX() + 1);
    Check(s == nEv, Form("%s total incl. under/overflow (%.0f) == events", nm, s));
  }
  const char *vtxH[] = {"hn_vtx", "hN_vtx", "hnW_vtx", "hNW_vtx", "hZtrue_vtx"};
  for (const char *nm : vtxH)
  {
    TH1D *h = H1(nm);
    const double s = h->Integral(0, h->GetNbinsX() + 1);
    Check(s == nV, Form("%s total (%.0f) == vertexed events", nm, s));
  }
  Check(H1("hn_novtx")->Integral(0, -1) == nNV, "hn_novtx total == vertexless events");
  // central-z subsample: totals equal the tree count of |zTrue| < zcut
  {
    double zcut = 10.0;
    if (auto *cfgObj = dynamic_cast<TNamed *>(f->Get("config")))
    {
      const std::string c = cfgObj->GetTitle();
      const auto p = c.find("|z_true|<");
      if (p != std::string::npos) zcut = std::stod(c.substr(p + 9));
    }
    const long aZ = t->GetEntries(Form("abs(zTrue)<%g", zcut));
    const long vZ = t->GetEntries(Form("abs(zTrue)<%g && hasVtx", zcut));
    const char *zA[] = {"hn_all_zc", "hN_all_zc", "hnW_all_zc", "hNW_all_zc"};
    const char *zV[] = {"hn_vtx_zc", "hN_vtx_zc", "hnW_vtx_zc", "hNW_vtx_zc"};
    for (const char *nm : zA) Check(H1(nm)->Integral(0, -1) == aZ, Form("%s total == tree |zTrue|<%g count (%ld)", nm, zcut, aZ));
    for (const char *nm : zV) Check(H1(nm)->Integral(0, -1) == vZ, Form("%s total == tree vertexed |zTrue|<%g count (%ld)", nm, zcut, vZ));
  }
  // numerator <= denominator in every bin -> efficiencies in [0,1]
  const std::pair<const char *, const char *> pairs[] = {
      {"hn_vtx", "hn_all"}, {"hN_vtx", "hN_all"}, {"hnW_vtx", "hnW_all"}, {"hNW_vtx", "hNW_all"}, {"hZtrue_vtx", "hZtrue_all"},
      {"hn_vtx_zc", "hn_all_zc"}, {"hN_vtx_zc", "hN_all_zc"}, {"hnW_vtx_zc", "hnW_all_zc"}, {"hNW_vtx_zc", "hNW_all_zc"}};
  for (auto &p : pairs)
  {
    TH1D *a = H1(p.first), *b = H1(p.second);
    bool ok = a->GetNbinsX() == b->GetNbinsX();
    for (int i = 0; ok && i <= a->GetNbinsX() + 1; ++i)
      ok = a->GetBinContent(i) <= b->GetBinContent(i) && a->GetXaxis()->GetBinLowEdge(i) == b->GetXaxis()->GetBinLowEdge(i);
    Check(ok, Form("%s <= %s bin by bin, identical binning (efficiency in [0,1])", p.first, p.second));
  }
  // vertexed + vertexless = all, bin by bin
  {
    TH1D *a = H1("hn_all"), *v = H1("hn_vtx"), *nv = H1("hn_novtx");
    bool ok = true;
    for (int i = 0; i <= a->GetNbinsX() + 1; ++i) ok &= (v->GetBinContent(i) + nv->GetBinContent(i) == a->GetBinContent(i));
    Check(ok, "hn_vtx + hn_novtx == hn_all, bin by bin");
  }
  // response-matrix sums: projection on true N == N distribution; on reco n == n distribution
  auto checkProj = [&](const char *r, const char *hN, const char *hn) {
    TH2D *R = H2(r);
    TH1D *a = H1(hN), *b = H1(hn);
    bool okN = true, okn = true;
    for (int ix = 0; ix <= R->GetNbinsX() + 1; ++ix)
    {
      double s = 0;
      for (int iy = 0; iy <= R->GetNbinsY() + 1; ++iy) s += R->GetBinContent(ix, iy);
      okN &= (s == a->GetBinContent(ix));
    }
    for (int iy = 0; iy <= R->GetNbinsY() + 1; ++iy)
    {
      double s = 0;
      for (int ix = 0; ix <= R->GetNbinsX() + 1; ++ix) s += R->GetBinContent(ix, iy);
      okn &= (s == b->GetBinContent(iy));
    }
    Check(okN, Form("%s: sum over n for each N (row sum) == %s", r, hN));
    Check(okn, Form("%s: sum over N for each n == %s", r, hn));
  };
  checkProj("hResp_all", "hN_all", "hn_all");
  checkProj("hResp_vtx", "hN_vtx", "hn_vtx");
  // composition histograms sum to the number of good tracks
  {
    const double comp = H1("hComp_good_vtxEvt")->Integral() + H1("hComp_good_novtxEvt")->Integral();
    const double sil = H1("hSilicon_good_vtxEvt")->Integral() + H1("hSilicon_good_novtxEvt")->Integral();
    const double good = stat("good tracks");
    Check(comp == good, Form("truth composition sums to good tracks (%.0f vs %.0f)", comp, good));
    Check(sil == good, Form("silicon-content categories sum to good tracks (%.0f vs %.0f)", sil, good));
    Check(H2("hComp_vs_n")->Integral() == good, "composition-vs-n sums to good tracks");
    const double ptGood = H1("hPt_good_vtxEvt")->Integral(0, -1) + H1("hPt_good_novtxEvt")->Integral(0, -1);
    Check(ptGood == good, Form("per-track good pT histograms count every good track (%.0f)", ptGood));
    const double ptAll = H1("hPt_all_vtxEvt")->Integral(0, -1) + H1("hPt_all_novtxEvt")->Integral(0, -1);
    Check(ptAll == stat("tracks"), Form("per-track pT histograms count every track (%.0f vs %.0f)", ptAll, stat("tracks")));
  }
  Check(stat("n=0") + stat("n=1") + stat("n=2") + stat("n=3") + stat("n=4") + stat("n=5") ==
            H1("hn_all")->Integral(1, 6),
        "hStats n=0..5 == hn_all bins 0..5");
  Check(stat("n>N") + stat("n<N") + stat("n=N") == nEv, "n>N + n<N + n=N == events");

  // ---------------------------------------------------------------------------
  // D: NaN / Inf and unique names
  // ---------------------------------------------------------------------------
  Out("\n[D] NaN/Inf and object names");
  // hadd cannot merge TNamed objects, so a merged production holds one "config"
  // key per chunk. These repeats are allowed if every copy is identical apart
  // from its per-chunk "files=" count. An older cycle of a TTree (ROOT autosave
  // snapshot of a large tree, e.g. evt;1 next to the complete evt;2) is skipped and
  // reported; any other repeated name is an error.
  std::set<std::string> names, cfgTitles;
  long dupNames = 0, nCfg = 0, badBins = 0, nHist = 0, nOldTreeCycles = 0;
  for (TObject *o : *f->GetListOfKeys())
  {
    TKey *k = static_cast<TKey *>(o);
    if (std::string(k->GetClassName()) == "TTree")
    {
      const TKey *newest = f->GetKey(k->GetName());  // highest cycle of this name
      if (newest && k->GetCycle() < newest->GetCycle() && std::string(newest->GetClassName()) == "TTree")
      {
        Out(Form("  info: %s;%d is an autosave snapshot (%s;%d is the complete tree), skipped", k->GetName(), k->GetCycle(),
                 newest->GetName(), newest->GetCycle()));
        nOldTreeCycles++;
        continue;
      }
    }
    if (std::string(k->GetName()) == "config")
    {
      nCfg++;
      if (auto *c = dynamic_cast<TNamed *>(k->ReadObj()))
      {
        std::string t = c->GetTitle();
        t = t.substr(0, t.find("files="));
        cfgTitles.insert(t);
      }
      names.insert("config");
      continue;
    }
    if (!names.insert(k->GetName()).second) dupNames++;
    TObject *obj = k->ReadObj();
    if (TH1 *h = dynamic_cast<TH1 *>(obj))
    {
      nHist++;
      const int nb = (h->GetNbinsX() + 2) * (h->GetNbinsY() + 2) * (h->GetNbinsZ() + 2);
      for (int b = 0; b < nb; ++b)
        if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) badBins++;
    }
  }
  Check(dupNames == 0, Form("unique object names in output (%zu objects, %ld duplicates)", names.size(), dupNames));
  Check(nCfg >= 1 && cfgTitles.size() == 1,
        Form("%ld \"config\" copies (one per merged chunk), %zu distinct selection(s): %s", nCfg, cfgTitles.size(),
             cfgTitles.empty() ? "-" : cfgTitles.begin()->c_str()));
  Check(badBins == 0, Form("no NaN/Inf bin contents or errors in %ld histograms", nHist));

  // ---------------------------------------------------------------------------
  // E: new diagnostics (silicon, truth classification, response, efficiencies)
  // ---------------------------------------------------------------------------
  Out("\n[E] Combined counts, response matrices, diagnostics and confidence intervals");
  if (got.count(2642) && got.count(0) && got[2642].events == 1000 && got[0].events == 1000 && got.size() == 2)
  {
    Check(nEv == 2000 && nV == 429 && nNV == 1571,
          Form("combined: %.0f collisions, %.0f vertexed, %.0f vertexless (expected 2000/429/1571)", nEv, nV, nNV));
  }
  auto total2 = [&](const char *nm) { TH2D *h = H2(nm); return h->Integral(0, h->GetNbinsX() + 1, 0, h->GetNbinsY() + 1); };
  Check(total2("hResp_all") == nEv, Form("hResp_all integral incl. under/overflow (%.0f) == collisions", total2("hResp_all")));
  Check(total2("hResp_vtx") == nV, Form("hResp_vtx integral incl. under/overflow (%.0f) == vertexed", total2("hResp_vtx")));
  Check(total2("hResp_uprim") == nEv, "hResp_uprim (diagnostic) integral == collisions");
  Check(total2("hResp_uprimAcc") == nEv, "hResp_uprimAcc (diagnostic) integral == collisions");
  Check(total2("hEff2D_all") == nEv && total2("hEff2D_vtx") == nV, "2D (n, n_MVTX2) maps: totals == collisions / vertexed");
  {
    // P(n|N): every populated true-N column sums to 1 after normalisation
    TH2D *R = H2("hResp_all");
    int populated = 0, empty = 0;
    bool ok = true;
    for (int ix = 0; ix <= R->GetNbinsX() + 1; ++ix)
    {
      double s = 0;
      for (int iy = 0; iy <= R->GetNbinsY() + 1; ++iy) s += R->GetBinContent(ix, iy);
      if (s <= 0)
      {
        if (ix >= 1 && ix <= R->GetNbinsX()) empty++;
        continue;
      }
      populated++;
      double norm = 0;
      for (int iy = 0; iy <= R->GetNbinsY() + 1; ++iy) norm += R->GetBinContent(ix, iy) / s;
      ok &= std::fabs(norm - 1.) < 1e-12;
    }
    Check(ok, Form("P(n|N): all %d populated N columns sum to 1 (%d empty N columns in 0..%d, left empty)", populated,
                   empty, R->GetNbinsX() - 1));
  }
  {
    TH2D *a = H2("hEff2D_all"), *v = H2("hEff2D_vtx");
    bool ok = true;
    for (int i = 0; i <= a->GetNbinsX() + 1; ++i)
      for (int j = 0; j <= a->GetNbinsY() + 1; ++j) ok &= v->GetBinContent(i, j) <= a->GetBinContent(i, j);
    Check(ok, "2D map: vertexed <= all in every (n, n_MVTX2) cell");
  }
  {
    const char *mA[] = {"hnM1_all", "hnM2_all", "hnM3_all"};
    for (const char *nm : mA) Check(H1(nm)->Integral(0, -1) == nEv, Form("%s total == collisions", nm));
    const char *mV[] = {"hnM1_vtx", "hnM2_vtx", "hnM3_vtx"};
    for (const char *nm : mV) Check(H1(nm)->Integral(0, -1) == nV, Form("%s total == vertexed", nm));
    TH1D *a0 = H1("hn_all_m0"), *a1 = H1("hn_all_m1"), *a2 = H1("hn_all_m2"), *all = H1("hn_all");
    TH1D *v0 = H1("hn_vtx_m0"), *v1 = H1("hn_vtx_m1"), *v2 = H1("hn_vtx_m2"), *vall = H1("hn_vtx");
    bool okSum = true, okSub = true;
    for (int b = 0; b <= all->GetNbinsX() + 1; ++b)
    {
      okSum &= (a0->GetBinContent(b) + a1->GetBinContent(b) == all->GetBinContent(b)) &&
               (v0->GetBinContent(b) + v1->GetBinContent(b) == vall->GetBinContent(b));
      okSub &= a2->GetBinContent(b) <= a1->GetBinContent(b) && v2->GetBinContent(b) <= v1->GetBinContent(b) &&
               v0->GetBinContent(b) <= a0->GetBinContent(b) && v1->GetBinContent(b) <= a1->GetBinContent(b);
    }
    Check(okSum, "n_MVTX1 classes: (=0) + (>=1) == all, bin by bin, for all and vertexed");
    Check(okSub, "n_MVTX1 classes: (>=2) subset of (>=1); vertexed <= all in every class");
  }
  {
    // per-event consistency of the classifications (from the tree)
    int nGoodMVTX, nGoodTPCOnly, nINTTonly, nMVTX2, nMVTX3, nUniquePrim, nUniqueSec, nAddPrim, nAddSec, nUniquePrimAcc;
    t->SetBranchAddress("nGoodMVTX", &nGoodMVTX);
    t->SetBranchAddress("nGoodTPCOnly", &nGoodTPCOnly);
    t->SetBranchAddress("nINTTonly", &nINTTonly);
    t->SetBranchAddress("nMVTX2", &nMVTX2);
    t->SetBranchAddress("nMVTX3", &nMVTX3);
    t->SetBranchAddress("nUniquePrim", &nUniquePrim);
    t->SetBranchAddress("nUniqueSec", &nUniqueSec);
    t->SetBranchAddress("nAddPrim", &nAddPrim);
    t->SetBranchAddress("nAddSec", &nAddSec);
    t->SetBranchAddress("nUniquePrimAcc", &nUniquePrimAcc);
    long badSil = 0, badTruth = 0, badAcc = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (nGoodMVTX + nINTTonly + nGoodTPCOnly != n || nMVTX3 > nMVTX2 || nMVTX2 > nGoodMVTX) badSil++;
      if (nUniquePrim + nAddPrim != nGoodPrim || nUniqueSec + nAddSec != nGoodSec || nAddPrim + nAddSec != nGoodDup ||
          nUniquePrim + nUniqueSec + nAddPrim + nAddSec + nGoodUnm != n)
        badTruth++;
      if (nUniquePrimAcc > N || nUniquePrimAcc > nUniquePrim) badAcc++;
    }
    Check(badSil == 0, Form("per event: n_MVTX1 + n_INTTonly + n_TPConly == n and n_MVTX3 <= n_MVTX2 <= n_MVTX1 (%ld bad)", badSil));
    Check(badTruth == 0, Form("per event: unique + additional + unmatched == n, origin counts consistent (%ld bad)", badTruth));
    Check(badAcc == 0, Form("per event: unique in-acceptance primary matches <= N (%ld bad)", badAcc));
  }
  {
    TH2D *ou = H2("hOriginUnique");
    const double good = stat("good tracks");
    Check(ou->Integral() == good, Form("origin x uniqueness table sums to good tracks (%.0f)", ou->Integral()));
    const double unmOff = ou->GetBinContent(3, 1) + ou->GetBinContent(3, 2) + ou->GetBinContent(1, 3) + ou->GetBinContent(2, 3);
    Check(unmOff == 0, "unmatched tracks appear only in the (unmatched, unmatched) cell");
    Check(H1("hTrkClass_all")->Integral() == good, "exclusive track classes sum to good tracks");
    Check(H1("hTrkClass_N0")->GetBinContent(5) == 0,
          Form("N = 0 collisions contain no good track matched to an in-acceptance primary (%.0f)",
               H1("hTrkClass_N0")->GetBinContent(5)));
    Check(stat("uniquePrimAcc>N") == 0 && stat("dup G4 primary rows") == 0,
          "no duplicated G4-primary rows; unique in-acceptance matches never exceed N");
    Check(H1("hPurity_single")->GetBinContent(H1("hPurity_single")->GetNbinsX() + 1) == 0 &&
              H1("hPurity_multiAdd")->GetBinContent(H1("hPurity_multiAdd")->GetNbinsX() + 1) == 0,
          "truth purity nfromtruth/nhits never exceeds 1");
    const double addTracks = H1("hPurity_multiAdd")->GetEntries();
    Check(addTracks == ou->GetBinContent(1, 2) + ou->GetBinContent(2, 2) || H1("hPurity_multiAdd")->GetEntries() == 0,
          Form("additional-match diagnostics count every additional track (%.0f)", addTracks));
  }
  {
    // confidence intervals: lower <= eff <= upper within [0,1]; zero-success bins have upper > 0
    const std::pair<const char *, const char *> effs[] = {{"hn_vtx", "hn_all"},       {"hN_vtx", "hN_all"},
                                                          {"hnW_vtx", "hnW_all"},     {"hNW_vtx", "hNW_all"},
                                                          {"hZtrue_vtx", "hZtrue_all"}, {"hnM1_vtx", "hnM1_all"},
                                                          {"hnM2_vtx", "hnM2_all"},   {"hnM3_vtx", "hnM3_all"},
                                                          {"hn_vtx_m0", "hn_all_m0"}, {"hn_vtx_m1", "hn_all_m1"},
                                                          {"hn_vtx_m2", "hn_all_m2"}};
    bool ok = true;
    long zeroBins = 0;
    for (auto &p : effs)
    {
      TEfficiency e(*H1(p.first), *H1(p.second));
      e.SetStatisticOption(TEfficiency::kFCP);
      e.SetConfidenceLevel(0.683);
      for (int b = 1; b <= H1(p.second)->GetNbinsX(); ++b)
      {
        if (H1(p.second)->GetBinContent(b) <= 0) continue;
        const double eff = e.GetEfficiency(b), lo = eff - e.GetEfficiencyErrorLow(b), up = eff + e.GetEfficiencyErrorUp(b);
        ok &= (lo >= 0 && lo <= eff && eff <= up && up <= 1 && up > lo);
        if (H1(p.first)->GetBinContent(b) == 0)
        {
          zeroBins++;
          ok &= up > 0;
        }
      }
    }
    Check(ok, Form("Clopper-Pearson 68.3%%: 0 <= low <= eff <= up <= 1 in all populated bins; %ld zero-success bins have up > 0", zeroBins));
    // one-sided 95% upper limit used in the tables (k = 0): must equal 1 - 0.05^(1/N)
    bool okUL = true;
    for (int Nt : {1, 5, 41, 156, 303})
      okUL &= std::fabs(TEfficiency::ClopperPearson(Nt, 0, 0.90, true) - (1. - std::pow(0.05, 1. / Nt))) < 1e-9;
    Check(okUL, "table UL95 (k = 0) equals the analytic one-sided 95% limit 1 - 0.05^(1/N)");
  }

  Out(Form("\nSUMMARY: %d checks passed, %d failed", gPass, gFail));
  Out(gFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED - do not proceed to larger samples");
  gLog.close();
}
