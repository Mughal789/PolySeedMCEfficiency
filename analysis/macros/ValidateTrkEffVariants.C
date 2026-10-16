// =============================================================================
// ValidateTrkEffVariants.C - consistency checks across the efficiency variants
//
// Reads the 12 outputs of analysis/run_trkeff_variants.sh,
//   TrkEff_mvtx{0,2,3}_{ztrue,zreco}_{incl,track}.root,
// and checks:
//   R  the baseline variant (mvtx0, ztrue, incl) reproduces the Step-2 output
//      exactly (every histogram bin and error, every TEfficiency);
//   B  NumQ <= Num <= Den in every bin of every variant;
//   M  the denominator does not depend on the MVTX requirement, and
//      Num(MVTX3) <= Num(MVTX2) <= Num(MVTX0) (also for NumQ), bin by bin;
//   T  trackable <= inclusive for Den, Num and NumQ, bin by bin;
//   Z  the selected-event count depends only on the vertex mode;
//   L  every output is labelled with its vertex mode, MVTX cut and denominator;
//   G  (if PLOTDIR is given) every graph drawn by CompareTrkEffVariants.C equals the
//      variant's TEfficiency (points only for defined bins, values, Clopper-Pearson
//      errors) and TrkEffVariants_numbers.txt equals the integrated counts.
// The comparison with the original evaluator (event selection, no duplicates,
// D/N/NQ bin by bin) is done per variant by ValidateTrkEff.C.
//
// Usage:
//   root -l -b -q 'ValidateTrkEffVariants.C+("VARIANTDIR","Step2_TrkEff_baseline.root","report.log"[,"PLOTDIR"])'
// =============================================================================
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TNamed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <memory>
#include <string>
#include <vector>

namespace vv
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
const char *kTags[6] = {"pip", "pim", "kp", "km", "p", "pbar"};
const int kMvtx[3] = {0, 2, 3};
const char *kZ[2] = {"ztrue", "zreco"};
const char *kDen[2] = {"incl", "track"};

std::string Name(int m, const char *z, const char *d) { return Form("mvtx%d_%s_%s", m, z, d); }

// a <= b in every cell (including under/overflow)
bool LessEqual(const TH1 *a, const TH1 *b)
{
  for (int g = 0; g < a->GetNcells(); ++g)
    if (a->GetBinContent(g) > b->GetBinContent(g)) return false;
  return true;
}
bool Equal(const TH1 *a, const TH1 *b)
{
  if (a->GetNcells() != b->GetNcells()) return false;
  for (int g = 0; g < a->GetNcells(); ++g)
    if (a->GetBinContent(g) != b->GetBinContent(g) || std::fabs(a->GetBinError(g) - b->GetBinError(g)) > 1e-12) return false;
  return true;
}
}  // namespace vv

void ValidateTrkEffVariants(const char *dir, const char *step2File, const char *report, const char *plotDir = "")
{
  using namespace vv;
  gOut = fopen(report, "w");
  Log("=== ValidateTrkEffVariants ===");
  std::map<std::string, std::unique_ptr<TFile>> F;
  for (int m : kMvtx)
    for (const char *z : kZ)
      for (const char *d : kDen)
      {
        const std::string n = Name(m, z, d);
        F[n].reset(TFile::Open(Form("%s/TrkEff_%s.root", dir, n.c_str()), "READ"));
        Test("F all 12 variant outputs readable", F[n] && !F[n]->IsZombie(), n);
      }
  auto H = [&](const std::string &v, const char *kind, const char *tag, const char *proj) {
    return static_cast<TH1 *>(F[v]->Get(Form("%s_%s_%s", kind, tag, proj)));
  };

  // ---- R: baseline variant == Step 2 ------------------------------------------------
  {
    std::unique_ptr<TFile> S(TFile::Open(step2File, "READ"));
    TFile *B = F[Name(0, "ztrue", "incl")].get();
    TIter next(S->GetListOfKeys());
    TKey *key;
    while ((key = static_cast<TKey *>(next())))
    {
      TObject *a = key->ReadObj(), *b = B->Get(key->GetName());
      if (a->InheritsFrom(TH1::Class()))
        Test("R baseline variant == Step 2 (every histogram, bin and error)", b && Equal(static_cast<TH1 *>(a), static_cast<TH1 *>(b)),
             key->GetName());
      else if (a->InheritsFrom(TEfficiency::Class()))
      {
        auto *ea = static_cast<TEfficiency *>(a), *eb = static_cast<TEfficiency *>(b);
        Test("R baseline variant == Step 2 (every TEfficiency passed/total)",
             eb && Equal(ea->GetPassedHistogram(), eb->GetPassedHistogram()) && Equal(ea->GetTotalHistogram(), eb->GetTotalHistogram()),
             key->GetName());
      }
      else if (std::string(key->GetName()) == "config")
      {
        // intended differences only: keys added by later program versions, present on one side
        // and then at their "off" default (truth_min_tpc_layers = 0; v2: trk_nintt_min = 0,
        // trk_dca_xy_max = 0, trk_dca_z_max = 0); every key present on both sides must agree
        auto kv = [](const std::string &text) {
          std::map<std::string, std::string> m;
          std::istringstream in(text);
          std::string line;
          while (std::getline(in, line))
          {
            const auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            auto trim = [](std::string x) {
              x.erase(0, x.find_first_not_of(" \t"));
              x.erase(x.find_last_not_of(" \t") + 1);
              return x;
            };
            m[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
          }
          return m;
        };
        const std::map<std::string, std::string> offDefault = {
            {"truth_min_tpc_layers", "0"}, {"trk_nintt_min", "0"}, {"trk_dca_xy_max", "0"}, {"trk_dca_z_max", "0"}};
        const auto ma = kv(static_cast<TNamed *>(a)->GetTitle()), mb = kv(static_cast<TNamed *>(b)->GetTitle());
        std::string bad;
        for (const auto *pair : {&ma, &mb})
          for (const auto &[k, v] : *pair)
          {
            const auto &other = pair == &ma ? mb : ma;
            if (other.count(k)) { if (other.at(k) != v) bad = k; }
            else if (!offDefault.count(k) || offDefault.at(k) != v) bad = k;
          }
        Test("R baseline config == Step 2 config (keys added later only at their off default)", bad.empty(), bad);
      }
    }
  }

  for (const char *tag : kTags)
    for (const char *proj : {"ypt", "pt", "y"})
    {
      // ---- B: NumQ <= Num <= Den -----------------------------------------------------
      for (auto &kv : F)
        Test("B NumQ <= Num <= Den in every bin (all variants)",
             LessEqual(H(kv.first, "hNumQ", tag, proj), H(kv.first, "hNum", tag, proj)) &&
                 LessEqual(H(kv.first, "hNum", tag, proj), H(kv.first, "hDen", tag, proj)),
             kv.first + " " + tag + " " + proj);
      for (const char *z : kZ)
        for (const char *d : kDen)
        {
          // ---- M: MVTX ordering -------------------------------------------------------
          const std::string v0 = Name(0, z, d), v2 = Name(2, z, d), v3 = Name(3, z, d);
          Test("M denominator independent of the MVTX requirement",
               Equal(H(v0, "hDen", tag, proj), H(v2, "hDen", tag, proj)) && Equal(H(v0, "hDen", tag, proj), H(v3, "hDen", tag, proj)),
               v0 + " " + tag + " " + proj);
          for (const char *kind : {"hNum", "hNumQ"})
            Test(std::string("M ") + kind + "(MVTX3) <= " + kind + "(MVTX2) <= " + kind + "(MVTX0) in every bin",
                 LessEqual(H(v3, kind, tag, proj), H(v2, kind, tag, proj)) && LessEqual(H(v2, kind, tag, proj), H(v0, kind, tag, proj)),
                 v0 + " " + tag + " " + proj);
        }
      for (int m : kMvtx)
        for (const char *z : kZ)
          // ---- T: trackable <= inclusive -------------------------------------------------
          for (const char *kind : {"hDen", "hNum", "hNumQ"})
            Test(std::string("T trackable ") + kind + " <= inclusive " + kind + " in every bin",
                 LessEqual(H(Name(m, z, "track"), kind, tag, proj), H(Name(m, z, "incl"), kind, tag, proj)),
                 Name(m, z, "track") + " " + tag + " " + proj);
    }

  // ---- Z: selected events depend only on the vertex mode; L: labels -------------------
  std::map<std::string, double> selEv;
  for (auto &kv : F)
  {
    auto *st = static_cast<TH1 *>(kv.second->Get("hTrkEffStats"));
    selEv[kv.first] = st->GetBinContent(st->GetXaxis()->FindBin("selected events"));
    auto *sel = static_cast<TNamed *>(kv.second->Get("selection"));
    const std::string label = sel ? sel->GetTitle() : "";
    const bool reco = kv.first.find("zreco") != std::string::npos, track = kv.first.find("track") != std::string::npos;
    const std::string mv = kv.first.substr(4, 1);
    Test("L output labelled with vertex mode, MVTX cut and denominator",
         label.find(reco ? "reco vertex, |z_reco| < 8" : "|z_true| < 8") != std::string::npos &&
             label.find("MVTX >= " + mv) != std::string::npos &&
             label.find(track ? "trackable (>= 18 TPC layers)" : "inclusive") != std::string::npos,
         kv.first + ": " + label);
  }
  for (const char *z : kZ)
    for (auto &kv : selEv)
      if (kv.first.find(z) != std::string::npos)
        Test("Z selected events identical for all variants of a vertex mode", kv.second == selEv[Name(0, z, "incl")], kv.first);
  Test("Z reco-z selection differs from true-z selection", selEv[Name(0, "zreco", "incl")] != selEv[Name(0, "ztrue", "incl")]);

  // ---- G: comparison plots and numbers table ---------------------------------------------
  if (plotDir && *plotDir)
  {
    std::unique_ptr<TFile> P(TFile::Open(Form("%s/TrkEffVariants_plots.root", plotDir), "READ"));
    Test("G plot-object file readable", P && !P->IsZombie());
    long nGraphs = 0;
    for (auto &kv : F)
      for (const char *tag : kTags)
      {
        auto *g = dynamic_cast<TGraphAsymmErrors *>(P->Get(Form("g_%s_eff_%s_pt", kv.first.c_str(), tag)));
        if (!g) continue;  // only the variants used in the comparisons are drawn
        nGraphs++;
        auto *e = static_cast<TEfficiency *>(kv.second->Get(Form("eff_%s_pt", tag)));
        const TH1 *tot = e->GetTotalHistogram(), *pas = e->GetPassedHistogram();
        int defined = 0;
        for (int b = 1; b <= tot->GetNbinsX(); ++b) defined += tot->GetBinContent(b) > 0;
        bool ok = g->GetN() == defined;
        for (int i = 0; ok && i < g->GetN(); ++i)
        {
          const int b = tot->GetXaxis()->FindFixBin(g->GetX()[i]);
          const double d = tot->GetBinContent(b), n = pas->GetBinContent(b);
          ok = d > 0 && std::fabs(g->GetY()[i] - n / d) < 1e-9 &&
               std::fabs(g->GetErrorYlow(i) - (n / d - TEfficiency::ClopperPearson(d, n, 0.683, false))) < 1e-6 &&
               std::fabs(g->GetErrorYhigh(i) - (TEfficiency::ClopperPearson(d, n, 0.683, true) - n / d)) < 1e-6;
        }
        Test("G plotted graph == variant TEfficiency (defined bins only, values, CP errors)", ok, kv.first + " " + tag);
      }
    Test("G 5 distinct variants x 6 species drawn (MVTX0/2/3, trackable, reco-z)", nGraphs == 30, Form("%ld graphs", nGraphs));
    std::ifstream in(Form("%s/TrkEffVariants_numbers.txt", plotDir));
    std::string line;
    long nRows = 0;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string v, tag;
      double D, N, NQ;
      ss >> v >> tag >> D >> N >> NQ;
      nRows++;
      auto sum = [&](const char *kind) {
        const TH1 *h = H(v, kind, tag.c_str(), "ypt");
        return static_cast<const TH2 *>(h)->Integral(1, h->GetNbinsX(), 1, h->GetNbinsY());
      };
      Test("G numbers table D, N, NQ == in-range sums of the variant histograms",
           F.count(v) && D == sum("hDen") && N == sum("hNum") && NQ == sum("hNumQ"), v + " " + tag);
    }
    Test("G numbers table has 12 variants x 6 species", nRows == 72, Form("%ld rows", nRows));
  }

  Log("");
  Log("--- checks (tests / failures) ---");
  long nFail = 0;
  for (auto &kv : gChk)
  {
    const bool ok = kv.second.second == 0;
    if (!ok) nFail++;
    Log(Form("%s  %-86s %5ld tests %4ld failed%s", ok ? "PASS" : "FAIL", kv.first.c_str(), kv.second.first, kv.second.second,
             ok ? "" : Form("  (first: %s)", gFirst[kv.first].c_str())));
  }
  Log(Form("selected events: true-z %.0f, reco-z %.0f", selEv[Name(0, "ztrue", "incl")], selEv[Name(0, "zreco", "incl")]));
  Log("");
  Log(Form("SUMMARY: %zu checks, %ld failed", gChk.size(), nFail));
  Log(nFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  fclose(gOut);
  gOut = nullptr;
}
