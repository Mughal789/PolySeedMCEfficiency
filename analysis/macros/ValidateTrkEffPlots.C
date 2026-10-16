// =============================================================================
// ValidateTrkEffPlots.C - checks that the plots show exactly the Step-2 results
//
// Compares every object drawn by PlotTrkEff.C (stored in
// TrackEfficiency_plots.root) and every number it printed
// (TrackEfficiency_numbers.txt, LowPtDiagnostics.txt) with the TrkEffAnalysis
// output and the Step-2 text summary. Clopper-Pearson intervals are recomputed
// independently from the raw counts.
//
// Usage:
//   root -l -b -q 'ValidateTrkEffPlots.C+("TrkEff_baseline.root","TrkEff_baseline_summary.txt","PLOTDIR","report.log")'
// =============================================================================
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace vp
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
bool Close(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }
const char *kTags[6] = {"pip", "pim", "kp", "km", "p", "pbar"};
}  // namespace vp

void ValidateTrkEffPlots(const char *trkEffFile, const char *step2Summary, const char *plotDir, const char *report)
{
  using namespace vp;
  gOut = fopen(report, "w");
  Log("=== ValidateTrkEffPlots ===");
  TFile *R = TFile::Open(trkEffFile, "READ");
  TFile *P = TFile::Open(Form("%s/TrackEfficiency_plots.root", plotDir), "READ");
  Test("I plot-object file readable", P && !P->IsZombie());

  // ---- 1D graphs: one point per defined bin, value and CP errors -------------------
  for (const char *tag : kTags)
    for (const char *kind : {"eff", "effQ"})
      for (const char *proj : {"pt", "y"})
      {
        const std::string name = Form("%s_%s_%s", kind, tag, proj);
        auto *e = dynamic_cast<TEfficiency *>(R->Get(name.c_str()));
        auto *g = dynamic_cast<TGraphAsymmErrors *>(P->Get(("g_" + name).c_str()));
        const bool plotted = g != nullptr;
        // eff_*_pt / eff_*_y are all plotted (pages 5, 6), effQ_*_pt on page 7; effQ_*_y is not plotted
        if (std::string(kind) == "effQ" && std::string(proj) == "y")
        {
          Test("G effQ vs y is not plotted (by design)", !plotted, name);
          continue;
        }
        Test("G every plotted 1D efficiency exists", plotted, name);
        if (!plotted) continue;
        const TH1 *tot = e->GetTotalHistogram(), *pas = e->GetPassedHistogram();
        int defined = 0;
        for (int b = 1; b <= tot->GetNbinsX(); ++b) defined += tot->GetBinContent(b) > 0;
        Test("G number of points == number of defined bins (no point for empty bins)", g->GetN() == defined,
             Form("%s: %d vs %d", name.c_str(), g->GetN(), defined));
        for (int i = 0; i < g->GetN(); ++i)
        {
          const int b = tot->GetXaxis()->FindFixBin(g->GetX()[i]);
          const double d = tot->GetBinContent(b), n = pas->GetBinContent(b);
          Test("G point lies in a bin with total > 0", d > 0, name);
          if (d <= 0) continue;
          Test("G point value == passed/total of the Step-2 histograms", Close(g->GetY()[i], n / d), name);
          Test("G x error == bin half width", Close(g->GetErrorXlow(i), tot->GetBinWidth(b) / 2) &&
                                                   Close(g->GetErrorXhigh(i), tot->GetBinWidth(b) / 2), name);
          const double lo = TEfficiency::ClopperPearson(d, n, 0.683, false), up = TEfficiency::ClopperPearson(d, n, 0.683, true);
          Test("G y errors == Clopper-Pearson 68.3% recomputed from the counts",
               Close(g->GetErrorYlow(i), n / d - lo, 1e-6) && Close(g->GetErrorYhigh(i), up - n / d, 1e-6),
               Form("%s point %d", name.c_str(), i));
        }
      }

  // ---- 2D maps: colour value per defined bin; undefined bins counted ----------------
  std::map<std::string, int> undefinedYpt;
  for (const char *tag : kTags)
    for (const char *kind : {"eff", "effQ"})
    {
      const std::string name = Form("%s_%s_ypt", kind, tag);
      auto *e = dynamic_cast<TEfficiency *>(R->Get(name.c_str()));
      auto *m = dynamic_cast<TH2 *>(P->Get(("map_" + name).c_str()));
      Test("M every 2D map exists", m != nullptr, name);
      if (!m) continue;
      const TH2 *tot = static_cast<const TH2 *>(e->GetTotalHistogram()), *pas = static_cast<const TH2 *>(e->GetPassedHistogram());
      Test("M map binning == Step-2 binning", m->GetNbinsX() == tot->GetNbinsX() && m->GetNbinsY() == tot->GetNbinsY());
      int undef = 0;
      for (int ix = 1; ix <= tot->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= tot->GetNbinsY(); ++iy)
        {
          const double d = tot->GetBinContent(ix, iy), n = pas->GetBinContent(ix, iy);
          Test("M bin edges identical", Close(m->GetXaxis()->GetBinLowEdge(ix), tot->GetXaxis()->GetBinLowEdge(ix)) &&
                                            Close(m->GetYaxis()->GetBinLowEdge(iy), tot->GetYaxis()->GetBinLowEdge(iy)));
          if (d > 0) Test("M defined bin colour value == passed/total", Close(m->GetBinContent(ix, iy), n / d), name);
          else
          {
            undef++;
            // drawn content is 0, which PlotTrkEff covers with a hatched box (MaskUndefinedBins:
            // exactly the bins with total = 0, the same condition as here)
            Test("M undefined bin carries no efficiency value (content 0, masked)", m->GetBinContent(ix, iy) == 0, name);
          }
        }
      undefinedYpt[name] = undef;
    }

  // ---- numbers file vs Step-2 summary and recomputation ------------------------------
  std::map<std::string, std::vector<double>> num;  // tag -> D N NQ eff lo up effQ lo up empty_ypt empty_pt empty_y nbins
  {
    std::ifstream in(Form("%s/TrackEfficiency_numbers.txt", plotDir));
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      double x;
      while (ss >> x) num[tag].push_back(x);
    }
  }
  std::map<std::string, std::vector<double>> sum2;  // Step-2 summary: tag -> D N NQ eps epsQ
  {
    std::ifstream in(step2Summary);
    std::string line;
    while (std::getline(in, line))
    {
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      bool known = false;
      for (const char *t : kTags) known = known || tag == t;
      if (!known || line.find("+-") == std::string::npos) continue;
      double D, N, NQ, eps, de, epsq;
      std::string pm;
      ss >> D >> N >> NQ >> eps >> pm >> de >> epsq;
      sum2[tag] = {D, N, NQ, eps, epsq};
    }
  }
  for (const char *tag : kTags)
  {
    const auto &v = num[tag];
    Test("N numbers file has a complete row per species", v.size() == 13, tag);
    if (v.size() != 13) continue;
    auto *e = static_cast<TEfficiency *>(R->Get(Form("eff_%s_ypt", tag)));
    const TH1 *tot = e->GetTotalHistogram(), *pas = e->GetPassedHistogram();
    const TH1 *pasQ = static_cast<TEfficiency *>(R->Get(Form("effQ_%s_ypt", tag)))->GetPassedHistogram();
    double D = 0, N = 0, NQ = 0;
    for (int ix = 1; ix <= tot->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= tot->GetNbinsY(); ++iy)
      {
        const int b = tot->GetBin(ix, iy);
        D += tot->GetBinContent(b);
        N += pas->GetBinContent(b);
        NQ += pasQ->GetBinContent(b);
      }
    Test("N D, N, NQ == sums of the Step-2 histograms", v[0] == D && v[1] == N && v[2] == NQ, tag);
    // The Step-2 summary totals include particles outside the (y, pT) range (under/overflow
    // cells, e.g. pT > last edge); the plotted numbers are in-range only. Compare the summary
    // with the all-cell sums of the same histograms, and report the out-of-range particles.
    const int nx = tot->GetNbinsX(), ny = tot->GetNbinsY();
    const double Dall = static_cast<const TH2 *>(tot)->Integral(0, nx + 1, 0, ny + 1),
                 Nall = static_cast<const TH2 *>(pas)->Integral(0, nx + 1, 0, ny + 1),
                 NQall = static_cast<const TH2 *>(pasQ)->Integral(0, nx + 1, 0, ny + 1);
    if (Dall != D) Log(Form("  info: %s has %.0f particle(s) outside the (y, pT) range (in summary totals, not in the plots)", tag, Dall - D));
    Test("N D, N, NQ == Step-2 report (summary, all cells incl. out-of-range)",
         sum2.count(tag) && Dall == sum2[tag][0] && Nall == sum2[tag][1] && NQall == sum2[tag][2], tag);
    Test("N summary eff, effQ == all-cell ratio (3 decimals)",
         sum2.count(tag) && Dall > 0 && std::fabs(Nall / Dall - sum2[tag][3]) < 5e-4 && std::fabs(NQall / Dall - sum2[tag][4]) < 5e-4, tag);
    Test("N plotted integrated eff == in-range ratio", D > 0 && std::fabs(v[3] - N / D) < 5e-4, tag);
    Test("N integrated CP interval recomputed", Close(v[4], TEfficiency::ClopperPearson(D, N, 0.683, false), 1e-3) &&
                                                     Close(v[5], TEfficiency::ClopperPearson(D, N, 0.683, true), 1e-3), tag);
    Test("N empty (y,pT) bins == masked bins in the map", v[9] == undefinedYpt[Form("eff_%s_ypt", tag)], tag);
    // projections stored in Step 2 == 2D sums (the plotted pT/y efficiencies use them)
    auto *dPt = static_cast<TH1 *>(R->Get(Form("hDen_%s_pt", tag)));
    auto *dY = static_cast<TH1 *>(R->Get(Form("hDen_%s_y", tag)));
    bool ok = true;
    for (int b = 1; b <= dPt->GetNbinsX(); ++b) ok = ok && dPt->GetBinContent(b) == static_cast<const TH2 *>(tot)->Integral(1, tot->GetNbinsX(), b, b);
    for (int b = 1; b <= dY->GetNbinsX(); ++b) ok = ok && dY->GetBinContent(b) == static_cast<const TH2 *>(tot)->Integral(b, b, 1, tot->GetNbinsY());
    Test("P stored projections == 2D sums", ok, tag);
  }

  // ---- low-pT diagnostics: "found" and totals == Step-2 numerator / denominator ------
  {
    std::ifstream in(Form("%s/LowPtDiagnostics.txt", plotDir));
    if (in)
    {
      const std::map<std::string, std::vector<std::string>> groups = {
          {"pions", {"pip", "pim"}}, {"kaons", {"kp", "km"}}, {"protons", {"p", "pbar"}}};
      std::string line, group;
      int bin = 0;
      while (std::getline(in, line))
      {
        if (line.find(": categories") != std::string::npos)
        {
          group = line.substr(0, line.find(':'));
          bin = 0;
          continue;
        }
        if (line.find("N=") == std::string::npos) continue;
        bin++;
        const long total = std::stol(line.substr(line.find("N=") + 2));
        std::string rest = line.substr(line.find('|') + 1);
        const long found = std::stol(rest);
        double D = 0, N = 0;
        for (const auto &t : groups.at(group))
        {
          auto *hD = static_cast<TH2 *>(R->Get(Form("hDen_%s_ypt", t.c_str())));
          auto *hN = static_cast<TH2 *>(R->Get(Form("hNum_%s_ypt", t.c_str())));
          D += hD->Integral(0, hD->GetNbinsX() + 1, bin, bin);
          N += hN->Integral(0, hN->GetNbinsX() + 1, bin, bin);
        }
        Test("L low-pT category 'found' == Step-2 numerator (per group, pT bin)", found == N,
             Form("%s bin %d: %ld vs %g", group.c_str(), bin, found, N));
        Test("L low-pT total == Step-2 denominator (per group, pT bin)", total == D,
             Form("%s bin %d: %ld vs %g", group.c_str(), bin, total, D));
      }
    }
    else
      Test("L LowPtDiagnostics.txt present", false);
  }

  Log("");
  Log("--- checks (tests / failures) ---");
  long nFail = 0;
  for (auto &kv : gChk)
  {
    const bool ok = kv.second.second == 0;
    if (!ok) nFail++;
    Log(Form("%s  %-82s %6ld tests %5ld failed%s", ok ? "PASS" : "FAIL", kv.first.c_str(), kv.second.first, kv.second.second,
             ok ? "" : Form("  (first: %s)", gFirst[kv.first].c_str())));
  }
  Log("");
  Log(Form("SUMMARY: %zu checks, %ld failed", gChk.size(), nFail));
  Log(nFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  fclose(gOut);
  gOut = nullptr;
}
