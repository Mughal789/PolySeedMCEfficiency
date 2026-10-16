// =============================================================================
// MakeFinalTrkEffMaps.C - per-species tracking-efficiency maps in the analysis
// windows, from the fine-grid raw counts of a TrkEffAnalysis output
// (config/trkeff_final_mvtx3.cfg), plus the Corral hook format and plots.
//
// Species windows (true y, true pT) and map binning (subsets of the fine grid):
//   pi+-  : |y| < 1.0 in 0.2 steps;  pT 0.10 ... 5.00 GeV/c (14 bins)
//   K+-   : |y| < 0.7 (-0.7,-0.4,0,0.4,0.7); pT 0.1-0.6 in 0.1 steps; the bin
//           0.6-1.0 is a QA bin (the dE/dx PID of the analysis stops at p = 0.6)
//   p,pbar: |y| < 0.6 (-0.6,-0.3,0,0.3,0.6); pT 0.4, 0.5, 0.6, 0.7, 0.8, 1.0
// Every map edge must be an edge of the fine grid (checked; otherwise the
// macro stops), so the rebinning is an exact sum of raw counts.
//
// Output OUTDIR/TrkEffFinalMaps.root, per species <sp> (pip pim kp km p pbar):
//   hDen_<sp>_map, hNum_<sp>_map, hNumQ_<sp>_map   TH2D raw counts (x = y_true, y = pT_true)
//   eff_<sp>_map, effQ_<sp>_map                    TEfficiency, Clopper-Pearson 68.3%
//   eff_<sp>_pt, eff_<sp>_y                        1D TEfficiency inside the window
//   heffpty_<sp>                                   Corral hook format: TH2D(y, pT) with
//                                                  content = efficiency, error = mean of the
//                                                  CP half-intervals; content = -1 where D = 0
//   heffpty_<sp>_errlo, heffpty_<sp>_errhi         CP lower / upper half-intervals
//   window_<sp>, selection, source                 TNamed
// Optional reference (e.g. the same selection with MVTX >= 2): its maps are
// written with the suffix _ref and compared on the plots.
// Checks (OUTDIR/maps_validation.log): edges on the fine grid; rebinned sums
// equal the fine-bin sums in the window; hook content = passed/total; NQ <= N <= D.
//
// Usage: root -l -b -q 'MakeFinalTrkEffMaps.C+("TrkEff_final_mvtx3.root","OUTDIR","label"[,"ref.root","ref label"])'
// =============================================================================
#include "TrkEffPlotStyle.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace finmaps
{
struct SpeciesWindow
{
  const char *sp, *label, *text;
  std::vector<double> y, pt;
  int color, marker;
  double ptAnaMax = 99;  // integrated numbers use pT < ptAnaMax (kaons: QA bin above 0.6 excluded)
};
const std::vector<SpeciesWindow> &Windows()
{
  static const std::vector<SpeciesWindow> w = {
      {"pip", "#pi^{+}", "|y| < 1.0, pT >= 0.1 GeV/c", {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0},
       {0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 3.00, 5.00}, kAzure + 2, 20},
      {"pim", "#pi^{-}", "|y| < 1.0, pT >= 0.1 GeV/c", {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0},
       {0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 3.00, 5.00}, kAzure + 2, 24},
      {"kp", "K^{+}", "|y| < 0.7, 0.1 <= pT < 0.6 GeV/c (0.6-1.0: QA bin)", {-0.7, -0.4, 0.0, 0.4, 0.7},
       {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0}, kRed + 1, 21, 0.6},
      {"km", "K^{-}", "|y| < 0.7, 0.1 <= pT < 0.6 GeV/c (0.6-1.0: QA bin)", {-0.7, -0.4, 0.0, 0.4, 0.7},
       {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0}, kRed + 1, 25, 0.6},
      {"p", "p", "|y| < 0.6, 0.4 <= pT < 1.0 GeV/c", {-0.6, -0.3, 0.0, 0.3, 0.6}, {0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
       kGreen + 2, 22},
      {"pbar", "#bar{p}", "|y| < 0.6, 0.4 <= pT < 1.0 GeV/c", {-0.6, -0.3, 0.0, 0.3, 0.6}, {0.4, 0.5, 0.6, 0.7, 0.8, 1.0},
       kGreen + 2, 26}};
  return w;
}

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

// index of edge e on axis a (exact grid edge), or -1
int EdgeIndex(const TAxis *a, double e)
{
  for (int i = 1; i <= a->GetNbins() + 1; ++i)
    if (std::fabs(a->GetBinLowEdge(i) - e) < 1e-6) return i;
  return -1;
}

// exact rebin of a fine (y, pT) histogram onto the window edges
TH2D *Rebin(const TH2D *fine, const SpeciesWindow &w, const char *name)
{
  auto *h = new TH2D(name, Form("%s;y_{true};p_{T,true} (GeV/c)", fine->GetTitle()), w.y.size() - 1, w.y.data(),
                     w.pt.size() - 1, w.pt.data());
  h->SetDirectory(nullptr);
  h->Sumw2();
  for (size_t i = 0; i + 1 < w.y.size(); ++i)
    for (size_t j = 0; j + 1 < w.pt.size(); ++j)
    {
      const int x1 = EdgeIndex(fine->GetXaxis(), w.y[i]), x2 = EdgeIndex(fine->GetXaxis(), w.y[i + 1]) - 1;
      const int y1 = EdgeIndex(fine->GetYaxis(), w.pt[j]), y2 = EdgeIndex(fine->GetYaxis(), w.pt[j + 1]) - 1;
      double err = 0;
      const double c = fine->IntegralAndError(x1, x2, y1, y2, err);
      h->SetBinContent(i + 1, j + 1, c);
      h->SetBinError(i + 1, j + 1, err);
    }
  return h;
}

TEfficiency *Eff(const TH1 &pass, const TH1 &tot, const char *name)
{
  auto *e = new TEfficiency(pass, tot);
  e->SetName(name);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(0.683);
  e->SetDirectory(nullptr);
  return e;
}

// integrated efficiency inside the analysis window (pT < ptAnaMax), Clopper-Pearson 68.3%
trkeffplot::Integrated IntegrateWindow(const TEfficiency *e, double ptMax)
{
  const auto *t = static_cast<const TH2 *>(e->GetTotalHistogram()), *p = static_cast<const TH2 *>(e->GetPassedHistogram());
  trkeffplot::Integrated r;
  for (int i = 1; i <= t->GetNbinsX(); ++i)
    for (int j = 1; j <= t->GetNbinsY(); ++j)
      if (t->GetYaxis()->GetBinUpEdge(j) <= ptMax + 1e-9)
      {
        r.total += t->GetBinContent(i, j);
        r.passed += p->GetBinContent(i, j);
      }
  if (r.total > 0)
  {
    r.eff = r.passed / r.total;
    r.low = TEfficiency::ClopperPearson(r.total, r.passed, 0.683, false);
    r.up = TEfficiency::ClopperPearson(r.total, r.passed, 0.683, true);
  }
  return r;
}

struct SpeciesMaps
{
  TH2D *den = nullptr, *num = nullptr, *numQ = nullptr;
  TEfficiency *map = nullptr, *mapQ = nullptr, *pt = nullptr, *y = nullptr;
  TH2D *hook = nullptr, *lo = nullptr, *hi = nullptr;
};

SpeciesMaps Build(TFile &in, const SpeciesWindow &w, const std::string &suffix, bool check)
{
  SpeciesMaps m;
  const TH2D *fD = in.Get<TH2D>(Form("hDen_%s_ypt", w.sp)), *fN = in.Get<TH2D>(Form("hNum_%s_ypt", w.sp)),
             *fQ = in.Get<TH2D>(Form("hNumQ_%s_ypt", w.sp));
  if (!fD || !fN || !fQ)
  {
    Out(Form("FATAL: %s histograms missing in %s", w.sp, in.GetName()));
    std::exit(2);
  }
  bool edgesOk = true;
  for (double e : w.y) edgesOk = edgesOk && EdgeIndex(fD->GetXaxis(), e) > 0;
  for (double e : w.pt) edgesOk = edgesOk && EdgeIndex(fD->GetYaxis(), e) > 0;
  if (check) Check(edgesOk, Form("%s%s map edges are edges of the fine grid", w.sp, suffix.c_str()));
  if (!edgesOk)
  {
    Out("FATAL: map edges not on the fine grid - the rebinning would not be exact");
    std::exit(2);
  }
  const std::string s = std::string(w.sp) + suffix;
  m.den = Rebin(fD, w, Form("hDen_%s_map", s.c_str()));
  m.num = Rebin(fN, w, Form("hNum_%s_map", s.c_str()));
  m.numQ = Rebin(fQ, w, Form("hNumQ_%s_map", s.c_str()));
  m.map = Eff(*m.num, *m.den, Form("eff_%s_map", s.c_str()));
  m.mapQ = Eff(*m.numQ, *m.den, Form("effQ_%s_map", s.c_str()));
  std::unique_ptr<TH1D> dPt(m.den->ProjectionY("_dpt", 1, m.den->GetNbinsX())), nPt(m.num->ProjectionY("_npt", 1, m.num->GetNbinsX()));
  std::unique_ptr<TH1D> dY(m.den->ProjectionX("_dy", 1, m.den->GetNbinsY())), nY(m.num->ProjectionX("_ny", 1, m.num->GetNbinsY()));
  m.pt = Eff(*nPt, *dPt, Form("eff_%s_pt", s.c_str()));
  m.y = Eff(*nY, *dY, Form("eff_%s_y", s.c_str()));
  // hook format
  m.hook = static_cast<TH2D *>(m.den->Clone(Form("heffpty_%s", s.c_str())));
  m.lo = static_cast<TH2D *>(m.den->Clone(Form("heffpty_%s_errlo", s.c_str())));
  m.hi = static_cast<TH2D *>(m.den->Clone(Form("heffpty_%s_errhi", s.c_str())));
  for (TH2D *h : {m.hook, m.lo, m.hi}) h->Reset();
  m.hook->SetTitle(Form("%s tracking efficiency (content -1 = undefined, no particles);y_{true};p_{T,true} (GeV/c)", w.label));
  long badHook = 0, badOrder = 0, undefined = 0;
  for (int i = 1; i <= m.den->GetNbinsX(); ++i)
    for (int j = 1; j <= m.den->GetNbinsY(); ++j)
    {
      const double d = m.den->GetBinContent(i, j), n = m.num->GetBinContent(i, j), q = m.numQ->GetBinContent(i, j);
      badOrder += !(q <= n && n <= d);
      if (d <= 0)
      {
        ++undefined;
        m.hook->SetBinContent(i, j, -1);
        continue;
      }
      const int b = m.map->GetGlobalBin(i, j);
      const double e = n / d, el = m.map->GetEfficiencyErrorLow(b), eh = m.map->GetEfficiencyErrorUp(b);
      m.hook->SetBinContent(i, j, e);
      m.hook->SetBinError(i, j, 0.5 * (el + eh));
      m.lo->SetBinContent(i, j, el);
      m.hi->SetBinContent(i, j, eh);
      badHook += std::fabs(m.map->GetEfficiency(b) - e) > 1e-12;
    }
  if (check)
  {
    // rebinned totals = fine sums inside the window
    const int x1 = EdgeIndex(fD->GetXaxis(), w.y.front()), x2 = EdgeIndex(fD->GetXaxis(), w.y.back()) - 1;
    const int y1 = EdgeIndex(fD->GetYaxis(), w.pt.front()), y2 = EdgeIndex(fD->GetYaxis(), w.pt.back()) - 1;
    const bool sumsOk = fD->Integral(x1, x2, y1, y2) == m.den->Integral() && fN->Integral(x1, x2, y1, y2) == m.num->Integral() &&
                        fQ->Integral(x1, x2, y1, y2) == m.numQ->Integral();
    Check(sumsOk, Form("%s%s rebinned D, N, NQ = fine-grid sums in the window (D %.0f, N %.0f)", w.sp, suffix.c_str(),
                       m.den->Integral(), m.num->Integral()));
    Check(badHook == 0 && badOrder == 0,
          Form("%s%s hook content = reconstructed/generated; reco. (correct charge) <= reco. <= generated (%ld undefined bins, content -1)", w.sp, suffix.c_str(),
               undefined));
  }
  return m;
}
}  // namespace finmaps

void MakeFinalTrkEffMaps(const char *inPath, const char *outDir, const char *sampleLabel, const char *refPath = "",
                         const char *refLabel = "reference")
{
  using namespace finmaps;
  using namespace trkeffplot;
  const std::string dir = outDir;
  gSystem->mkdir((dir + "/figures").c_str(), true);
  gLog = fopen((dir + "/maps_validation.log").c_str(), "w");
  TH1::AddDirectory(false);
  std::unique_ptr<TFile> in(TFile::Open(inPath, "READ"));
  if (!in || in->IsZombie())
  {
    Out(std::string("FATAL: cannot open ") + inPath);
    return;
  }
  std::unique_ptr<TFile> ref;
  if (refPath && *refPath) ref.reset(TFile::Open(refPath, "READ"));
  const std::string selection = in->Get<TNamed>("selection") ? in->Get<TNamed>("selection")->GetTitle() : "?";
  Out("=== MakeFinalTrkEffMaps ===");
  Out(std::string("input    : ") + inPath + "\nselection: " + selection);
  if (ref) Out(std::string("reference: ") + refPath + " (" + refLabel + ")");

  std::map<std::string, SpeciesMaps> M, R;
  for (const auto &w : Windows())
  {
    M[w.sp] = Build(*in, w, "", true);
    if (ref) R[w.sp] = Build(*ref, w, "_ref", true);
  }

  // ---- ROOT output ---------------------------------------------------------------
  TFile out((dir + "/TrkEffFinalMaps.root").c_str(), "RECREATE");
  TNamed("selection", selection.c_str()).Write();
  TNamed("source", inPath).Write();
  TNamed("hook_convention",
         "heffpty_<sp>: TH2D x = y_true, y = pT_true; content = efficiency (Clopper-Pearson central value passed/total), "
         "error = mean of the 68.3% CP half-intervals (heffpty_<sp>_errlo/_errhi); content = -1 where no generated particle "
         "(undefined). Corral (commented hook) reads heffpty_cent[ParticleID][cent]->FindBin(y, min(pT, MaxPtWithEff-0.05)).")
      .Write();
  for (const auto &w : Windows())
  {
    TNamed(Form("window_%s", w.sp), w.text).Write();
    for (auto *mm : {&M[w.sp], ref ? &R[w.sp] : nullptr})
    {
      if (!mm) continue;
      for (TObject *o : std::vector<TObject *>{mm->den, mm->num, mm->numQ, mm->map, mm->mapQ, mm->pt, mm->y, mm->hook, mm->lo, mm->hi})
        o->Write();
    }
  }
  out.Close();

  // ---- integrated numbers --------------------------------------------------------
  Out("");
  Out("  integrated in the analysis window (Clopper-Pearson 68.3%):");
  for (const auto &w : Windows())
  {
    const Integrated I = IntegrateWindow(M[w.sp].map, w.ptAnaMax), IQ = IntegrateWindow(M[w.sp].mapQ, w.ptAnaMax);
    std::string line = Form("    %-5s %-52s eps %s   epsQ %.3f", w.sp, w.text, FormatIntegrated(I).c_str(), IQ.eff);
    if (ref) line += Form("   | %s %.3f", refLabel, IntegrateWindow(R[w.sp].map, w.ptAnaMax).eff);
    Out(line);
  }

  // ---- plots ---------------------------------------------------------------------
  SetTrkEffStyle();
  std::vector<TCanvas *> pages;
  auto page = [&](const char *name, int nx, int ny, int w = 1400, int h = 650) {
    auto *c = new TCanvas(name, name, w, h);
    c->Divide(nx, ny);
    pages.push_back(c);
    return c;
  };
  // page 1: summary
  {
    auto *c = new TCanvas("summary", "summary", 1100, 800);
    pages.push_back(c);
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.030);
    double y = 0.93;
    t.DrawLatex(0.05, y, "#bf{#it{sPHENIX}} Simulation   final tracking-efficiency maps");
    y -= 0.05;
    t.DrawLatex(0.05, y, sampleLabel);
    y -= 0.05;
    t.SetTextSize(0.022);
    {
      // the selection label is long: one line per ';'-separated part
      std::string rest = selection;
      bool first = true;
      while (!rest.empty())
      {
        const auto cut = rest.find(';');
        std::string part = rest.substr(0, cut);
        rest = cut == std::string::npos ? "" : rest.substr(cut + 1);
        t.DrawLatex(0.05, y, ((first ? "selection: " : "                 ") + part).c_str());
        first = false;
        y -= 0.035;
      }
    }
    t.DrawLatex(0.05, y, "#varepsilon_{trk} = P(selected track | generated primary, final reco vertex, |z_{reco}| < 8 cm); vertex and PID efficiency NOT included");
    y -= 0.06;
    t.SetTextSize(0.026);
    for (const auto &w : Windows())
    {
      const Integrated I = IntegrateWindow(M[w.sp].map, w.ptAnaMax);
      std::string s = Form("%s:  %s   #rightarrow   %s", w.label, w.text, FormatIntegrated(I).c_str());
      if (ref) s += Form("   (%s: %.3f)", refLabel, IntegrateWindow(R[w.sp].map, w.ptAnaMax).eff);
      t.DrawLatex(0.05, y, s.c_str());
      y -= 0.055;
    }
    t.SetTextSize(0.020);
    t.DrawLatex(0.05, 0.08, "Maps: colour = efficiency; text = reconstructed/generated per bin; hatched grey = no generated particle (undefined).");
    t.DrawLatex(0.05, 0.05, "Errors: Clopper-Pearson 68.3%. No smoothing, fitting or interpolation.");
    t.DrawLatex(0.05, 0.11, "Numbers above: analysis window only. Legends of the 1D plots integrate all plotted bins (kaons: incl. the 0.6-1.0 GeV/c QA bin).");
  }
  // pages 2-4: maps per charge pair
  for (int pair = 0; pair < 3; ++pair)
  {
    auto *c = page(Form("maps_%d", pair), 2, 1);
    for (int k = 0; k < 2; ++k)
    {
      const SpeciesWindow &w = Windows()[2 * pair + k];
      c->cd(k + 1);
      gPad->SetRightMargin(0.18);
      gPad->SetLogy(pair == 0);
      TH2D *m = EfficiencyMap(M[w.sp].map, Form("map_%s", w.sp));
      m->SetTitle(Form(";y_{true};p_{T,true} (GeV/c);%s efficiency", w.label));
      m->SetMinimum(-1e-6);
      m->SetMaximum(1.0);
      m->GetZaxis()->SetTitleOffset(1.25);
      if (pair == 0) m->GetYaxis()->SetMoreLogLabels();
      m->Draw("COLZ");
      MaskUndefinedBins(M[w.sp].map);
      DrawBinCounts(M[w.sp].map, pair == 0 ? 0.018 : 0.028);
      DrawHeader(std::string(sampleLabel) + "  " + w.label);
    }
  }
  // page 5/6: efficiency vs pT / y, + and - overlaid
  for (int isPt = 1; isPt >= 0; --isPt)
  {
    auto *c = page(isPt ? "eff_vs_pt" : "eff_vs_y", 3, 1, 1800, 600);
    for (int pair = 0; pair < 3; ++pair)
    {
      c->cd(pair + 1);
      std::vector<Curve> curves;
      for (int k = 0; k < 2; ++k)
      {
        const SpeciesWindow &w = Windows()[2 * pair + k];
        curves.push_back({isPt ? M[w.sp].pt : M[w.sp].y, w.label, w.label, w.color, w.marker, 1, ""});
      }
      DrawEfficiencyCurves(curves, isPt, sampleLabel, true);
      DrawHeader(sampleLabel);
    }
  }
  // page 7: reference comparison vs pT
  if (ref)
  {
    auto *c = page("final_vs_ref_pt", 3, 1, 1800, 600);
    for (int pair = 0; pair < 3; ++pair)
    {
      c->cd(pair + 1);
      std::vector<Curve> curves;
      for (int k = 0; k < 2; ++k)
      {
        const SpeciesWindow &w = Windows()[2 * pair + k];
        curves.push_back({M[w.sp].pt, std::string(w.label) + " final", w.label, w.color, w.marker, 1, ""});
        curves.push_back({R[w.sp].pt, std::string(w.label) + " " + refLabel, std::string(w.label) + "ref", kGray + 2 - k, w.marker, 2, ""});
      }
      DrawEfficiencyCurves(curves, true, sampleLabel, true);
      DrawHeader(sampleLabel);
    }
  }
  const std::string pdf = dir + "/TrkEffFinalMaps.pdf";
  pages.front()->Print((pdf + "[").c_str());
  for (auto *c : pages) c->Print(pdf.c_str(), Form("Title:%s", c->GetName()));
  pages.front()->Print((pdf + "]").c_str());
  for (auto *c : pages)
  {
    c->SaveAs((dir + "/figures/" + c->GetName() + ".pdf").c_str());
    c->SaveAs((dir + "/figures/" + c->GetName() + ".png").c_str());
  }
  Out(Form("%d checks, %d failed", gN, gFail));
  Out(gFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  Out("wrote " + dir + "/TrkEffFinalMaps.root and " + pdf);
  if (gLog) fclose(gLog);
}
