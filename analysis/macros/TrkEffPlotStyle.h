// =============================================================================
// TrkEffPlotStyle.h - common style and helpers for the tracking-efficiency plots
//
// Used by PlotTrkEff.C and LowPtDiagnostics.C (ROOT macros, compiled with
// ACLiC). Only ROOT is needed.
//
// Statistical conventions (identical to TrkEffAnalysis):
//   - efficiencies come from the stored TEfficiency objects
//     (Clopper-Pearson, 68.3% CL);
//   - a bin whose denominator is 0 has an UNDEFINED efficiency. It is never
//     drawn as a measured 0: it is left out of 1D graphs and covered by a
//     hatched grey box in 2D maps.
// =============================================================================
#ifndef TRKEFFPLOTSTYLE_H
#define TRKEFFPLOTSTYLE_H

#include "../src/EffTypes.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace trkeffplot
{
// ---- species ------------------------------------------------------------------
struct SpeciesStyle
{
  const char *tag;    // histogram-name tag in the TrkEffAnalysis output
  const char *label;  // ROOT LaTeX label
  int color;
  int marker;  // filled for positive, open for negative particles
};

// The identified species of EffTypes.h. Particle / antiparticle share a colour;
// positive = filled, negative = open marker.
inline const std::vector<SpeciesStyle> &Species()
{
  static const std::vector<SpeciesStyle> s = []
  {
    std::vector<SpeciesStyle> v;
    for (const auto *e : efftypes::MeasuredSpecies(false)) v.push_back({e->tag, e->label, e->color, e->marker});
    return v;
  }();
  return s;
}

// ---- global style ---------------------------------------------------------------
inline void SetTrkEffStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetTitleSize(0.050, "xyz");
  gStyle->SetTitleOffset(1.10, "x");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.10, "z");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetEndErrorSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(kViridis);  // perceptually uniform, colour-blind friendly
  gStyle->SetNumberContours(50);
}

// "sPHENIX Simulation" header plus a free sample label, top left of the current pad.
inline void DrawHeader(const std::string &sampleLabel, double x = 0.0, double y = 0.0)
{
  const double px = x > 0 ? x : gPad->GetLeftMargin() + 0.02;
  const double py = y > 0 ? y : 1 - gPad->GetTopMargin() + 0.015;
  const std::string text = "#bf{#it{sPHENIX}} Simulation   " + sampleLabel;
  // shrink the text in narrow pads so that the header always fits (approx. 0.5 em per character)
  const double wpx = gPad->GetWw() * gPad->GetAbsWNDC(), hpx = gPad->GetWh() * gPad->GetAbsHNDC();
  const double fit = (1 - px - 0.02) * wpx / (0.56 * (text.size() - 10) * hpx);  // ~0.56 em per character; 10 markup characters
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(std::min(0.040, fit));
  t.SetTextAlign(11);
  t.DrawLatex(px, py, text.c_str());
}

// ---- efficiency helpers -----------------------------------------------------------
// Graph of a 1D TEfficiency; bins with total = 0 are omitted by TEfficiency
// (no fake zero points). Points carry the Clopper-Pearson errors.
inline TGraphAsymmErrors *EfficiencyGraph(const TEfficiency *e, const SpeciesStyle &s)
{
  auto *g = const_cast<TEfficiency *>(e)->CreateGraph();
  g->SetName(Form("g_%s", e->GetName()));
  g->SetMarkerStyle(s.marker);
  g->SetMarkerColor(s.color);
  g->SetLineColor(s.color);
  g->SetMarkerSize(1.3);
  g->SetLineWidth(2);
  return g;
}

// Number of bins of a 1D or 2D TEfficiency (inside the axis range) with total = 0.
inline int CountUndefinedBins(const TEfficiency *e)
{
  const TH1 *t = e->GetTotalHistogram();
  int n = 0;
  const int nx = t->GetNbinsX(), ny = t->GetDimension() > 1 ? t->GetNbinsY() : 1;
  for (int ix = 1; ix <= nx; ++ix)
    for (int iy = 1; iy <= ny; ++iy)
      if (t->GetBinContent(t->GetDimension() > 1 ? t->GetBin(ix, iy) : ix) <= 0) ++n;
  return n;
}

// Map of a 2D efficiency for drawing: content = passed / total where total > 0,
// content = 0 where total = 0 (these bins are covered by MaskUndefinedBins, so
// the 0 is never visible). The returned histogram is a new, detached object.
inline TH2D *EfficiencyMap(const TEfficiency *e, const char *name)
{
  const auto *tot = dynamic_cast<const TH2 *>(e->GetTotalHistogram());
  const auto *pas = dynamic_cast<const TH2 *>(e->GetPassedHistogram());
  auto *m = new TH2D(name, "", tot->GetNbinsX(), tot->GetXaxis()->GetXbins()->GetArray(), tot->GetNbinsY(),
                     tot->GetYaxis()->GetXbins()->GetArray());
  m->SetDirectory(nullptr);
  for (int ix = 1; ix <= tot->GetNbinsX(); ++ix)
    for (int iy = 1; iy <= tot->GetNbinsY(); ++iy)
    {
      const double d = tot->GetBinContent(ix, iy);
      m->SetBinContent(ix, iy, d > 0 ? pas->GetBinContent(ix, iy) / d : 0.);
    }
  return m;
}

// Cover every bin with total = 0 by an opaque light-grey box with grey hatching.
// The boxes are owned by the pad (kCanDelete).
inline void MaskUndefinedBins(const TEfficiency *e)
{
  const auto *tot = dynamic_cast<const TH2 *>(e->GetTotalHistogram());
  for (int ix = 1; ix <= tot->GetNbinsX(); ++ix)
    for (int iy = 1; iy <= tot->GetNbinsY(); ++iy)
    {
      if (tot->GetBinContent(ix, iy) > 0) continue;
      const double x1 = tot->GetXaxis()->GetBinLowEdge(ix), x2 = tot->GetXaxis()->GetBinUpEdge(ix);
      const double y1 = tot->GetYaxis()->GetBinLowEdge(iy), y2 = tot->GetYaxis()->GetBinUpEdge(iy);
      auto *fill = new TBox(x1, y1, x2, y2);
      fill->SetFillColor(kGray);
      fill->SetFillStyle(1001);
      fill->SetBit(kCanDelete);
      fill->Draw();
      auto *hatch = new TBox(x1, y1, x2, y2);
      hatch->SetFillColor(kGray + 2);
      hatch->SetFillStyle(3354);
      hatch->SetBit(kCanDelete);
      hatch->Draw();
    }
}

// Write "passed/total" in every defined bin of a 2D TEfficiency, so that the
// statistics behind each colour are visible (low-statistics bins are not hidden).
// On a logarithmic y axis the text is placed at the geometric bin centre.
inline void DrawBinCounts(const TEfficiency *e, double textSize = 0.020)
{
  const auto *tot = dynamic_cast<const TH2 *>(e->GetTotalHistogram());
  const auto *pas = dynamic_cast<const TH2 *>(e->GetPassedHistogram());
  TLatex t;
  t.SetTextFont(42);
  t.SetTextSize(textSize);
  t.SetTextAlign(22);
  for (int ix = 1; ix <= tot->GetNbinsX(); ++ix)
    for (int iy = 1; iy <= tot->GetNbinsY(); ++iy)
    {
      const double d = tot->GetBinContent(ix, iy);
      if (d <= 0) continue;
      const double eff = pas->GetBinContent(ix, iy) / d;
      t.SetTextColor(eff > 0.6 ? kBlack : kWhite);  // readable on the viridis palette
      const double lo = tot->GetYaxis()->GetBinLowEdge(iy), hi = tot->GetYaxis()->GetBinUpEdge(iy);
      const double yc = (gPad->GetLogy() && lo > 0) ? std::sqrt(lo * hi) : 0.5 * (lo + hi);
      t.DrawLatex(tot->GetXaxis()->GetBinCenter(ix), yc,
                  Form("%.0f/%.0f", pas->GetBinContent(ix, iy), d));
    }
}

// Integrated efficiency of a TEfficiency over all in-range bins, with the
// Clopper-Pearson 68.3% interval of the summed counts.
struct Integrated
{
  double passed = 0, total = 0, eff = 0, low = 0, up = 0;
  bool defined() const { return total > 0; }
};
inline Integrated Integrate(const TEfficiency *e)
{
  const TH1 *t = e->GetTotalHistogram(), *p = e->GetPassedHistogram();
  Integrated r;
  const int nx = t->GetNbinsX(), ny = t->GetDimension() > 1 ? t->GetNbinsY() : 1;
  for (int ix = 1; ix <= nx; ++ix)
    for (int iy = 1; iy <= ny; ++iy)
    {
      const int b = t->GetDimension() > 1 ? t->GetBin(ix, iy) : ix;
      r.total += t->GetBinContent(b);
      r.passed += p->GetBinContent(b);
    }
  if (r.total > 0)
  {
    r.eff = r.passed / r.total;
    r.low = TEfficiency::ClopperPearson(r.total, r.passed, 0.683, false);
    r.up = TEfficiency::ClopperPearson(r.total, r.passed, 0.683, true);
  }
  return r;
}
inline std::string FormatIntegrated(const Integrated &r)
{
  if (!r.defined()) return "undefined (no particles)";
  return Form("%.0f/%.0f = %.3f^{+%.3f}_{-%.3f}", r.passed, r.total, r.eff, r.up - r.eff, r.eff - r.low);
}


// ---- generic 1D efficiency comparison ---------------------------------------------
// One curve = one TEfficiency with its style and labels.
struct Curve
{
  const TEfficiency *eff;
  std::string label;       // legend label
  std::string shortLabel;  // prefix of the statistics line in compact mode
  int color, marker, lineStyle = 1;
  std::string graphName;   // name of the drawn graph (default "g_" + efficiency name)
};

// Draws the curves versus pT (log x) or y into the current pad. The legend sits in
// an enlarged right margin so it never covers data. Every curve shows its counts
// passed/total, the integrated efficiency (Clopper-Pearson 68.3%) and its number
// of empty bins. compact = true (narrow pads): the legend holds only the labels
// and the statistics are written as short lines below it.
// Returns the drawn graphs (owned by the pad).
inline std::vector<TGraphAsymmErrors *> DrawEfficiencyCurves(const std::vector<Curve> &curves, bool isPt,
                                                             const std::string &sample, bool compact = false)
{
  static int frameCounter = 0;
  gPad->SetLogx(isPt);
  gPad->SetRightMargin(compact ? 0.36 : 0.42);
  const TH1 *axisFrom = curves.front().eff->GetTotalHistogram();
  const double lo = axisFrom->GetXaxis()->GetXmin(), hi = axisFrom->GetXaxis()->GetXmax();
  auto *frame = new TH1D(Form("frame_%d", ++frameCounter), "", 100, lo, hi);
  frame->SetDirectory(nullptr);
  frame->SetBit(kCanDelete);
  frame->SetMinimum(0.0);
  frame->SetMaximum(1.25);
  frame->GetXaxis()->SetTitle(isPt ? "p_{T,true} (GeV/c)" : "y_{true}");
  frame->GetYaxis()->SetTitle("tracking efficiency #varepsilon");
  if (isPt)
  {
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();
  }
  frame->Draw("AXIS");
  TLine one;
  one.SetLineStyle(2);
  one.SetLineColor(kGray + 1);
  one.DrawLine(lo, 1.0, hi, 1.0);

  const double legX = 1 - gPad->GetRightMargin() + 0.015;
  const double legBottom = 0.88 - (compact ? 0.075 : 0.115) * curves.size();
  auto *leg = new TLegend(legX, legBottom, 0.995, 0.88);
  leg->SetBit(kCanDelete);
  leg->SetTextSize(curves.size() > 2 ? 0.036 : 0.040);
  leg->SetMargin(0.12);
  leg->SetFillStyle(0);
  std::vector<std::string> statLines;
  std::vector<TGraphAsymmErrors *> graphs;
  for (const auto &c : curves)
  {
    SpeciesStyle st{"", "", c.color, c.marker};
    TGraphAsymmErrors *g = EfficiencyGraph(c.eff, st);
    if (!c.graphName.empty()) g->SetName(c.graphName.c_str());
    g->SetLineStyle(c.lineStyle);
    g->SetBit(kCanDelete);
    g->Draw("P SAME");
    graphs.push_back(g);
    const Integrated I = Integrate(c.eff);
    const int nEmpty = CountUndefinedBins(c.eff);
    const std::string empty = Form("%d empty bin%s (no point)", nEmpty, nEmpty == 1 ? "" : "s");
    if (compact)
    {
      leg->AddEntry(g, c.label.c_str(), "pe");
      statLines.push_back(c.shortLabel + ": " + FormatIntegrated(I));
      statLines.push_back("      " + empty);
    }
    else
    {
      leg->AddEntry(g, (c.label + ": " + FormatIntegrated(I)).c_str(), "pe");
      leg->AddEntry(static_cast<TObject *>(nullptr), ("   " + empty).c_str(), "");
    }
  }
  leg->Draw();
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.031);
  double y = legBottom - 0.06;  // (TLegend NDC coordinates are only set when painted)
  for (const auto &l : statLines)
  {
    t.DrawLatex(legX, y, l.c_str());
    y -= 0.058;
  }
  DrawHeader(sample);
  return graphs;
}

}  // namespace trkeffplot

#endif
