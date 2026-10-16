// =============================================================================
// PlotVtxEff.C
//
// Makes all plots and efficiency tables from the output of bin/VtxEffAnalysis
// (or from an hadd of several such outputs).
//
// Usage (from VertexEfficiencyStudy/):
//   root -l -b -q 'macros/PlotVtxEff.C+("output/VtxEff_validation.root","plots","summary/VtxEff_tables.txt")'
//
// Produces in <plotDir>:
//   VertexEfficiency.pdf, MultiplicityResponse.pdf, TrackQA.pdf,
//   TruthMatchingQA.pdf and the combined VertexEfficiency_CompleteQA.pdf.
// Writes efficiency tables to <tableFile>.
//
// Statistics: all efficiencies use TEfficiency with Clopper-Pearson ("exact")
// binomial intervals at 68.3% CL. Clopper-Pearson never gives a zero-width
// interval: a bin with 0 vertexed out of N_all events gets [0, upper limit],
// and the tables also give the one-sided 95% CL upper limit. Efficiencies are
// only computed as passed/total from raw counts; nothing is divided by an
// efficiency.
//
// Layout rules used on every page (so legends never cover data):
//   - the sample description is written in the top margin, above the frame;
//   - each frame gets explicit head-room (linear: x1.6, log: extra decades,
//     efficiency: y up to 1.4) and legends/text are placed in that empty band;
//   - overlaid efficiency series are shifted slightly in x so that points
//     with identical x do not hide each other.
// =============================================================================

#include <TCanvas.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
TFile *gIn = nullptr;
double gNevents = 0;
const double kCL = 0.683;  // 1-sigma equivalent confidence level

// Sample description shown on every plot.
const char *kSample1 = "PYTHIA8 p+p #sqrt{s} = 200 GeV MB (Detroit), GEANT4";
const char *kSample2 = "pileup-free sample, polyseed tracking";

// Colour palette (colour-blind friendly, also readable in grey scale).
int cBlue = kBlue, cRed = kRed, cGreen = kGreen + 2, cOrange = kOrange + 1, cGray = kGray + 1;
const int cBlack = kBlack;
int cLightBlue = kAzure - 9;
void InitColors()
{
  cBlue = TColor::GetColor("#1f5fa8");
  cRed = TColor::GetColor("#c0392b");
  cGreen = TColor::GetColor("#2e8b57");
  cOrange = TColor::GetColor("#e08a1e");
  cGray = TColor::GetColor("#95a5a6");
  cLightBlue = TColor::GetColor("#a9c6e8");
}

template <class T>
T *Get(const char *name)
{
  T *o = dynamic_cast<T *>(gIn->Get(name));
  if (!o)
  {
    std::cerr << "FATAL: object " << name << " not found" << std::endl;
    gSystem->Exit(1);
  }
  return o;
}

// -----------------------------------------------------------------------------
// Style
// -----------------------------------------------------------------------------
void SetStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetTitleSize(0.043, "XYZ");
  gStyle->SetTitleOffset(1.15, "X");
  gStyle->SetTitleOffset(1.40, "Y");
  gStyle->SetTitleOffset(1.05, "Z");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetEndErrorSize(0);
  gStyle->SetErrorX(0.);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(99);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetHistLineWidth(2);
  // PDF page size = canvas aspect ratio (800 x 650).
  gStyle->SetPaperSize(20., 16.25);
  // Histograms read from the analysis file carry the style they were created
  // with; force the current style so that all pages use identical fonts.
  gROOT->ForceStyle();
}

// Page title (bold) and the two-line sample description in the top margin.
void Label(const char *title)
{
  if (gPad && gPad->GetCanvas()) gPad->GetCanvas()->SetTitle(title);  // used by the figure manifest
  TLatex l;
  l.SetNDC();
  l.SetTextFont(62);
  l.SetTextSize(0.036);
  l.DrawLatex(0.14, 0.945, title);
  l.SetTextFont(42);
  l.SetTextSize(0.027);
  l.SetTextColor(kGray + 2);
  l.DrawLatex(0.14, 0.905, kSample1);
  l.DrawLatex(0.14, 0.872, Form("%s;  %.0f collisions", kSample2, gNevents));
}

// Text lines inside the frame (top-left of the head-room band by default).
void Info(const std::vector<std::string> &lines, double x = 0.18, double y = 0.80, double size = 0.031)
{
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(size);
  for (size_t i = 0; i < lines.size(); ++i) l.DrawLatex(x, y - i * 1.45 * size, lines[i].c_str());
}

TLegend *Leg(double x1, double y1, double x2, double y2, double size = 0.033)
{
  auto *lg = new TLegend(x1, y1, x2, y2);
  lg->SetTextFont(42);
  lg->SetTextSize(size);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  return lg;
}

// Global marker-size scale applied to every histogram/graph marker.
const double kMarkerScale = 1.3;

void Style(TH1 *h, int color, int marker = 20, double msize = 1.0)
{
  msize *= kMarkerScale;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(msize);
  h->SetLineWidth(2);
}

void Style(TGraph *g, int color, int marker = 20, double msize = 1.1)
{
  msize *= kMarkerScale;
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(msize);
  g->SetLineWidth(2);
}

// Set a common y range with free head-room for legends.
//   linear: max * linFactor
//   log   : the data occupy the lower (1 - freeFrac) of the decades
void Headroom(const std::vector<TH1 *> &hs, bool logy, double linFactor = 1.6, double freeFrac = 0.33,
              double logMin = 0.5)
{
  double mx = 0;
  for (TH1 *h : hs)
    for (int b = h->GetXaxis()->GetFirst(); b <= h->GetXaxis()->GetLast(); ++b)
      mx = std::max(mx, h->GetBinContent(b) + h->GetBinError(b));
  if (mx <= 0) mx = 1;
  for (TH1 *h : hs)
  {
    if (logy)
    {
      const double dec = std::log10(mx / logMin);
      h->SetMinimum(logMin);
      h->SetMaximum(logMin * std::pow(10., dec / (1. - freeFrac)));
    }
    else
    {
      h->SetMinimum(0.);
      h->SetMaximum(mx * linFactor);
    }
  }
}

// passed / total -> TEfficiency (Clopper-Pearson, 68.3% CL)
TEfficiency *Eff(const TH1 *pass, const TH1 *tot, const char *name)
{
  if (!TEfficiency::CheckConsistency(*pass, *tot))
  {
    std::cerr << "FATAL: inconsistent histograms for " << name << std::endl;
    gSystem->Exit(1);
  }
  auto *e = new TEfficiency(*pass, *tot);
  e->SetName(name);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(kCL);
  return e;
}

// Graph of an efficiency, for bins with entries in the denominator.
//   unitBins: no horizontal error bar for bins of width 1 (single integers)
//   shift   : x offset to separate overlaid series
TGraphAsymmErrors *EffGraph(TEfficiency *e, double xmax = 1e9, double shift = 0., bool unitBins = true)
{
  auto *g = new TGraphAsymmErrors();
  const TH1 *tot = e->GetTotalHistogram();
  int ip = 0;
  for (int b = 1; b <= tot->GetNbinsX(); ++b)
  {
    if (tot->GetBinContent(b) <= 0) continue;
    const double x = tot->GetBinCenter(b);
    if (x > xmax) continue;
    const double w = tot->GetBinWidth(b);
    const double ex = (unitBins && w <= 1.0 + 1e-9) ? 0. : w / 2.;
    g->SetPoint(ip, x + shift, e->GetEfficiency(b));
    g->SetPointError(ip, ex, ex, e->GetEfficiencyErrorLow(b), e->GetEfficiencyErrorUp(b));
    ++ip;
  }
  return g;
}

TH1F *Frame(TCanvas *c, double x0, double x1, double y0, double y1, const char *xt, const char *yt)
{
  c->cd();
  TH1F *f = c->DrawFrame(x0, y0, x1, y1);
  f->GetXaxis()->SetTitle(xt);
  f->GetYaxis()->SetTitle(yt);
  f->GetXaxis()->SetNdivisions(510);
  f->GetYaxis()->SetNdivisions(507);
  return f;
}

// Efficiency frame: y from 0 to 1.4; data live below 1, the band above ~1.05
// is reserved for legends/text. A dotted line marks the physical limit 1.
TH1F *EffFrame(TCanvas *c, double x0, double x1, const char *xt, const char *yt)
{
  TH1F *f = Frame(c, x0, x1, 0., 1.4, xt, yt);
  TLine one(x0, 1., x1, 1.);
  one.SetLineStyle(3);
  one.SetLineColor(kGray + 2);
  one.DrawClone();
  return f;
}

// Collect canvases per section; printed at the end.
struct Book
{
  std::string file;
  std::vector<TCanvas *> pages;
};

TCanvas *NewPage(Book &b, const char *name, bool colz = false)
{
  auto *c = new TCanvas(name, name, 800, 650);
  if (colz) c->SetRightMargin(0.17);
  b.pages.push_back(c);
  return c;
}

std::ofstream gManifest;  // figure manifest: pdf, page, canvas name, page title

void PrintBook(const Book &b, const std::string &dir, bool manifest = true)
{
  const std::string f = dir + "/" + b.file;
  if (manifest && gManifest.is_open())
    for (size_t i = 0; i < b.pages.size(); ++i)
      gManifest << b.file << "\t" << i + 1 << "\t" << b.pages[i]->GetName() << "\t" << b.pages[i]->GetTitle() << "\n";
  for (size_t i = 0; i < b.pages.size(); ++i)
  {
    std::string opt = f;
    if (b.pages.size() > 1)
    {
      if (i == 0) opt = f + "(";
      else if (i + 1 == b.pages.size()) opt = f + ")";
    }
    b.pages[i]->Print(opt.c_str(), Form("pdf Title:%s", b.pages[i]->GetName()));
  }
}

// P(n|N) = N_events(N, n) / N_events(N): each true-N column normalised to 1.
TH2D *ColumnNormalised(const TH2D *h, const char *name)
{
  auto *p = (TH2D *) h->Clone(name);
  for (int ix = 1; ix <= p->GetNbinsX(); ++ix)
  {
    double s = 0;
    for (int iy = 0; iy <= p->GetNbinsY() + 1; ++iy) s += p->GetBinContent(ix, iy);
    for (int iy = 0; iy <= p->GetNbinsY() + 1; ++iy)
    {
      p->SetBinContent(ix, iy, s > 0 ? p->GetBinContent(ix, iy) / s : 0.);
      p->SetBinError(ix, iy, 0.);
    }
  }
  return p;
}

// Fraction of each y category per x bin (for composition-vs-n plots).
std::vector<TH1D *> CategoryFractions(const TH2D *h, const char *base, int xmaxBin)
{
  std::vector<TH1D *> out;
  const int ny = h->GetNbinsY();
  for (int iy = 1; iy <= ny; ++iy)
  {
    auto *f = new TH1D(Form("%s_%d", base, iy), "", xmaxBin, h->GetXaxis()->GetBinLowEdge(1),
                       h->GetXaxis()->GetBinUpEdge(xmaxBin));
    f->SetDirectory(nullptr);
    out.push_back(f);
  }
  for (int ix = 1; ix <= xmaxBin; ++ix)
  {
    double s = 0;
    for (int iy = 1; iy <= ny; ++iy) s += h->GetBinContent(ix, iy);
    if (s <= 0) continue;
    for (int iy = 1; iy <= ny; ++iy) out[iy - 1]->SetBinContent(ix, h->GetBinContent(ix, iy) / s);
  }
  return out;
}

// Unit-area copy (for shape comparisons).
TH1D *UnitArea(const TH1 *h, const char *name)
{
  auto *c = (TH1D *) h->Clone(name);
  c->SetDirectory(nullptr);
  const double s = c->Integral(0, c->GetNbinsX() + 1);
  if (s > 0) c->Scale(1. / s);
  return c;
}

// Write an efficiency table: one row per bin with entries.
// continuous = true: label bins by their edges [lo, hi) (e.g. z in cm)
void WriteEffTable(std::ofstream &t, TEfficiency *e, const char *var, bool wide, bool continuous = false)
{
  const TH1 *tot = e->GetTotalHistogram();
  const TH1 *pas = e->GetPassedHistogram();
  t << Form("  %-10s %10s %10s %9s %9s %9s %12s  %s\n", var, "N_all", "N_vtx", "eps", "CP68 low", "CP68 up",
            "UL95 (if 0)", "note");
  for (int b = 1; b <= tot->GetNbinsX(); ++b)
  {
    const double nt = tot->GetBinContent(b);
    if (nt <= 0) continue;
    const double np = pas->GetBinContent(b);
    const double eff = e->GetEfficiency(b);
    const double lo = eff - e->GetEfficiencyErrorLow(b);
    const double up = eff + e->GetEfficiencyErrorUp(b);
    std::string lab;
    const double l = tot->GetXaxis()->GetBinLowEdge(b) + 0.5;
    const double h = tot->GetXaxis()->GetBinUpEdge(b) - 0.5;
    if (continuous) lab = Form("[%.0f,%.0f)", tot->GetXaxis()->GetBinLowEdge(b), tot->GetXaxis()->GetBinUpEdge(b));
    else if (wide && h > l) lab = Form("%.0f-%.0f", l, h);
    else lab = Form("%.0f", l);
    std::string ul = "";
    // One-sided 95% CL upper limit for 0 successes in nt trials: 1 - 0.05^(1/nt).
    // TEfficiency::ClopperPearson(level) is two-sided (alpha = (1-level)/2 per side),
    // so the one-sided 95% limit needs level = 0.90. (Before 2026-10-12 level 0.95
    // was used here by mistake, which gives the one-sided 97.5% limit.)
    if (np == 0) ul = Form("%.4f", TEfficiency::ClopperPearson(nt, 0, 0.90, true));
    t << Form("  %-10s %10.0f %10.0f %9.4f %9.4f %9.4f %12s  %s\n", lab.c_str(), nt, np, eff, lo, up, ul.c_str(),
              nt < 10 ? "low stats (N_all < 10)" : "");
  }
}

// Bar charts: drop the per-bin errors so that ROOT draws plain bars
// (with Sumw2 set, ROOT would otherwise add error ticks on top of the bars).
void NoErrors(TH1 *h) { h->Sumw2(false); }

// Rebin factor for sparse residual plots, so that shapes stay readable.
int RebinFor(double entries)
{
  if (entries < 1000) return 5;
  if (entries < 10000) return 2;
  return 1;
}
}  // namespace

void PlotVtxEff(const char *inFile = "output/VtxEff_validation.root", const char *plotDir = "plots",
                const char *tableFile = "summary/VtxEff_tables.txt")
{
  SetStyle();
  InitColors();
  gIn = TFile::Open(inFile, "READ");
  if (!gIn || gIn->IsZombie())
  {
    std::cerr << "FATAL: cannot open " << inFile << std::endl;
    return;
  }
  gSystem->mkdir(plotDir, true);

  TH1D *hStats = Get<TH1D>("hStats");
  auto stat = [&](const char *lab) { return hStats->GetBinContent(hStats->GetXaxis()->FindBin(lab)); };
  gNevents = stat("events");

  // ---------------------------------------------------------------------------
  // Efficiencies (fine and wide bins; same binning for numerator/denominator)
  // ---------------------------------------------------------------------------
  TEfficiency *en_fine = Eff(Get<TH1D>("hn_vtx"), Get<TH1D>("hn_all"), "eff_n_fine");
  TEfficiency *eN_fine = Eff(Get<TH1D>("hN_vtx"), Get<TH1D>("hN_all"), "eff_N_fine");
  TEfficiency *en_wide = Eff(Get<TH1D>("hnW_vtx"), Get<TH1D>("hnW_all"), "eff_n_wide");
  TEfficiency *eN_wide = Eff(Get<TH1D>("hNW_vtx"), Get<TH1D>("hNW_all"), "eff_N_wide");
  TEfficiency *eZ = Eff(Get<TH1D>("hZtrue_vtx"), Get<TH1D>("hZtrue_all"), "eff_ztrue");
  // central true-z window (the vertex-finder acceptance falls to zero at |z| ~ 20 cm)
  TEfficiency *en_fine_zc = Eff(Get<TH1D>("hn_vtx_zc"), Get<TH1D>("hn_all_zc"), "eff_n_fine_zc");
  TEfficiency *eN_fine_zc = Eff(Get<TH1D>("hN_vtx_zc"), Get<TH1D>("hN_all_zc"), "eff_N_fine_zc");
  TEfficiency *en_wide_zc = Eff(Get<TH1D>("hnW_vtx_zc"), Get<TH1D>("hnW_all_zc"), "eff_n_wide_zc");
  TEfficiency *eN_wide_zc = Eff(Get<TH1D>("hNW_vtx_zc"), Get<TH1D>("hNW_all_zc"), "eff_N_wide_zc");
  const std::string zcTitle = Get<TH1D>("hn_all_zc")->GetTitle();  // "All collisions, |z_{true}| < X cm;..."
  std::string zcText = zcTitle.substr(zcTitle.find('|'));
  zcText = zcText.substr(0, zcText.find(';'));
  double zcut = 10.;
  {
    const auto p = zcText.find("< ");
    if (p != std::string::npos) zcut = std::atof(zcText.c_str() + p + 2);
  }
  // Optional MVTX requirement (config key good_nmvtx_min; " nmvtx>=K" in the config string).
  std::string mvtxText;
  if (auto *cfgObj = dynamic_cast<TNamed *>(gIn->Get("config")))
  {
    const std::string c = cfgObj->GetTitle();
    const auto p = c.find("nmvtx>=");
    if (p != std::string::npos) mvtxText = Form(", N_{MVTX clus} #geq %d", std::stoi(c.substr(p + 7)));
  }
  const std::string cutText = "good track: p_{T} > 0.1 GeV/c, |#eta| < 1.1, N_{TPC clus} #geq 18" + mvtxText;
  const std::string vtxText = "vertex: #geq 1 reconstructed tracking vertex";
  const std::string clText = "error bars: Clopper-Pearson 68.3% CL";
  const char *xn = "reconstructed good-track multiplicity n";
  const char *xN = "true primary charged multiplicity N";

  // ===========================================================================
  // 1. VertexEfficiency.pdf
  // ===========================================================================
  Book bEff{"VertexEfficiency.pdf", {}};

  // Upper x limit of a wide-bin efficiency plot: just past the last populated bin.
  auto wideMax = [](TEfficiency *e) {
    const TH1 *tot = e->GetTotalHistogram();
    double up = 10.5;
    for (int b = 1; b <= tot->GetNbinsX(); ++b)
      if (tot->GetBinContent(b) > 0) up = tot->GetXaxis()->GetBinUpEdge(b);
    return std::min(tot->GetXaxis()->GetXmax(), up + 5.);
  };

  // Efficiency page with one or two series (second = central z or true N).
  auto effPage = [&](const char *cname, const char *title, double xmax, const char *xt, const char *yt,
                     TEfficiency *e1, int c1, int m1, const char *l1, TEfficiency *e2, int c2, int m2, const char *l2,
                     bool unitBins, const std::vector<std::string> &info) {
    TCanvas *c = NewPage(bEff, cname);
    EffFrame(c, -0.5, xmax, xt, yt);
    const double sh = (e2 && unitBins) ? 0.15 : 0.;
    auto *g1 = EffGraph(e1, xmax, -sh, unitBins);
    Style(g1, c1, m1);
    g1->Draw("P same");
    TGraphAsymmErrors *g2 = nullptr;
    if (e2)
    {
      g2 = EffGraph(e2, xmax, +sh, unitBins);
      Style(g2, c2, m2);
      g2->Draw("P same");
    }
    // legend in the upper band (right), info text in the upper band (left)
    auto *lg = Leg(0.56, e2 ? 0.70 : 0.745, 0.93, 0.82, 0.032);
    lg->AddEntry(g1, l1, "pe");
    if (g2) lg->AddEntry(g2, l2, "pe");
    lg->Draw();
    Info(info, 0.18, 0.80, 0.027);
    Label(title);
  };
  const std::vector<std::string> infoN = {"good track:", "  p_{T} > 0.1 GeV/c, |#eta| < 1.1,", "  N_{TPC clus} #geq 18" + mvtxText};
  const std::vector<std::string> infoTrueN = {"N: charged G4 primaries,", "  p_{T} > 0.1 GeV/c, |#eta| < 1.1"};
  effPage("eff_n_fine", "Vertex-finding efficiency vs reconstructed multiplicity", 30.5, xn, "#varepsilon_{vtx}(n)",
          en_fine, cBlue, 20, "all true z", nullptr, 0, 0, "", true, infoN);
  effPage("eff_n_zc", "#varepsilon_{vtx}(n): all z vs central z", 30.5, xn, "#varepsilon_{vtx}(n)", en_fine, cBlue, 24,
          "all true z", en_fine_zc, cBlack, 20, zcText.c_str(), true, infoN);
  effPage("eff_n_wide_zc", "#varepsilon_{vtx}(n), wide bins", wideMax(en_wide), xn, "#varepsilon_{vtx}(n)", en_wide, cBlue, 25,
          "all true z", en_wide_zc, cBlack, 21, zcText.c_str(), false, infoN);
  effPage("eff_N_fine", "Vertex-finding efficiency vs true multiplicity", 30.5, xN, "#varepsilon_{vtx}(N)", eN_fine,
          cRed, 20, "all true z", nullptr, 0, 0, "", true, infoTrueN);
  effPage("eff_N_zc", "#varepsilon_{vtx}(N): all z vs central z", 30.5, xN, "#varepsilon_{vtx}(N)", eN_fine, cRed, 24,
          "all true z", eN_fine_zc, cBlack, 20, zcText.c_str(), true, infoTrueN);
  effPage("eff_N_wide_zc", "#varepsilon_{vtx}(N), wide bins", wideMax(eN_wide), xN, "#varepsilon_{vtx}(N)", eN_wide, cRed, 25,
          "all true z", eN_wide_zc, cBlack, 21, zcText.c_str(), false, infoTrueN);
  effPage("eff_n_vs_N", "#varepsilon_{vtx}: reconstructed n vs true N (all z)", 30.5, "multiplicity (n or N)",
          "#varepsilon_{vtx}", en_fine, cBlue, 20, "vs reconstructed n", eN_fine, cRed, 24, "vs true N", true,
          {"all true z", "", ""});

  // --- epsilon vs true z
  {
    TCanvas *c = NewPage(bEff, "eff_ztrue");
    const double zm = Get<TH1D>("hZtrue_all")->GetXaxis()->GetXmax();
    Frame(c, -zm, zm, 0., 1.4, "true collision z (cm)", "#varepsilon_{vtx} (all multiplicities)");
    TLine one(-zm, 1., zm, 1.);
    one.SetLineStyle(3);
    one.SetLineColor(kGray + 2);
    one.DrawClone();
    // lines owned by the canvas (heap), because the legend refers to them
    auto *lz = new TLine(-zcut, 0., -zcut, 1.);
    lz->SetLineStyle(2);
    lz->SetLineWidth(2);
    lz->SetLineColor(cGreen);
    lz->Draw();
    auto *lz2 = (TLine *) lz->Clone();
    lz2->SetX1(zcut);
    lz2->SetX2(zcut);
    lz2->Draw();
    auto *g = EffGraph(eZ, 1e9, 0., false);
    Style(g, cBlack, 20);
    g->Draw("P same");
    auto *lg = Leg(0.18, 0.70, 0.70, 0.82, 0.032);
    lg->AddEntry(g, "#varepsilon_{vtx} in 2 cm bins of true z", "pe");
    lg->AddEntry(lz, Form("central window |z_{true}| < %.0f cm", zcut), "l");
    lg->Draw();
    Label("Vertex-finding efficiency vs true collision z");
  }
  // --- multiplicity distributions, vertexed vs vertexless
  auto threeMult = [&](const char *cname, TH1 *a0, TH1 *v0, TH1 *nv0, const char *xt, const char *title) {
    TCanvas *c = NewPage(bEff, cname);
    c->SetLogy();
    auto *a = (TH1D *) a0->Clone(Form("%s_a", cname));
    auto *v = (TH1D *) v0->Clone(Form("%s_v", cname));
    auto *nv = (TH1D *) nv0->Clone(Form("%s_nv", cname));
    Style(a, cBlack, 20, 0.9);
    Style(v, cBlue, 21, 0.85);
    Style(nv, cRed, 24, 0.9);
    for (TH1 *h : {(TH1 *) a, (TH1 *) v, (TH1 *) nv}) h->GetXaxis()->SetRange(1, 51);
    Headroom({a, v, nv}, true, 1.6, 0.36);
    a->GetXaxis()->SetTitle(xt);
    a->GetYaxis()->SetTitle("collisions");
    a->Draw("E1");
    v->Draw("E1 same");
    nv->Draw("E1 same");
    auto *lg = Leg(0.52, 0.67, 0.93, 0.82, 0.031);
    lg->AddEntry(a, Form("all collisions (mean %.2f)", a->GetMean()), "pe");
    lg->AddEntry(v, Form("vertexed (mean %.2f)", v->GetMean()), "pe");
    lg->AddEntry(nv, Form("vertexless (mean %.2f)", nv->GetMean()), "pe");
    lg->Draw();
    Label(title);
  };
  {
    auto *Nnv = (TH1D *) Get<TH1D>("hN_all")->Clone("hN_novtx_tmp");
    Nnv->Add(Get<TH1D>("hN_vtx"), -1);
    threeMult("dist_n", Get<TH1D>("hn_all"), Get<TH1D>("hn_vtx"), Get<TH1D>("hn_novtx"), xn,
              "Reconstructed n: vertexed vs vertexless collisions");
    threeMult("dist_N", Get<TH1D>("hN_all"), Get<TH1D>("hN_vtx"), Nnv, xN, "True N: vertexed vs vertexless collisions");
  }
  // --- true z, vertexed vs vertexless (unit area)
  {
    TCanvas *c = NewPage(bEff, "ztrue_shapes");
    auto *v = UnitArea(Get<TH1D>("hZtrue_vtx"), "zt_v");
    auto *nv = UnitArea(Get<TH1D>("hZtrue_novtx"), "zt_nv");
    Style(v, cBlue, 21, 0.85);
    Style(nv, cRed, 24, 0.9);
    Headroom({v, nv}, false, 1.65);
    v->GetXaxis()->SetTitle("true collision z (cm)");
    v->GetYaxis()->SetTitle("fraction of collisions / 2 cm");
    v->Draw("E1");
    nv->Draw("E1 same");
    auto *lg = Leg(0.18, 0.72, 0.55, 0.82, 0.032);
    lg->AddEntry(v, "vertexed collisions", "pe");
    lg->AddEntry(nv, "vertexless collisions", "pe");
    lg->Draw();
    Label("True collision z: vertexed vs vertexless (unit area)");
  }
  // --- vertex residuals
  {
    const char *names[3] = {"hDz", "hDx", "hDy"};
    const char *coord[3] = {"z", "x", "y"};
    for (int i = 0; i < 3; ++i)
    {
      TCanvas *c = NewPage(bEff, Form("res_%s", coord[i]));
      auto *h = (TH1D *) Get<TH1D>(names[i])->Clone(Form("%s_p", names[i]));
      h->SetDirectory(nullptr);
      const double mean = h->GetMean(), rms = h->GetRMS(), ent = h->GetEntries();
      h->Rebin(RebinFor(ent));
      h->SetLineColor(cBlue);
      h->SetFillColor(cLightBlue);
      h->SetMarkerStyle(1);
      Headroom({h}, false, 1.6);
      h->GetXaxis()->SetTitle(Form("%s_{reco} - %s_{true} (cm)", coord[i], coord[i]));
      h->GetYaxis()->SetTitle(Form("vertexed collisions / %.2f mm", 10. * h->GetBinWidth(1)));
      h->Draw("hist");
      h->Draw("E1 same");
      TLine l0(0., 0., 0., h->GetMaximum() / 1.6 * 1.05);
      l0.SetLineStyle(2);
      l0.SetLineColor(kGray + 2);
      l0.DrawClone();
      Info({Form("mean = %+.3f mm", 10. * mean), Form("RMS = %.3f mm", 10. * rms), Form("entries = %.0f", ent)},
           0.18, 0.79, 0.032);
      Label(Form("Vertex %s residual (reconstructed - true)", coord[i]));
    }
  }
  // --- vertices per collision and tracks per vertex
  {
    TCanvas *c = NewPage(bEff, "nvtx");
    c->SetLogy();
    auto *h = (TH1D *) Get<TH1D>("hNvtx")->Clone("nvtx_p");
    h->SetLineColor(cBlue);
    h->SetFillColor(cLightBlue);
    Headroom({h}, true, 1.6, 0.36);
    h->GetYaxis()->SetTitle("collisions");
    h->GetXaxis()->SetRangeUser(-0.5, 3.5);
    h->GetXaxis()->SetNdivisions(4);
    h->Draw("hist");
    Info({Form("0 vertices: %.0f", h->GetBinContent(1)), Form("1 vertex: %.0f", h->GetBinContent(2)),
          Form("> 1 vertex: %.0f", stat("nVtx>1"))},
         0.62, 0.79, 0.032);
    Label("Reconstructed tracking vertices per collision");
  }
  {
    TCanvas *c = NewPage(bEff, "vtx_ntracks");
    auto *h = (TH1D *) Get<TH1D>("hVtxNtracks")->Clone("vtxntrk_p");
    int last = 1;
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetBinContent(b) > 0) last = b;
    h->GetXaxis()->SetRangeUser(-0.5, std::max(10.5, h->GetBinCenter(last) + 5.5));
    h->SetLineColor(cBlue);
    h->SetFillColor(cLightBlue);
    Headroom({h}, false, 1.5);
    h->GetYaxis()->SetTitle("vertexed collisions");
    h->Draw("hist");
    h->Draw("E1 same");
    Info({Form("mean = %.2f tracks", h->GetMean())}, 0.62, 0.79, 0.032);
    Label("Tracks attached to the reconstructed vertex");
  }

  // ===========================================================================
  // 2. MultiplicityResponse.pdf
  // ===========================================================================
  Book bResp{"MultiplicityResponse.pdf", {}};
  TH2D *hR = Get<TH2D>("hResp_all");
  TH2D *hRv = Get<TH2D>("hResp_vtx");
  {
    TCanvas *c = NewPage(bResp, "dist_n_N");
    c->SetLogy();
    auto *n = (TH1D *) Get<TH1D>("hn_all")->Clone("n_dist");
    auto *N = (TH1D *) Get<TH1D>("hN_all")->Clone("N_dist");
    Style(n, cBlue, 20, 0.9);
    Style(N, cRed, 24, 0.9);
    for (TH1 *h : {(TH1 *) n, (TH1 *) N}) h->GetXaxis()->SetRange(1, 51);
    Headroom({n, N}, true, 1.6, 0.36);
    n->GetXaxis()->SetTitle("multiplicity");
    n->GetYaxis()->SetTitle("collisions");
    n->Draw("E1");
    N->Draw("E1 same");
    auto *lg = Leg(0.48, 0.715, 0.93, 0.82, 0.031);
    lg->AddEntry(n, Form("reconstructed n (mean %.2f)", n->GetMean()), "pe");
    lg->AddEntry(N, Form("true N (mean %.2f)", N->GetMean()), "pe");
    lg->Draw();
    Label("Reconstructed n and true N (all collisions)");
  }
  auto colzPage = [&](const char *cname, TH2D *h2, const char *ztitle, const char *title, bool norm,
                      const char *ytitle = "reconstructed n") {
    TCanvas *c = NewPage(bResp, cname, true);
    c->SetLogz();
    auto *h = norm ? ColumnNormalised(h2, Form("%s_h", cname)) : (TH2D *) h2->Clone(Form("%s_h", cname));
    h->GetXaxis()->SetRangeUser(-0.5, 30.5);
    h->GetYaxis()->SetRangeUser(-0.5, 50.5);
    h->GetXaxis()->SetTitle("true N");
    h->GetYaxis()->SetTitle(ytitle);
    h->GetZaxis()->SetTitle(ztitle);
    h->SetMinimum(norm ? 1e-3 : 0.99);
    h->GetZaxis()->SetMoreLogLabels();
    h->Draw("colz");
    auto *d = new TLine(-0.5, -0.5, 30.5, 30.5);
    d->SetLineStyle(2);
    d->SetLineColor(cRed);
    d->SetLineWidth(2);
    d->Draw();
    // legend below the diagonal (lower right is empty: n >= ~N always holds roughly)
    auto *lg = Leg(0.58, 0.17, 0.83, 0.23, 0.030);
    lg->SetFillStyle(1001);
    lg->SetFillColor(kWhite);
    lg->AddEntry(d, "n = N", "l");
    lg->Draw();
    Label(title);
  };
  colzPage("resp_counts", hR, "collisions", "Response matrix: collisions per (N, n), all collisions", false);
  colzPage("resp_norm", hR, "P(n | N)", "Response P(n|N): each true-N column normalised to 1", true);
  colzPage("resp_vtx", hRv, "vertexed collisions", "Response matrix: vertexed collisions only", false);
  {
    TCanvas *c = NewPage(bResp, "mean_n_vs_N");
    auto *pa = hR->ProfileX("prof_all");
    auto *pv = hRv->ProfileX("prof_vtx");
    auto *gpa = new TGraphAsymmErrors();
    auto *gpv = new TGraphAsymmErrors();
    double mx = 0;
    for (int b = 1; b <= pa->FindBin(25.); ++b)
    {
      if (pa->GetBinEntries(b) > 1)
      {
        const int ip = gpa->GetN();
        gpa->SetPoint(ip, pa->GetBinCenter(b) - 0.12, pa->GetBinContent(b));
        gpa->SetPointError(ip, 0, 0, pa->GetBinError(b), pa->GetBinError(b));
        mx = std::max(mx, pa->GetBinContent(b) + pa->GetBinError(b));
      }
      if (pv->GetBinEntries(b) > 1)
      {
        const int ip = gpv->GetN();
        gpv->SetPoint(ip, pv->GetBinCenter(b) + 0.12, pv->GetBinContent(b));
        gpv->SetPointError(ip, 0, 0, pv->GetBinError(b), pv->GetBinError(b));
        mx = std::max(mx, pv->GetBinContent(b) + pv->GetBinError(b));
      }
    }
    Frame(c, -0.5, 25.5, 0., 1.55 * mx, "true N", "#LT n #GT");
    auto *d = new TLine(0, 0, 25.5, 25.5);
    d->SetLineStyle(2);
    d->SetLineColor(kGray + 2);
    d->Draw();
    Style(gpa, cBlue, 20);
    Style(gpv, cRed, 24);
    gpa->Draw("P same");
    gpv->Draw("P same");
    auto *lg = Leg(0.18, 0.66, 0.55, 0.82, 0.032);
    lg->AddEntry(gpa, "all collisions", "pe");
    lg->AddEntry(gpv, "vertexed collisions", "pe");
    lg->AddEntry(d, "n = N", "l");
    lg->Draw();
    Label("Mean reconstructed n vs true N");
  }
  {
    TCanvas *c = NewPage(bResp, "n_minus_N");
    c->SetLogy();
    auto *a = (TH1D *) Get<TH1D>("hnMinusN_all")->Clone("nmN_all");
    auto *v = (TH1D *) Get<TH1D>("hnMinusN_vtx")->Clone("nmN_vtx");
    Style(a, cBlue, 20, 0.9);
    Style(v, cRed, 24, 0.9);
    for (TH1 *h : {(TH1 *) a, (TH1 *) v}) h->GetXaxis()->SetRangeUser(-15.5, 40.5);
    Headroom({a, v}, true, 1.6, 0.40);
    a->GetYaxis()->SetTitle("collisions");
    a->Draw("E1");
    v->Draw("E1 same");
    TLine l0(0., 0.5, 0., a->GetMaximum() * 0.02);
    l0.SetLineStyle(2);
    l0.SetLineColor(kGray + 2);
    l0.DrawClone();
    const double s = stat("events");
    auto *lg = Leg(0.55, 0.715, 0.93, 0.82, 0.031);
    lg->AddEntry(a, Form("all collisions (mean %.2f)", a->GetMean()), "pe");
    lg->AddEntry(v, Form("vertexed (mean %.2f)", v->GetMean()), "pe");
    lg->Draw();
    Info({Form("n > N: %.1f%%", 100. * stat("n>N") / s), Form("n < N: %.1f%%", 100. * stat("n<N") / s),
          Form("N = 0, n > 0: %.1f%%", 100. * stat("N=0 & n>0") / s)},
         0.18, 0.80, 0.029);
    Label("Reconstructed minus true multiplicity, n - N");
  }
  {
    TCanvas *c = NewPage(bResp, "n_over_N");
    auto *h = (TH1D *) Get<TH1D>("hnOverN_all")->Clone("noN");
    h->SetDirectory(nullptr);
    if (h->GetEntries() < 5000) h->Rebin(2);
    h->GetXaxis()->SetRangeUser(0., 8.);
    h->SetLineColor(cBlue);
    h->SetFillColor(cLightBlue);
    Headroom({h}, false, 1.6);
    h->GetYaxis()->SetTitle("collisions with N > 0");
    h->Draw("hist");
    h->Draw("E1 same");
    auto *l1 = new TLine(1., 0., 1., h->GetMaximum() / 1.6 * 1.05);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->SetLineColor(cRed);
    l1->Draw();
    auto *lg = Leg(0.58, 0.745, 0.93, 0.82, 0.031);
    lg->AddEntry(l1, "n = N", "l");
    lg->Draw();
    Info({Form("mean n / N = %.2f", h->GetMean())}, 0.60, 0.72, 0.031);
    Label("Ratio n / N (collisions with N > 0)");
  }
  {
    TCanvas *c = NewPage(bResp, "species");
    c->SetLogy();
    const char *nm[6] = {"hNpip", "hNpim", "hNKp", "hNKm", "hNp", "hNpbar"};
    const char *lb[6] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};
    const int col[6] = {cBlue, cBlue, cRed, cRed, cGreen, cGreen};
    const int mk[6] = {20, 24, 21, 25, 22, 26};
    std::vector<TH1 *> hs;
    for (int i = 0; i < 6; ++i)
    {
      auto *h = (TH1D *) Get<TH1D>(nm[i])->Clone(Form("%s_p", nm[i]));
      Style(h, col[i], mk[i], 0.9);
      h->GetXaxis()->SetRangeUser(-0.5, 15.5);
      hs.push_back(h);
    }
    Headroom(hs, true, 1.6, 0.38);
    auto *lg = Leg(0.42, 0.67, 0.93, 0.82, 0.030);
    lg->SetNColumns(2);
    for (int i = 0; i < 6; ++i)
    {
      hs[i]->GetXaxis()->SetTitle("primaries per collision (p_{T} > 0.1 GeV/c, |#eta| < 1.1)");
      hs[i]->GetYaxis()->SetTitle("collisions");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
      lg->AddEntry(hs[i], Form("%s (mean %.2f)", lb[i], hs[i]->GetMean()), "pe");
    }
    lg->Draw();
    Label("Identified primary multiplicities (generator level)");
  }

  // ===========================================================================
  // 3. TrackQA.pdf
  // ===========================================================================
  Book bTrk{"TrackQA.pdf", {}};
  // Four-curve shape comparison: all/good tracks x vertexed/vertexless events.
  auto overlay = [&](const char *cname, const char *base, const char *title, bool logy, double cutLine,
                     const char *cutLabel, double legX1, double legX2, double xlo, double xhi,
                     double freeFrac = 0.42, double logMin = 1e-5) {
    TCanvas *c = NewPage(bTrk, cname);
    if (logy) c->SetLogy();
    const char *var[4] = {"all_vtxEvt", "all_novtxEvt", "good_vtxEvt", "good_novtxEvt"};
    const char *lab[4] = {"all tracks, vertexed coll.", "all tracks, vertexless coll.", "good tracks, vertexed coll.",
                          "good tracks, vertexless coll."};
    const int col[4] = {cBlue, cRed, cBlue, cRed};
    const int mk[4] = {24, 25, 20, 21};
    std::vector<TH1 *> hs;
    for (int i = 0; i < 4; ++i)
    {
      auto *h = UnitArea(Get<TH1D>(Form("%s_%s", base, var[i])), Form("%s_%s_p", base, var[i]));
      Style(h, col[i], mk[i], 0.8);
      if (xhi > xlo) h->GetXaxis()->SetRangeUser(xlo, xhi);
      hs.push_back(h);
    }
    Headroom(hs, logy, 1.75, freeFrac, logMin);
    for (int i = 0; i < 4; ++i)
    {
      hs[i]->GetYaxis()->SetTitle("fraction of tracks");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
    }
    auto *lg = Leg(legX1, 0.635, legX2, 0.82, 0.029);
    for (int i = 0; i < 4; ++i) lg->AddEntry(hs[i], lab[i], "pe");
    TLine *cl = nullptr;
    if (cutLine > -1e8)
    {
      const double ymin = hs[0]->GetMinimum();
      const double ytop = logy ? ymin * std::pow(hs[0]->GetMaximum() / ymin, 0.56) : hs[0]->GetMaximum() * 0.56;
      cl = new TLine(cutLine, ymin, cutLine, ytop);
      cl->SetLineStyle(2);
      cl->SetLineWidth(2);
      cl->SetLineColor(cGreen);
      cl->Draw();
      lg->AddEntry(cl, cutLabel, "l");
      lg->SetY1NDC(0.60);
    }
    lg->Draw();
    Label(title);
  };
  // log pages: the data sit at the left, so the legend goes to the right with little extra head-room
  overlay("pt", "hPt", "Track p_{T} (unit area)", true, 0.1, "cut: p_{T} > 0.1 GeV/c", 0.50, 0.93, 0., 4., 0.12, 5e-5);
  overlay("eta", "hEta", "Track #eta (unit area)", false, -1e9, "", 0.52, 0.93, -1.5, 1.5);
  overlay("ntpc", "hNtpc", "TPC clusters per track (unit area)", false, 17.5, "cut: N_{TPC clus} #geq 18", 0.18, 0.60,
          0., 60.);
  overlay("nmvtx", "hNmvtx", "MVTX clusters per track (unit area)", true, -1e9, "", 0.52, 0.93, 0., 0., 0.10, 1e-4);
  overlay("nintt", "hNintt", "INTT clusters per track (unit area)", true, -1e9, "", 0.52, 0.93, 0., 0., 0.10, 1e-4);
  // grouped bars: silicon content of good tracks
  {
    TCanvas *c = NewPage(bTrk, "silicon_content");
    auto *v = UnitArea(Get<TH1D>("hSilicon_good_vtxEvt"), "sil_v");
    auto *nv = UnitArea(Get<TH1D>("hSilicon_good_novtxEvt"), "sil_nv");
    NoErrors(v);
    NoErrors(nv);
    v->SetFillColor(cBlue);
    v->SetLineColor(cBlue);
    nv->SetFillColor(cRed);
    nv->SetLineColor(cRed);
    v->SetBarWidth(0.36);
    v->SetBarOffset(0.13);
    nv->SetBarWidth(0.36);
    nv->SetBarOffset(0.51);
    v->SetMinimum(0.);
    v->SetMaximum(1.4);
    v->GetYaxis()->SetTitle("fraction of good tracks");
    v->GetXaxis()->SetTitle("");
    v->GetXaxis()->SetLabelSize(0.050);
    v->Draw("bar");
    nv->Draw("bar same");
    TLine one(0.5, 1., 3.5, 1.);
    one.SetLineStyle(3);
    one.SetLineColor(kGray + 2);
    one.DrawClone();
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.028);
    t.SetTextAlign(21);
    for (int b = 1; b <= 3; ++b)
    {
      t.SetTextColor(cBlue);
      t.DrawLatex(b - 0.5 + 0.13 + 0.18, v->GetBinContent(b) + 0.02, Form("%.2f", v->GetBinContent(b)));
      t.SetTextColor(cRed);
      t.DrawLatex(b - 0.5 + 0.51 + 0.18, nv->GetBinContent(b) + 0.02, Form("%.2f", nv->GetBinContent(b)));
    }
    auto *lg = Leg(0.56, 0.715, 0.93, 0.82, 0.032);
    lg->AddEntry(v, "vertexed collisions", "f");
    lg->AddEntry(nv, "vertexless collisions", "f");
    lg->Draw();
    Label("Silicon content of good tracks");
  }
  // stacked silicon fractions vs n
  {
    TCanvas *c = NewPage(bTrk, "silicon_vs_n");
    auto fr = CategoryFractions(Get<TH2D>("hSilicon_vs_n"), "silfrac", 31);
    auto *st = new THStack("st_sil", Form(";%s;fraction of good tracks", xn));
    const int col[3] = {cBlue, cOrange, cGray};
    const char *lab[3] = {"with MVTX clusters", "INTT, no MVTX", "TPC only (no silicon)"};
    for (int i = 0; i < 3; ++i)
    {
      fr[i]->SetFillColor(col[i]);
      fr[i]->SetLineColor(col[i]);
      st->Add(fr[i], "hist");
    }
    st->SetMaximum(1.4);
    st->SetMinimum(0.);
    st->Draw();
    st->GetXaxis()->SetRangeUser(0.5, 30.5);
    st->GetXaxis()->SetTitleSize(0.043);
    st->GetYaxis()->SetTitleSize(0.043);
    st->GetXaxis()->SetLabelSize(0.038);
    st->GetYaxis()->SetLabelSize(0.038);
    st->GetYaxis()->SetTitleOffset(1.40);
    st->GetXaxis()->SetTitleOffset(1.15);
    auto *lg = Leg(0.18, 0.73, 0.93, 0.82, 0.030);
    lg->SetNColumns(3);
    for (int i = 0; i < 3; ++i) lg->AddEntry(fr[i], lab[i], "f");
    lg->Draw();
    TLine one(0.5, 1., 30.5, 1.);
    one.SetLineStyle(3);
    one.SetLineColor(kGray + 2);
    one.DrawClone();
    c->Modified();
    Label("Silicon content of good tracks vs n");
  }
  auto perEvent = [&](const char *cname, std::vector<TH1 *> hs0, std::vector<std::string> labs, const char *xt,
                      const char *title, int xmaxBin) {
    TCanvas *c = NewPage(bTrk, cname);
    c->SetLogy();
    const int col[3] = {cBlack, cBlue, cRed};
    const int mk[3] = {20, 21, 24};
    std::vector<TH1 *> hs;
    for (size_t i = 0; i < hs0.size(); ++i)
    {
      auto *h = (TH1D *) hs0[i]->Clone(Form("%s_%zu", cname, i));
      Style(h, col[i], mk[i], 0.85);
      h->GetXaxis()->SetRange(1, xmaxBin);
      hs.push_back(h);
    }
    Headroom(hs, true, 1.6, 0.36);
    auto *lg = Leg(0.46, 0.67, 0.93, 0.82, 0.030);
    for (size_t i = 0; i < hs.size(); ++i)
    {
      hs[i]->GetXaxis()->SetTitle(xt);
      hs[i]->GetYaxis()->SetTitle("collisions");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
      lg->AddEntry(hs[i], Form("%s (mean %.2f)", labs[i].c_str(), hs[i]->GetMean()), "pe");
    }
    lg->Draw();
    Label(title);
  };
  perEvent("trk_per_evt", {Get<TH1D>("hnAll_all"), Get<TH1D>("hnAll_vtx"), Get<TH1D>("hnAll_novtx")},
           {"all collisions", "vertexed", "vertexless"}, "reconstructed tracks per collision (no cuts)",
           "All reconstructed tracks per collision", 61);
  perEvent("good_per_evt", {Get<TH1D>("hn_all"), Get<TH1D>("hnGoodMVTX_all"), Get<TH1D>("hnGoodNoMVTX_all")},
           {"all good tracks (n)", "good, with MVTX", "good, without MVTX"}, "good tracks per collision",
           "Good tracks per collision", 51);

  // ===========================================================================
  // 4. TruthMatchingQA.pdf
  // ===========================================================================
  Book bTru{"TruthMatchingQA.pdf", {}};
  {
    TCanvas *c = NewPage(bTru, "composition");
    auto *v = UnitArea(Get<TH1D>("hComp_good_vtxEvt"), "comp_v");
    auto *nv = UnitArea(Get<TH1D>("hComp_good_novtxEvt"), "comp_nv");
    auto *a0 = (TH1D *) Get<TH1D>("hComp_good_vtxEvt")->Clone("comp_a0");
    a0->Add(Get<TH1D>("hComp_good_novtxEvt"));
    auto *a = UnitArea(a0, "comp_a");
    TH1D *hs[3] = {a, v, nv};
    for (TH1D *h : hs) NoErrors(h);
    const int col[3] = {cGray, cBlue, cRed};
    const double off[3] = {0.08, 0.36, 0.64};
    double mx = 0;
    for (TH1D *h : hs) mx = std::max(mx, h->GetMaximum());
    for (int i = 0; i < 3; ++i)
    {
      hs[i]->SetFillColor(col[i]);
      hs[i]->SetLineColor(col[i]);
      hs[i]->SetBarWidth(0.28);
      hs[i]->SetBarOffset(off[i]);
      hs[i]->SetMinimum(0.);
      hs[i]->SetMaximum(1.8 * mx);
    }
    a->GetYaxis()->SetTitle("fraction of good tracks");
    a->GetXaxis()->SetTitle("");
    a->GetXaxis()->SetLabelSize(0.050);
    a->Draw("bar");
    v->Draw("bar same");
    nv->Draw("bar same");
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.024);
    t.SetTextAlign(21);
    for (int i = 0; i < 3; ++i)
      for (int b = 1; b <= 4; ++b)
      {
        t.SetTextColor(col[i] == cGray ? kGray + 2 : col[i]);
        t.DrawLatex(b - 0.5 + off[i] + 0.14, hs[i]->GetBinContent(b) + 0.01 * mx, Form("%.2f", hs[i]->GetBinContent(b)));
      }
    auto *lg = Leg(0.58, 0.69, 0.93, 0.82, 0.030);
    lg->AddEntry(a, "all collisions", "f");
    lg->AddEntry(v, "vertexed collisions", "f");
    lg->AddEntry(nv, "vertexless collisions", "f");
    lg->Draw();
    Info({"primary / secondary: first good", "  track of a G4 particle", "duplicate: further good tracks",
          "  of the same G4 particle"},
         0.18, 0.80, 0.025);
    Label("Truth origin of good reconstructed tracks");
  }
  {
    TCanvas *c = NewPage(bTru, "composition_vs_n");
    auto fr = CategoryFractions(Get<TH2D>("hComp_vs_n"), "compfrac", 31);
    auto *st = new THStack("st_comp", Form(";%s;fraction of good tracks", xn));
    const int col[4] = {cBlue, cOrange, cRed, cGray};
    const char *lab[4] = {"G4 primary", "G4 secondary", "duplicate", "unmatched"};
    for (int i = 0; i < 4; ++i)
    {
      fr[i]->SetFillColor(col[i]);
      fr[i]->SetLineColor(col[i]);
      st->Add(fr[i], "hist");
    }
    st->SetMaximum(1.4);
    st->SetMinimum(0.);
    st->Draw();
    st->GetXaxis()->SetRangeUser(0.5, 30.5);
    st->GetXaxis()->SetTitleSize(0.043);
    st->GetYaxis()->SetTitleSize(0.043);
    st->GetXaxis()->SetLabelSize(0.038);
    st->GetYaxis()->SetLabelSize(0.038);
    st->GetYaxis()->SetTitleOffset(1.40);
    st->GetXaxis()->SetTitleOffset(1.15);
    auto *lg = Leg(0.18, 0.73, 0.93, 0.82, 0.031);
    lg->SetNColumns(4);
    for (int i = 0; i < 4; ++i) lg->AddEntry(fr[i], lab[i], "f");
    lg->Draw();
    TLine one(0.5, 1., 30.5, 1.);
    one.SetLineStyle(3);
    one.SetLineColor(kGray + 2);
    one.DrawClone();
    c->Modified();
    Label("Truth composition of good tracks vs n");
  }
  {
    TCanvas *c = NewPage(bTru, "truth_mult");
    c->SetLogy();
    auto *p = (TH1D *) Get<TH1D>("hnGoodPrim_all")->Clone("gp");
    auto *s = (TH1D *) Get<TH1D>("hnGoodSec_all")->Clone("gs");
    auto *N = (TH1D *) Get<TH1D>("hN_all")->Clone("gN");
    Style(p, cBlue, 20, 0.9);
    Style(s, cOrange, 21, 0.85);
    Style(N, cRed, 24, 0.9);
    for (TH1 *h : {(TH1 *) p, (TH1 *) s, (TH1 *) N}) h->GetXaxis()->SetRange(1, 36);
    Headroom({p, s, N}, true, 1.6, 0.38);
    p->GetXaxis()->SetTitle("multiplicity per collision");
    p->GetYaxis()->SetTitle("collisions");
    p->Draw("E1");
    s->Draw("E1 same");
    N->Draw("E1 same");
    auto *lg = Leg(0.42, 0.67, 0.93, 0.82, 0.030);
    lg->AddEntry(p, Form("good, matched to G4 primary (mean %.2f)", p->GetMean()), "pe");
    lg->AddEntry(s, Form("good, matched to G4 secondary (mean %.2f)", s->GetMean()), "pe");
    lg->AddEntry(N, Form("true N (mean %.2f)", N->GetMean()), "pe");
    lg->Draw();
    Label("Truth-matched good-track multiplicities");
  }
  {
    TCanvas *c = NewPage(bTru, "prim_vs_n", true);
    c->SetLogz();
    auto *h = (TH2D *) Get<TH2D>("hGoodPrim_vs_n")->Clone("gpvn");
    h->GetXaxis()->SetRangeUser(-0.5, 40.5);
    h->GetYaxis()->SetRangeUser(-0.5, 40.5);
    h->GetXaxis()->SetTitle(xn);
    h->GetYaxis()->SetTitle("good tracks matched to G4 primaries");
    h->GetZaxis()->SetTitle("collisions");
    h->GetZaxis()->SetMoreLogLabels();
    h->SetMinimum(0.99);
    h->Draw("colz");
    auto *d = new TLine(-0.5, -0.5, 40.5, 40.5);
    d->SetLineStyle(2);
    d->SetLineColor(cRed);
    d->SetLineWidth(2);
    d->Draw();
    // upper-left triangle (more primaries than n) is kinematically empty
    auto *lg = Leg(0.18, 0.76, 0.52, 0.82, 0.030);
    lg->SetFillStyle(1001);
    lg->SetFillColor(kWhite);
    lg->AddEntry(d, "all good tracks primary", "l");
    lg->Draw();
    Label("Primary-matched good tracks vs n");
  }

  // ===========================================================================
  // Additional diagnostic pages (appended, so earlier page numbers are stable)
  // ===========================================================================
  // --- VertexEfficiency.pdf: residuals vs true z, transverse positions
  {
    TCanvas *c = NewPage(bEff, "res_vs_ztrue");
    Frame(c, -20., 20., -0.13, 0.20, "true collision z (cm)", "mean residual, reco - true (cm)");
    TLine z0(-20., 0., 20., 0.);
    z0.SetLineStyle(3);
    z0.SetLineColor(kGray + 2);
    z0.DrawClone();
    const char *nm[3] = {"hDz_vs_z", "hDx_vs_z", "hDy_vs_z"};
    const char *lb[3] = {"#Deltaz = z_{reco} - z_{true}", "#Deltax = x_{reco} - x_{true}", "#Deltay = y_{reco} - y_{true}"};
    const int col[3] = {cBlack, cBlue, cRed};
    const int mk[3] = {20, 21, 24};
    const double sh[3] = {-0.3, 0., 0.3};
    auto *lg = Leg(0.18, 0.66, 0.60, 0.82, 0.031);
    for (int i = 0; i < 3; ++i)
    {
      auto *pr = Get<TH2D>(nm[i])->ProfileX(Form("%s_pf", nm[i]));
      auto *g = new TGraphAsymmErrors();
      for (int b = 1; b <= pr->GetNbinsX(); ++b)
      {
        if (pr->GetBinEntries(b) < 3) continue;
        const int ip = g->GetN();
        g->SetPoint(ip, pr->GetBinCenter(b) + sh[i], pr->GetBinContent(b));
        g->SetPointError(ip, 0, 0, pr->GetBinError(b), pr->GetBinError(b));
      }
      Style(g, col[i], mk[i]);
      g->Draw("P same");
      lg->AddEntry(g, lb[i], "pe");
    }
    lg->Draw();
    Info({"mean per 2 cm bin of true z", "(bins with #geq 3 vertices)", "bars: error of the mean"}, 0.64, 0.80, 0.026);
    Label("Vertex residuals vs true collision z");
  }
  for (int xy = 0; xy < 2; ++xy)
  {
    TCanvas *c = NewPage(bEff, xy == 0 ? "beamspot_x" : "beamspot_y");
    auto *t0 = (TH1D *) Get<TH1D>(xy == 0 ? "hXtrue_vtx" : "hYtrue_vtx")->Clone(Form("bs_t_%d", xy));
    auto *r0 = (TH1D *) Get<TH1D>(xy == 0 ? "hXreco" : "hYreco")->Clone(Form("bs_r_%d", xy));
    t0->SetLineColor(cBlack);
    t0->SetLineWidth(2);
    r0->SetLineColor(cBlue);
    r0->SetFillColor(cLightBlue);
    Headroom({t0, r0}, false, 1.6);
    r0->GetXaxis()->SetTitle(xy == 0 ? "vertex x (cm)" : "vertex y (cm)");
    r0->GetYaxis()->SetTitle("vertexed collisions");
    r0->Draw("hist");
    t0->Draw("hist same");
    auto *lg = Leg(0.18, 0.70, 0.60, 0.82, 0.031);
    lg->AddEntry(t0, Form("true vertex (mean %.4f cm, RMS %.4f cm)", t0->GetMean(), t0->GetRMS()), "l");
    lg->AddEntry(r0, Form("reconstructed vertex (mean %.4f cm, RMS %.4f cm)", r0->GetMean(), r0->GetRMS()), "f");
    lg->Draw();
    Info({"simulated beam spot (G4_Input.C, pp_COLLISION): mean (x, y) = (-0.070, 0.141) cm, #sigma_{x,y} = 0.012 cm"},
         0.18, 0.645, 0.024);
    Label(xy == 0 ? "Transverse vertex position: x" : "Transverse vertex position: y");
  }

  // --- MultiplicityResponse.pdf: diagnostic response with unique primary matches
  colzPage("resp_uprim", Get<TH2D>("hResp_uprim"), "P(n_{uniq. prim.} | N)",
           "Diagnostic response: unique primary-matched good tracks", true,
           "unique primary-matched good tracks");
  {
    TCanvas *c = NewPage(bResp, "mean_resp_compare");
    const char *nm[3] = {"hResp_all", "hResp_uprim", "hResp_uprimAcc"};
    const char *lb[3] = {"baseline n (all good tracks)", "unique primary-matched good tracks",
                         "unique primary-matched, primary in N acceptance"};
    const int col[3] = {cBlue, cOrange, cGreen};
    const int mk[3] = {20, 21, 22};
    const double sh[3] = {-0.15, 0., 0.15};
    std::vector<TGraphAsymmErrors *> gs;
    double mx = 0;
    for (int i = 0; i < 3; ++i)
    {
      auto *pr = Get<TH2D>(nm[i])->ProfileX(Form("%s_cmp", nm[i]));
      auto *g = new TGraphAsymmErrors();
      for (int b = 1; b <= pr->FindBin(25.); ++b)
      {
        if (pr->GetBinEntries(b) < 2) continue;
        const int ip = g->GetN();
        g->SetPoint(ip, pr->GetBinCenter(b) + sh[i], pr->GetBinContent(b));
        g->SetPointError(ip, 0, 0, pr->GetBinError(b), pr->GetBinError(b));
        mx = std::max(mx, pr->GetBinContent(b) + pr->GetBinError(b));
      }
      Style(g, col[i], mk[i]);
      gs.push_back(g);
    }
    Frame(c, -0.5, 25.5, 0., 1.6 * mx, "true N", "mean reconstructed multiplicity");
    auto *d = new TLine(0, 0, 25.5, 25.5);
    d->SetLineStyle(2);
    d->SetLineColor(kGray + 2);
    d->Draw();
    auto *lg = Leg(0.18, 0.62, 0.70, 0.82, 0.029);
    for (int i = 0; i < 3; ++i)
    {
      gs[i]->Draw("P same");
      lg->AddEntry(gs[i], lb[i], "pe");
    }
    lg->AddEntry(d, "equal to N", "l");
    lg->Draw();
    Label("Mean multiplicity vs true N: baseline vs diagnostic definitions");
  }

  // --- TruthMatchingQA.pdf: independent classifications and duplicate diagnostics
  {
    TCanvas *c = NewPage(bTru, "origin_x_uniqueness", true);
    auto *h = (TH2D *) Get<TH2D>("hOriginUnique")->Clone("ou_p");
    const double tot = h->Integral();
    if (tot > 0) h->Scale(1. / tot);
    h->SetMarkerSize(2.0);
    h->SetMarkerColor(kBlack);
    gStyle->SetPaintTextFormat(".3f");
    h->GetZaxis()->SetTitle("fraction of good tracks");
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitle("A: truth origin");
    h->GetYaxis()->SetTitle("");
    c->SetLeftMargin(0.24);
    h->Draw("colz text");
    Label(Form("Origin (A, columns) x uniqueness (B, rows), %.0f good tracks", tot));
  }
  {
    TCanvas *c = NewPage(bTru, "track_classes");
    auto *a = UnitArea(Get<TH1D>("hTrkClass_all"), "cls_all");
    auto *z = UnitArea(Get<TH1D>("hTrkClass_N0"), "cls_N0");
    NoErrors(a);
    NoErrors(z);
    const char *shortLab[5] = {"unmatched", "additional", "secondary", "prim. out acc.", "prim. in acc."};
    for (int b = 1; b <= 5; ++b) a->GetXaxis()->SetBinLabel(b, shortLab[b - 1]);
    a->GetXaxis()->LabelsOption("h");
    a->SetFillColor(cBlue);
    a->SetLineColor(cBlue);
    z->SetFillColor(cRed);
    z->SetLineColor(cRed);
    a->SetBarWidth(0.38);
    a->SetBarOffset(0.11);
    z->SetBarWidth(0.38);
    z->SetBarOffset(0.51);
    const double mx = std::max(a->GetMaximum(), z->GetMaximum());
    a->SetMinimum(0.);
    a->SetMaximum(1.7 * mx);
    a->GetYaxis()->SetTitle("fraction of good tracks");
    a->GetXaxis()->SetTitle("");
    a->GetXaxis()->SetLabelSize(0.040);
    a->Draw("bar");
    z->Draw("bar same");
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.025);
    t.SetTextAlign(21);
    for (int b = 1; b <= 5; ++b)
    {
      t.SetTextColor(cBlue);
      t.DrawLatex(b - 0.5 + 0.11 + 0.19, a->GetBinContent(b) + 0.01 * mx, Form("%.2f", a->GetBinContent(b)));
      t.SetTextColor(cRed);
      t.DrawLatex(b - 0.5 + 0.51 + 0.19, z->GetBinContent(b) + 0.01 * mx, Form("%.2f", z->GetBinContent(b)));
    }
    auto *lg = Leg(0.56, 0.71, 0.93, 0.82, 0.030);
    lg->AddEntry(a, Form("all collisions (%.0f tracks)", Get<TH1D>("hTrkClass_all")->Integral()), "f");
    lg->AddEntry(z, Form("collisions with N = 0 (%.0f tracks)", Get<TH1D>("hTrkClass_N0")->Integral()), "f");
    lg->Draw();
    Info({"exclusive classes,", "priority left to right", "(acc. = N acceptance)"}, 0.18, 0.80, 0.026);
    Label("Exclusive origin classes of good tracks");
  }
  auto dupPage = [&](const char *cname, const char *base, const char *xt, const char *title, bool logy, double xlo,
                     double xhi) {
    TCanvas *c = NewPage(bTru, cname);
    if (logy) c->SetLogy();
    const char *nm[3] = {"single", "multiFirst", "multiAdd"};
    const char *lb[3] = {"particle with exactly one good track", "multiply matched: first track",
                         "multiply matched: additional tracks"};
    const int col[3] = {cBlack, cBlue, cRed};
    const int mk[3] = {20, 21, 24};
    std::vector<TH1 *> hs;
    for (int i = 0; i < 3; ++i)
    {
      auto *h = UnitArea(Get<TH1D>(Form("%s_%s", base, nm[i])), Form("%s_%s_p", base, nm[i]));
      Style(h, col[i], mk[i], 0.85);
      if (xhi > xlo) h->GetXaxis()->SetRangeUser(xlo, xhi);
      hs.push_back(h);
    }
    Headroom(hs, logy, 1.7, 0.40, 1e-4);
    auto *lg = Leg(0.18, 0.66, 0.70, 0.82, 0.029);
    for (int i = 0; i < 3; ++i)
    {
      hs[i]->GetXaxis()->SetTitle(xt);
      hs[i]->GetYaxis()->SetTitle("fraction of tracks");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
      lg->AddEntry(hs[i], Form("%s (mean %.2f, %.0f tracks)", lb[i], hs[i]->GetMean(),
                               Get<TH1D>(Form("%s_%s", base, nm[i]))->GetEntries()),
                   "pe");
    }
    lg->Draw();
    Label(title);
  };
  dupPage("dup_purity", "hPurity", "truth purity = nfromtruth / nhits", "Truth purity of good tracks by match multiplicity",
          false, 0., 1.0001);
  dupPage("dup_ntpc", "hNtpcDup", "TPC clusters per track", "TPC clusters by match multiplicity", false, 15., 55.);
  dupPage("dup_quality", "hQualDup", "track quality (#chi^{2}/ndf)", "Track-fit quality by match multiplicity", true, 0., 20.);
  dupPage("dup_mvtx", "hMvtxDup", "MVTX clusters per track", "MVTX clusters by match multiplicity", true, 0., 0.);
  {
    TCanvas *c = NewPage(bTru, "matches_per_particle");
    c->SetLogy();
    auto *h = (TH1D *) Get<TH1D>("hMatchesPerParticle")->Clone("mpp");
    h->SetLineColor(cBlue);
    h->SetFillColor(cLightBlue);
    Headroom({h}, true, 1.6, 0.36);
    h->GetXaxis()->SetRangeUser(0.5, 6.5);
    h->GetYaxis()->SetTitle("matched G4 particles");
    h->Draw("hist");
    double tot = h->Integral(), multi = h->Integral(2, h->GetNbinsX());
    Info({Form("particles with > 1 good track: %.1f%%", tot > 0 ? 100. * multi / tot : 0.),
          Form("mean good tracks per matched particle: %.3f", h->GetMean())},
         0.45, 0.79, 0.030);
    Label("Good reconstructed tracks per matched G4 particle");
  }
  {
    TCanvas *c = NewPage(bTru, "dup_nfromtruth_diff");
    auto *h = (TH1D *) Get<TH1D>("hAddMinusFirst_nfromtruth")->Clone("nftd");
    h->GetXaxis()->SetRangeUser(-50.5, 50.5);
    h->SetLineColor(cRed);
    h->SetFillColor(TColor::GetColor("#f0b7b0"));
    Headroom({h}, false, 1.6);
    h->GetYaxis()->SetTitle("additional tracks");
    h->Draw("hist");
    const double tot = h->Integral(0, h->GetNbinsX() + 1);
    const double neg = h->Integral(0, h->FindBin(-0.6));
    const double zer = h->GetBinContent(h->FindBin(0.));
    Info({Form("additional < first: %.1f%%", tot > 0 ? 100. * neg / tot : 0.),
          Form("additional = first: %.1f%%", tot > 0 ? 100. * zer / tot : 0.),
          Form("additional > first: %.1f%%", tot > 0 ? 100. * (tot - neg - zer) / tot : 0.)},
         0.18, 0.79, 0.030);
    Label("Truth clusters: additional track minus first track of the same particle");
  }
  {
    TCanvas *c = NewPage(bTru, "category_mult");
    c->SetLogy();
    const char *nm[5] = {"hn_all", "hnUniquePrim_all", "hnUniqueSec_all", "hnAdd_all", "hnUnm_all"};
    const char *lb[5] = {"all good tracks (n)", "unique primary matches", "unique secondary matches",
                         "additional reco matches", "unmatched"};
    const int col[5] = {cBlack, cBlue, cOrange, cRed, cGray};
    const int mk[5] = {20, 21, 22, 24, 25};
    std::vector<TH1 *> hs;
    for (int i = 0; i < 5; ++i)
    {
      auto *h = (TH1D *) Get<TH1D>(nm[i])->Clone(Form("%s_cat", nm[i]));
      Style(h, col[i], mk[i], 0.85);
      h->GetXaxis()->SetRange(1, 41);
      hs.push_back(h);
    }
    Headroom(hs, true, 1.6, 0.42);
    auto *lg = Leg(0.42, 0.60, 0.93, 0.82, 0.029);
    for (int i = 0; i < 5; ++i)
    {
      hs[i]->GetXaxis()->SetTitle("tracks per collision");
      hs[i]->GetYaxis()->SetTitle("collisions");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
      lg->AddEntry(hs[i], Form("%s (mean %.2f)", lb[i], hs[i]->GetMean()), "pe");
    }
    lg->Draw();
    Label("Per-collision multiplicity of good-track categories");
  }

  // --- SiliconVertexQA.pdf (new): vertex finding vs silicon content
  Book bSil{"SiliconVertexQA.pdf", {}};
  {
    TEfficiency *em1 = Eff(Get<TH1D>("hnM1_vtx"), Get<TH1D>("hnM1_all"), "eff_m1");
    TEfficiency *em2 = Eff(Get<TH1D>("hnM2_vtx"), Get<TH1D>("hnM2_all"), "eff_m2");
    TEfficiency *em3 = Eff(Get<TH1D>("hnM3_vtx"), Get<TH1D>("hnM3_all"), "eff_m3");
    TCanvas *c = NewPage(bSil, "eff_nMVTX");
    EffFrame(c, -0.5, 10.5, "number of good tracks with #geq k MVTX clusters", "#varepsilon_{vtx}");
    TGraphAsymmErrors *g1 = EffGraph(em1, 10.5, -0.15), *g2 = EffGraph(em2, 10.5, 0.), *g3 = EffGraph(em3, 10.5, 0.15);
    Style(g1, cBlue, 20);
    Style(g2, cRed, 21);
    Style(g3, cBlack, 24);
    g1->Draw("P same");
    g2->Draw("P same");
    g3->Draw("P same");
    auto *lg = Leg(0.56, 0.66, 0.93, 0.82, 0.031);
    lg->AddEntry(g1, "n_{MVTX1} (k = 1)", "pe");
    lg->AddEntry(g2, "n_{MVTX2} (k = 2)", "pe");
    lg->AddEntry(g3, "n_{MVTX3} (k = 3)", "pe");
    lg->Draw();
    Info({"diagnostic: good tracks (baseline cuts)", "with #geq k MVTX clusters; all true z"}, 0.18, 0.80, 0.027);
    Label("Vertex-finding efficiency vs MVTX-associated multiplicity");
  }
  {
    TEfficiency *e0 = Eff(Get<TH1D>("hn_vtx_m0"), Get<TH1D>("hn_all_m0"), "eff_n_m0");
    TEfficiency *e1 = Eff(Get<TH1D>("hn_vtx_m1"), Get<TH1D>("hn_all_m1"), "eff_n_m1");
    TEfficiency *e2 = Eff(Get<TH1D>("hn_vtx_m2"), Get<TH1D>("hn_all_m2"), "eff_n_m2");
    TCanvas *c = NewPage(bSil, "eff_n_by_mvtx_class");
    EffFrame(c, -0.5, 30.5, xn, "#varepsilon_{vtx}(n)");
    TGraphAsymmErrors *g0 = EffGraph(e0, 30.5, -0.2), *g1 = EffGraph(e1, 30.5, 0.), *g2 = EffGraph(e2, 30.5, 0.2);
    Style(g0, cGray, 25);
    Style(g1, cBlue, 20);
    Style(g2, cRed, 21);
    g0->Draw("P same");
    g1->Draw("P same");
    g2->Draw("P same");
    auto *lg = Leg(0.52, 0.66, 0.93, 0.82, 0.031);
    lg->AddEntry(g0, Form("n_{MVTX1} = 0 (%.0f coll.)", stat("n_MVTX1=0")), "pe");
    lg->AddEntry(g1, Form("n_{MVTX1} #geq 1 (%.0f coll.)", stat("n_MVTX1>=1")), "pe");
    lg->AddEntry(g2, Form("n_{MVTX1} #geq 2 (%.0f coll.)", stat("n_MVTX1>=2")), "pe");
    lg->Draw();
    Info({"n: baseline good-track multiplicity", "classes by number of good tracks", "with #geq 1 MVTX cluster"}, 0.18,
         0.80, 0.027);
    Label("#varepsilon_{vtx}(n) by MVTX content of the event");
  }
  {
    TCanvas *c = NewPage(bSil, "eff2D_n_nMVTX2", true);
    auto *a = Get<TH2D>("hEff2D_all");
    auto *v = Get<TH2D>("hEff2D_vtx");
    auto *e = (TH2D *) v->Clone("eff2D");
    e->Reset();
    for (int i = 1; i <= a->GetNbinsX(); ++i)
      for (int j = 1; j <= a->GetNbinsY(); ++j)
        if (a->GetBinContent(i, j) > 0) e->SetBinContent(i, j, v->GetBinContent(i, j) / a->GetBinContent(i, j));
        else e->SetBinContent(i, j, -1);  // empty cell (N_all = 0): below the minimum -> not drawn
    // a slightly negative minimum makes ROOT paint cells with efficiency exactly 0
    e->SetMinimum(-1e-6);
    e->SetMaximum(1.);
    e->GetXaxis()->SetRangeUser(-0.5, 25.5);
    e->GetYaxis()->SetRangeUser(-0.5, 6.5);
    e->GetXaxis()->SetTitle(xn);
    e->GetYaxis()->SetTitle("n_{MVTX2}: good tracks with #geq 2 MVTX clusters");
    e->GetZaxis()->SetTitle("#varepsilon_{vtx} = N_{vtx} / N_{all} per cell");
    e->Draw("colz");
    Info({"empty cells (N_{all} = 0) left white; numbers of", "collisions per cell in SiliconQA tables"}, 0.18, 0.80,
         0.026);
    Label("#varepsilon_{vtx} in bins of (n, n_{MVTX2})");
  }
  {
    TCanvas *c = NewPage(bSil, "mvtx_mult_dists");
    c->SetLogy();
    const char *nm[3] = {"hnM1_all", "hnM2_all", "hnM3_all"};
    const char *lb[3] = {"n_{MVTX1}", "n_{MVTX2}", "n_{MVTX3}"};
    const int col[3] = {cBlue, cRed, cBlack};
    const int mk[3] = {20, 21, 24};
    std::vector<TH1 *> hs;
    for (int i = 0; i < 3; ++i)
    {
      auto *h = (TH1D *) Get<TH1D>(nm[i])->Clone(Form("%s_d", nm[i]));
      Style(h, col[i], mk[i], 0.9);
      h->GetXaxis()->SetRange(1, 16);
      hs.push_back(h);
    }
    Headroom(hs, true, 1.6, 0.36);
    auto *lg = Leg(0.52, 0.67, 0.93, 0.82, 0.031);
    for (int i = 0; i < 3; ++i)
    {
      hs[i]->GetXaxis()->SetTitle("good tracks with #geq k MVTX clusters per collision");
      hs[i]->GetYaxis()->SetTitle("collisions");
      hs[i]->Draw(i == 0 ? "E1" : "E1 same");
      lg->AddEntry(hs[i], Form("%s (mean %.3f)", lb[i], hs[i]->GetMean()), "pe");
    }
    lg->Draw();
    Label("MVTX-associated good-track multiplicities");
  }

  // ---------------------------------------------------------------------------
  // Print the four PDFs and the combined PDF
  // ---------------------------------------------------------------------------
  std::vector<Book *> books = {&bEff, &bResp, &bTrk, &bTru, &bSil};
  gManifest.open(std::string(plotDir) + "/figure_manifest.tsv");
  gManifest << "pdf\tpage\tcanvas\ttitle\n";
  for (Book *b : books) PrintBook(*b, plotDir);
  Book all{"VertexEfficiency_CompleteQA.pdf", {}};
  for (Book *b : books)
    for (TCanvas *c : b->pages) all.pages.push_back(c);
  PrintBook(all, plotDir);
  gManifest.close();

  // ---------------------------------------------------------------------------
  // Tables
  // ---------------------------------------------------------------------------
  std::ofstream t(tableFile);
  const TNamed *cfg = dynamic_cast<TNamed *>(gIn->Get("config"));
  t << "Vertex-finding efficiency tables\n";
  t << "Input : " << inFile << "\n";
  t << "Sample: pileup-free PYTHIA8 MB (Detroit), polyseed tracking\n";
  if (cfg) t << "Config: " << cfg->GetTitle() << "\n";
  t << "Intervals: Clopper-Pearson, 68.3% CL. UL95 = one-sided 95% CL upper limit for bins with 0 vertexed.\n\n";
  t << "Event statistics (hStats):\n";
  for (int b = 1; b <= hStats->GetNbinsX(); ++b)
    t << Form("  %-20s %12.0f\n", hStats->GetXaxis()->GetBinLabel(b), hStats->GetBinContent(b));
  t << "\nepsilon_vtx(n), fine bins:\n";
  WriteEffTable(t, en_fine, "n", false);
  t << "\nepsilon_vtx(n), wide bins:\n";
  WriteEffTable(t, en_wide, "n", true);
  t << "\nepsilon_vtx(N), fine bins:\n";
  WriteEffTable(t, eN_fine, "N", false);
  t << "\nepsilon_vtx(N), wide bins:\n";
  WriteEffTable(t, eN_wide, "N", true);
  t << "\nepsilon_vtx(n), fine bins, " << zcText << ":\n";
  WriteEffTable(t, en_fine_zc, "n", false);
  t << "\nepsilon_vtx(n), wide bins, " << zcText << ":\n";
  WriteEffTable(t, en_wide_zc, "n", true);
  t << "\nepsilon_vtx(N), wide bins, " << zcText << ":\n";
  WriteEffTable(t, eN_wide_zc, "N", true);
  t << "\nepsilon_vtx versus true z (all multiplicities):\n";
  WriteEffTable(t, eZ, "z (cm)", false, true);
  // ---- silicon diagnostics
  t << "\n=== SILICON DIAGNOSTICS (n_MVTXk = good tracks with >= k MVTX clusters) ===\n";
  t << "\nepsilon_vtx(n_MVTX1):\n";
  WriteEffTable(t, Eff(Get<TH1D>("hnM1_vtx"), Get<TH1D>("hnM1_all"), "t_m1"), "n_MVTX1", false);
  t << "\nepsilon_vtx(n_MVTX2):\n";
  WriteEffTable(t, Eff(Get<TH1D>("hnM2_vtx"), Get<TH1D>("hnM2_all"), "t_m2"), "n_MVTX2", false);
  t << "\nepsilon_vtx(n_MVTX3):\n";
  WriteEffTable(t, Eff(Get<TH1D>("hnM3_vtx"), Get<TH1D>("hnM3_all"), "t_m3"), "n_MVTX3", false);
  t << "\nepsilon_vtx(n), collisions with n_MVTX1 = 0:\n";
  WriteEffTable(t, Eff(Get<TH1D>("hn_vtx_m0"), Get<TH1D>("hn_all_m0"), "t_n_m0"), "n", false);
  t << "\nepsilon_vtx(n), collisions with n_MVTX1 >= 1:\n";
  WriteEffTable(t, Eff(Get<TH1D>("hn_vtx_m1"), Get<TH1D>("hn_all_m1"), "t_n_m1"), "n", false);
  t << "\nepsilon_vtx(n), collisions with n_MVTX1 >= 2:\n";
  WriteEffTable(t, Eff(Get<TH1D>("hn_vtx_m2"), Get<TH1D>("hn_all_m2"), "t_n_m2"), "n", false);
  {
    t << "\n2D map N_vtx/N_all in cells (n, n_MVTX2), n = 0..20, n_MVTX2 = 0..4 (empty: -):\n";
    auto *a = Get<TH2D>("hEff2D_all");
    auto *v = Get<TH2D>("hEff2D_vtx");
    t << "  n\\nM2";
    for (int j = 1; j <= 5; ++j) t << Form("  %12d", j - 1);
    t << "\n";
    for (int i = 1; i <= 21; ++i)
    {
      t << Form("  %5d", i - 1);
      for (int j = 1; j <= 5; ++j)
      {
        const double na = a->GetBinContent(i, j), nv = v->GetBinContent(i, j);
        if (na > 0) t << Form("  %5.0f/%-6.0f", nv, na);
        else t << Form("  %12s", "-");
      }
      t << "\n";
    }
  }
  // ---- truth classifications
  t << "\n=== TRUTH CLASSIFICATIONS OF GOOD TRACKS (QA only) ===\n";
  {
    auto *ou = Get<TH2D>("hOriginUnique");
    const double tot = ou->Integral();
    t << "\nA (origin) x B (uniqueness), counts and fraction of all good tracks:\n";
    t << Form("  %-18s %22s %22s %22s %10s\n", "", "unique match", "additional match", "unmatched", "sum");
    for (int i = 1; i <= 3; ++i)
    {
      double row = 0;
      t << Form("  %-18s", ou->GetXaxis()->GetBinLabel(i));
      for (int j = 1; j <= 3; ++j)
      {
        row += ou->GetBinContent(i, j);
        t << Form(" %10.0f (%8.4f)  ", ou->GetBinContent(i, j), tot > 0 ? ou->GetBinContent(i, j) / tot : 0.);
      }
      t << Form(" %10.0f\n", row);
    }
    auto *ca = Get<TH1D>("hTrkClass_all");
    auto *c0 = Get<TH1D>("hTrkClass_N0");
    t << "\nExclusive classes (priority unmatched > additional > secondary > primary out of N acc. > primary in N acc.):\n";
    t << Form("  %-26s %12s %10s %14s %10s\n", "class", "all coll.", "fraction", "N=0 coll.", "fraction");
    for (int b = 1; b <= 5; ++b)
      t << Form("  %-26s %12.0f %10.4f %14.0f %10.4f\n", ca->GetXaxis()->GetBinLabel(b), ca->GetBinContent(b),
                ca->Integral() > 0 ? ca->GetBinContent(b) / ca->Integral() : 0., c0->GetBinContent(b),
                c0->Integral() > 0 ? c0->GetBinContent(b) / c0->Integral() : 0.);
    t << "\nMean per-collision multiplicities:\n";
    const char *nm[8] = {"hn_all", "hN_all", "hnGoodPrim_all", "hnGoodSec_all", "hnUniquePrim_all", "hnUniqueSec_all",
                         "hnAdd_all", "hnUnm_all"};
    const char *lb[8] = {"n (all good tracks)", "N (true)", "good, primary-matched (incl. add.)",
                         "good, secondary-matched (incl. add.)", "unique primary matches", "unique secondary matches",
                         "additional reco matches", "unmatched"};
    for (int i = 0; i < 8; ++i) t << Form("  %-40s %8.4f\n", lb[i], Get<TH1D>(nm[i])->GetMean());
    auto *mp = Get<TH1D>("hMatchesPerParticle");
    t << Form("\nMatched G4 particles: %.0f, with > 1 good track: %.0f (%.2f%%); mean good tracks per particle %.4f\n",
              mp->Integral(), mp->Integral(2, mp->GetNbinsX()),
              mp->Integral() > 0 ? 100. * mp->Integral(2, mp->GetNbinsX()) / mp->Integral() : 0., mp->GetMean());
    const char *dn[3] = {"single", "multiFirst", "multiAdd"};
    t << "\nDuplicate diagnostics (means):            single   multi-first   multi-additional\n";
    const char *qn[4] = {"hPurity", "hNtpcDup", "hQualDup", "hMvtxDup"};
    const char *ql[4] = {"truth purity nfromtruth/nhits", "TPC clusters", "quality chi2/ndf (range 0-20)", "MVTX clusters"};
    for (int q = 0; q < 4; ++q)
    {
      t << Form("  %-38s", ql[q]);
      for (int i = 0; i < 3; ++i) t << Form(" %10.3f", Get<TH1D>(Form("%s_%s", qn[q], dn[i]))->GetMean());
      t << "\n";
    }
    t << Form("  %-38s %10.0f %10.0f %10.0f\n", "tracks", Get<TH1D>("hPurity_single")->GetEntries(),
              Get<TH1D>("hPurity_multiFirst")->GetEntries(), Get<TH1D>("hPurity_multiAdd")->GetEntries());
    auto *df = Get<TH1D>("hAddMinusFirst_nfromtruth");
    const double tot2 = df->Integral(0, df->GetNbinsX() + 1);
    const double neg = df->Integral(0, df->FindBin(-0.6)), zer = df->GetBinContent(df->FindBin(0.));
    t << Form("  nfromtruth(additional) - nfromtruth(first): <0 %.2f%%, =0 %.2f%%, >0 %.2f%% (mean %.3f)\n",
              tot2 > 0 ? 100. * neg / tot2 : 0., tot2 > 0 ? 100. * zer / tot2 : 0.,
              tot2 > 0 ? 100. * (tot2 - neg - zer) / tot2 : 0., df->GetMean());
  }
  // ---- residuals
  t << "\n=== VERTEX RESIDUALS (vertexed collisions) ===\n";
  {
    const char *nm[3] = {"hDz", "hDx", "hDy"};
    for (int i = 0; i < 3; ++i)
    {
      auto *h = Get<TH1D>(nm[i]);
      t << Form("  %s: entries %.0f, mean %+.4f +- %.4f cm, RMS %.4f cm, underflow %.0f, overflow %.0f\n", nm[i],
                h->GetEntries(), h->GetMean(), h->GetMeanError(), h->GetRMS(), h->GetBinContent(0),
                h->GetBinContent(h->GetNbinsX() + 1));
    }
    t << Form("  true  x mean %.4f cm (RMS %.4f), y mean %.4f cm (RMS %.4f)\n", Get<TH1D>("hXtrue_vtx")->GetMean(),
              Get<TH1D>("hXtrue_vtx")->GetRMS(), Get<TH1D>("hYtrue_vtx")->GetMean(), Get<TH1D>("hYtrue_vtx")->GetRMS());
    t << Form("  reco  x mean %.4f cm (RMS %.4f), y mean %.4f cm (RMS %.4f)\n", Get<TH1D>("hXreco")->GetMean(),
              Get<TH1D>("hXreco")->GetRMS(), Get<TH1D>("hYreco")->GetMean(), Get<TH1D>("hYreco")->GetRMS());
    t << "\n  Mean residuals in 4 cm bins of true z:\n";
    t << Form("  %-12s %8s %16s %16s %16s\n", "z_true (cm)", "N", "dz (cm)", "dx (cm)", "dy (cm)");
    TProfile *pz = Get<TH2D>("hDz_vs_z")->ProfileX("tz"), *px = Get<TH2D>("hDx_vs_z")->ProfileX("tx"),
             *py = Get<TH2D>("hDy_vs_z")->ProfileX("ty");
    pz->Rebin(2);
    px->Rebin(2);
    py->Rebin(2);
    for (int b = 1; b <= pz->GetNbinsX(); ++b)
    {
      if (pz->GetBinEntries(b) <= 0) continue;
      t << Form("  [%+3.0f,%+3.0f)   %8.0f  %+.4f+-%.4f  %+.4f+-%.4f  %+.4f+-%.4f\n", pz->GetXaxis()->GetBinLowEdge(b),
                pz->GetXaxis()->GetBinUpEdge(b), pz->GetBinEntries(b), pz->GetBinContent(b), pz->GetBinError(b),
                px->GetBinContent(b), px->GetBinError(b), py->GetBinContent(b), py->GetBinError(b));
    }
  }

  t << "\nMean reco n per true N (all collisions):\n";
  {
    auto *pa = hR->ProfileX("prof_tab");
    t << Form("  %4s %10s %8s %8s\n", "N", "events", "<n>", "err");
    for (int b = 1; b <= pa->GetNbinsX(); ++b)
      if (pa->GetBinEntries(b) > 0)
        t << Form("  %4.0f %10.0f %8.3f %8.3f\n", pa->GetBinCenter(b), pa->GetBinEntries(b), pa->GetBinContent(b),
                  pa->GetBinError(b));
  }
  t.close();
  std::cout << "Plots written to " << plotDir << "/ (" << all.pages.size() << " pages), tables to " << tableFile
            << std::endl;
}
