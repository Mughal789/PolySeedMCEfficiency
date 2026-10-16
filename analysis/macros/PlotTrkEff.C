// =============================================================================
// PlotTrkEff.C - plots of the single-particle tracking efficiency
//
// Input : the output of TrkEffAnalysis (raw-count histograms and TEfficiency
//         objects, see analysis/TRKEFF.md). The file is opened READ-ONLY and
//         nothing is recomputed: every efficiency shown comes from the stored
//         TEfficiency objects (Clopper-Pearson, 68.3% CL).
// Output (in OUTDIR):
//   <pdfName>                     multipage PDF (summary, maps, comparisons, diagnostics)
//   figures/*.pdf, figures/*.png  every panel as a stand-alone figure (for slides)
//   TrackEfficiency_plots.root    every plotted object (graphs, maps) - used by
//                                 ValidateTrkEffPlots.C
//   TrackEfficiency_numbers.txt   integrated efficiencies and empty-bin counts
//
// Statistical treatment:
//   - 1D: TEfficiency::CreateGraph -> points only where the denominator > 0,
//     asymmetric Clopper-Pearson errors. Undefined bins are absent, never 0.
//   - 2D: colour = passed/total where total > 0; bins with total = 0 are
//     covered by a hatched grey box ("undefined"). Each defined bin shows its
//     counts passed/total, so low-statistics bins are visible as such.
//   - Integrated efficiency: summed counts over the in-range bins, with the
//     Clopper-Pearson 68.3% interval of the sums.
//   - No smoothing, fitting or interpolation.
//
// Usage:
//   root -l -b -q 'PlotTrkEff.C+("TrkEff_baseline.root","OUTDIR","K1, 1000 ev","TrackEfficiency_K1.pdf","skim_files.list")'
//   (last argument optional: skim list for the low-pT diagnostics page; "" = skip)
// =============================================================================
#include "TrkEffConfigUtil.h"

#include "LowPtDiagnostics.C"
#include "TrkEffPlotStyle.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace trkeffplot;

namespace plottrkeff
{
// Everything the plots need, read once from the input file and detached from it.
struct Inputs
{
  std::map<std::string, TEfficiency *> eff;  // e.g. "eff_pip_ypt", "effQ_km_pt"
  std::map<std::string, TH1 *> hist;         // e.g. "hDen_pip_pt"
  TH1 *stats = nullptr;
  std::string config, definitions;
};

Inputs ReadInputs(const char *path)
{
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error(std::string("cannot open ") + path);
  Inputs in;
  for (const auto &s : Species())
  {
    for (const char *kind : {"eff", "effQ"})
      for (const char *proj : {"ypt", "pt", "y"})
      {
        const std::string n = Form("%s_%s_%s", kind, s.tag, proj);
        auto *e = dynamic_cast<TEfficiency *>(f->Get(n.c_str()));
        if (!e) throw std::runtime_error("missing " + n + " in " + path);
        e = static_cast<TEfficiency *>(e->Clone());
        e->SetDirectory(nullptr);
        in.eff[n] = e;
      }
    for (const char *kind : {"hDen", "hNum", "hNumQ"})
      for (const char *proj : {"ypt", "pt", "y"})
      {
        const std::string n = Form("%s_%s_%s", kind, s.tag, proj);
        auto *h = dynamic_cast<TH1 *>(f->Get(n.c_str()));
        if (!h) throw std::runtime_error("missing " + n + " in " + path);
        h = static_cast<TH1 *>(h->Clone());
        h->SetDirectory(nullptr);
        in.hist[n] = h;
      }
  }
  in.stats = static_cast<TH1 *>(f->Get("hTrkEffStats")->Clone());
  in.stats->SetDirectory(nullptr);
  in.config = static_cast<TNamed *>(f->Get("config"))->GetTitle();
  in.definitions = static_cast<TNamed *>(f->Get("definitions"))->GetTitle();
  return in;  // f is closed here; the input file is never written
}

double Stat(const Inputs &in, const char *label) { return in.stats->GetBinContent(in.stats->GetXaxis()->FindBin(label)); }

// Output bookkeeping: the multipage PDF, single figures and the plotted objects.
struct Output
{
  std::string dir, pdf;
  TFile *objects = nullptr;
  int page = 0;
  bool pagesPass = true;  // pass 1: multipage PDF + plotted objects; pass 2: single figures only
  mutable std::set<std::string> written;  // each plotted object is stored once
  // Single figures are written only in pass 2: ROOT has one active PDF stream,
  // so they must not be written while the multipage PDF is open.
  void SaveFigure(TCanvas *c, const std::string &name) const
  {
    if (pagesPass) return;
    c->SaveAs((dir + "/figures/" + name + ".pdf").c_str());
    c->SaveAs((dir + "/figures/" + name + ".png").c_str());
  }
  void Keep(TObject *o) const
  {
    if (!pagesPass || !written.insert(o->GetName()).second) return;
    objects->cd();
    o->Write();
  }
};

// ---- panels: each draws into the current pad -------------------------------------
// 2D efficiency map (efficiency colour, masked undefined bins, counts per bin).
void DrawMapPanel(const TEfficiency *e, const SpeciesStyle &s, const char *what, const std::string &sample,
                  const Output &out)
{
  gPad->SetRightMargin(0.16);
  gPad->SetTopMargin(0.235);  // header and statistics above the frame, never on top of data
  gPad->SetLogy();
  TH2D *m = EfficiencyMap(e, Form("map_%s", e->GetName()));
  m->SetBit(kCanDelete);
  m->SetMinimum(-1e-6);  // a measured efficiency of 0 is painted (darkest colour)
  m->SetMaximum(1.0);
  m->GetXaxis()->SetTitle("y_{true}");
  m->GetYaxis()->SetTitle("p_{T,true} (GeV/c)");
  m->GetZaxis()->SetTitle(Form("%s", what));
  m->GetYaxis()->SetMoreLogLabels();
  m->GetYaxis()->SetNoExponent();
  m->Draw("COLZ");
  MaskUndefinedBins(e);
  DrawBinCounts(e, 0.022);
  gPad->RedrawAxis();
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.034);
  t.DrawLatex(0.14, 0.955, "#bf{#it{sPHENIX}} Simulation");
  t.SetTextSize(0.028);
  t.DrawLatex(0.14, 0.915, sample.c_str());
  t.SetTextSize(0.040);
  t.DrawLatex(0.14, 0.855, Form("%s:  %s", s.label, FormatIntegrated(Integrate(e)).c_str()));
  t.SetTextSize(0.027);
  t.DrawLatex(0.14, 0.795, Form("text: reconstructed/generated;  hatched: no generated particle (%d of %d bins)", CountUndefinedBins(e),
                                m->GetNbinsX() * m->GetNbinsY()));
  out.Keep(m);
}

// 1D comparison of several efficiencies of this file (e.g. pi+ and pi-) versus
// "pt" or "y"; drawing is done by the shared DrawEfficiencyCurves. Correct-charge
// curves (names starting with "effQ") are dashed with open markers.
void DrawComparisonPanel(const Inputs &in, const std::vector<std::pair<std::string, const SpeciesStyle *>> &curves,
                         const char *proj, const std::string &sample, const Output &out, bool compact = false)
{
  std::vector<Curve> cv;
  for (const auto &c : curves)
  {
    const bool q = c.first.rfind("effQ", 0) == 0;
    const int m = c.second->marker;
    Curve k{in.eff.at(c.first), "", q ? "corr." : "any", c.second->color, m, 1, ""};
    if (q)
    {
      k.marker = m == 20 ? 24 : (m == 21 ? 25 : (m == 22 ? 26 : m));  // open marker
      k.lineStyle = 2;
    }
    k.label = compact ? std::string(c.second->label) + (q ? " correct charge" : " any charge")
                      : std::string(c.second->label) + (q ? " (corr. charge)" : "");
    cv.push_back(k);
  }
  for (TGraphAsymmErrors *g : DrawEfficiencyCurves(cv, std::string(proj) == "pt", sample, compact)) out.Keep(g);
}

// Raw counts versus true pT: denominator, numerator, correct-charge numerator.
void DrawCountsPanel(const Inputs &in, const SpeciesStyle &s, const std::string &sample, const Output &out)
{
  gPad->SetLogx();
  gPad->SetLogy();
  TH1 *d = in.hist.at(Form("hDen_%s_pt", s.tag)), *n = in.hist.at(Form("hNum_%s_pt", s.tag)),
      *q = in.hist.at(Form("hNumQ_%s_pt", s.tag));
  // draw clones (the inputs stay untouched)
  TH1 *dc = static_cast<TH1 *>(d->Clone(Form("draw_%s", d->GetName())));
  TH1 *nc = static_cast<TH1 *>(n->Clone(Form("draw_%s", n->GetName())));
  TH1 *qc = static_cast<TH1 *>(q->Clone(Form("draw_%s", q->GetName())));
  for (TH1 *h : {dc, nc, qc}) h->SetBit(kCanDelete);
  dc->SetLineColor(kBlack);
  dc->SetLineWidth(2);
  nc->SetLineColor(s.color);
  nc->SetLineWidth(2);
  qc->SetLineColor(s.color);
  qc->SetLineStyle(2);
  qc->SetLineWidth(2);
  dc->SetMinimum(0.5);
  dc->SetMaximum(std::max(10., dc->GetMaximum() * 60));  // log axis: room above the data for the legend
  dc->GetXaxis()->SetTitle("p_{T,true} (GeV/c)");
  dc->GetYaxis()->SetTitle("particles per bin");
  dc->GetXaxis()->SetMoreLogLabels();
  dc->GetXaxis()->SetNoExponent();
  dc->Draw("HIST");
  nc->Draw("HIST SAME");
  qc->Draw("HIST SAME");
  auto *leg = new TLegend(0.45, 0.70, 0.93, 0.88);
  leg->SetBit(kCanDelete);
  leg->SetTextSize(0.036);
  leg->SetFillStyle(0);
  leg->SetHeader(s.label);
  leg->AddEntry(dc, Form("generated = %.0f", d->Integral()), "l");
  leg->AddEntry(nc, Form("reconstructed = %.0f", n->Integral()), "l");
  leg->AddEntry(qc, Form("reco., correct charge = %.0f", q->Integral()), "l");
  leg->Draw();
  DrawHeader(sample);
  out.Keep(dc);
}

// ---- pages ---------------------------------------------------------------------------
TCanvas *NewPage(Output &out, const char *name, int nx, int ny)
{
  auto *c = new TCanvas(name, name, 1600, 900);
  c->Divide(nx, ny, 0.002, 0.002);
  out.page++;
  return c;
}
void ClosePage(Output &out, TCanvas *c)
{
  if (out.pagesPass) c->Print(out.pdf.c_str(), Form("Title:%s", c->GetName()));
  out.SaveFigure(c, Form("page%02d_%s", out.page, c->GetName()));
  delete c;
}

void SummaryPage(Output &out, const Inputs &in, const std::string &sample)
{
  auto *c = new TCanvas("summary", "summary", 1600, 900);
  out.page++;
  auto *p = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
  p->SetBit(kCanDelete);
  p->SetFillColor(0);
  p->SetBorderSize(0);
  p->SetTextAlign(12);
  p->SetTextFont(42);
  p->SetTextSize(0.023);
  p->AddText("#bf{Single-particle tracking efficiency}  (#bf{#it{sPHENIX}} Simulation)");
  p->AddText(("Sample: " + sample).c_str());
  p->AddText(Form("Events: %.0f, selected %.0f.  Selected true particles %.0f, reconstructed %.0f (correct charge %.0f).",
                  Stat(in, "events"), Stat(in, "selected events"), Stat(in, "selected true particles"),
                  Stat(in, "found (any charge)"), Stat(in, "found (correct charge)")));
  const auto cfg = lowpt::ParseConfig(in.config.c_str());
  p->AddText(Form("Events: %s < %s cm%s.  True particles: G4 primaries, p_{T} > %s GeV/c, %s.",
                  cfg.at("use_reco_zvtx") == "1" ? "|z_{reco}|" : "|z_{true}|", cfg.at("zvtx_max").c_str(),
                  cfg.at("use_reco_zvtx") == "1" ? " (reconstructed vertex required)" : " (all collisions)",
                  cfg.at("truth_pt_min").c_str(),
                  trkeffcfg::EffectiveTruthEtaMax({cfg.begin(), cfg.end()}) > 0 ? Form("|#eta| < %g", trkeffcfg::EffectiveTruthEtaMax({cfg.begin(), cfg.end()})) : "no #eta cut"));
  p->AddText(Form("Tracks: p_{T} > %s GeV/c, |#eta| < %s, n_{TPC} #geq %s, n_{MVTX} #geq %s, no DCA cut; matched to the true particle "
                  "(most clusters), each true particle counted once.",
                  cfg.at("trk_pt_min").c_str(), cfg.at("trk_abseta_max").c_str(), cfg.at("trk_ntpc_min").c_str(),
                  cfg.at("trk_nmvtx_min").c_str()));
  p->AddText("#varepsilon = reconstructed / generated (any charge);  #varepsilon_{Q} = reconstructed with correct charge / generated.  "
             "Binned in TRUE y and TRUE p_{T}.");
  p->AddText("Errors: Clopper-Pearson 68.3%. Bins without generated particles are undefined (hatched in maps, no point in graphs).");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText(" ");
  p->AddText("Low statistics: most (y, p_{T}) bins of kaons and (anti)protons are empty or contain a few particles; "
             "fluctuations are shown, not smoothed.");
  p->Draw();
  // species table with fixed column positions
  const double col[7] = {0.06, 0.15, 0.22, 0.29, 0.37, 0.60, 0.83};
  const char *head[7] = {"species", "generated", "reco.", "reco. (q ok)", "#varepsilon = reco./gen. (integrated)", "#varepsilon_{Q} = reco. (q ok)/gen.",
                         "empty bins (y,p_{T}) / p_{T} / y"};
  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.024);
  double yrow = 0.44;
  for (int k = 0; k < 7; ++k) t.DrawLatex(col[k], yrow, head[k]);
  t.SetTextFont(42);
  for (const auto &s : Species())
  {
    yrow -= 0.050;
    const TEfficiency *e = in.eff.at(Form("eff_%s_ypt", s.tag));
    const Integrated a = Integrate(e), q = Integrate(in.eff.at(Form("effQ_%s_ypt", s.tag)));
    t.DrawLatex(col[0], yrow, s.label);
    t.DrawLatex(col[1], yrow, Form("%.0f", a.total));
    t.DrawLatex(col[2], yrow, Form("%.0f", a.passed));
    t.DrawLatex(col[3], yrow, Form("%.0f", q.passed));
    t.DrawLatex(col[4], yrow, FormatIntegrated(a).c_str());
    t.DrawLatex(col[5], yrow, FormatIntegrated(q).c_str());
    t.DrawLatex(col[6], yrow, Form("%d / %d / %d  (of %d / %d / %d)", CountUndefinedBins(e),
                                   CountUndefinedBins(in.eff.at(Form("eff_%s_pt", s.tag))),
                                   CountUndefinedBins(in.eff.at(Form("eff_%s_y", s.tag))),
                                   e->GetTotalHistogram()->GetNbinsX() * e->GetTotalHistogram()->GetNbinsY(),
                                   e->GetTotalHistogram()->GetNbinsY(), e->GetTotalHistogram()->GetNbinsX()));
  }
  if (out.pagesPass) c->Print(out.pdf.c_str(), "Title:summary");
  out.SaveFigure(c, Form("page%02d_summary", out.page));
  delete c;
}

void MapPages(Output &out, const Inputs &in, const std::string &sample, const char *kind, const char *what)
{
  const auto &S = Species();
  for (size_t i = 0; i < S.size(); i += 2)
  {
    TCanvas *c = NewPage(out, Form("%s_maps_%s_%s", kind, S[i].tag, S[i + 1].tag), 2, 1);
    for (int k = 0; k < 2; ++k)
    {
      c->cd(k + 1);
      DrawMapPanel(in.eff.at(Form("%s_%s_ypt", kind, S[i + k].tag)), S[i + k], what, sample, out);
    }
    ClosePage(out, c);
    // single-panel figures
    for (int k = 0; k < 2; ++k)
    {
      TCanvas one("one", "", 900, 750);
      DrawMapPanel(in.eff.at(Form("%s_%s_ypt", kind, S[i + k].tag)), S[i + k], what, sample, out);
      out.SaveFigure(&one, Form("%s_map_%s", kind, S[i + k].tag));
    }
  }
}

void ComparisonPage(Output &out, const Inputs &in, const std::string &sample, const char *proj)
{
  const auto &S = Species();
  TCanvas *c = NewPage(out, Form("eff_vs_%s", proj), 2, 2);
  for (int pair = 0; pair < 3; ++pair)
  {
    c->cd(pair + 1);
    std::vector<std::pair<std::string, const SpeciesStyle *>> curves = {
        {Form("eff_%s_%s", S[2 * pair].tag, proj), &S[2 * pair]}, {Form("eff_%s_%s", S[2 * pair + 1].tag, proj), &S[2 * pair + 1]}};
    DrawComparisonPanel(in, curves, proj, sample, out);
    TCanvas one("one", "", 900, 750);
    DrawComparisonPanel(in, curves, proj, sample, out);
    out.SaveFigure(&one, Form("eff_vs_%s_%s_%s", proj, S[2 * pair].tag, S[2 * pair + 1].tag));
  }
  c->cd(4);
  std::vector<std::pair<std::string, const SpeciesStyle *>> all;
  for (const auto &s : S) all.push_back({Form("eff_%s_%s", s.tag, proj), &s});
  DrawComparisonPanel(in, all, proj, sample, out);
  {
    TCanvas one("one", "", 900, 750);
    DrawComparisonPanel(in, all, proj, sample, out);
    out.SaveFigure(&one, Form("eff_vs_%s_all", proj));
  }
  ClosePage(out, c);
}

void ChargePage(Output &out, const Inputs &in, const std::string &sample)
{
  const auto &S = Species();
  TCanvas *c = NewPage(out, "any_vs_correct_charge", 3, 2);
  // layout: top row positive, bottom row negative (pi, K, p columns)
  const int padOf[6] = {1, 4, 2, 5, 3, 6};
  for (size_t i = 0; i < S.size(); ++i)
  {
    c->cd(padOf[i]);
    std::vector<std::pair<std::string, const SpeciesStyle *>> curves = {{Form("eff_%s_pt", S[i].tag), &S[i]},
                                                                        {Form("effQ_%s_pt", S[i].tag), &S[i]}};
    DrawComparisonPanel(in, curves, "pt", sample, out, true);
    TCanvas one("one", "", 900, 750);
    DrawComparisonPanel(in, curves, "pt", sample, out, true);
    out.SaveFigure(&one, Form("any_vs_correct_charge_%s", S[i].tag));
  }
  ClosePage(out, c);
}

void CountsPage(Output &out, const Inputs &in, const std::string &sample)
{
  const auto &S = Species();
  TCanvas *c = NewPage(out, "counts_vs_pt", 3, 2);
  const int padOf[6] = {1, 4, 2, 5, 3, 6};
  for (size_t i = 0; i < S.size(); ++i)
  {
    c->cd(padOf[i]);
    DrawCountsPanel(in, S[i], sample, out);
    TCanvas one("one", "", 900, 750);
    DrawCountsPanel(in, S[i], sample, out);
    out.SaveFigure(&one, Form("counts_vs_pt_%s", S[i].tag));
  }
  ClosePage(out, c);
}

// Low-pT diagnostics: fraction of each loss category per true-pT bin, per group.
void LowPtPage(Output &out, const lowpt::Result &R, const std::string &sample)
{
  const int nShow = std::min<int>(7, static_cast<int>(R.cuts.ptEdges.size()) - 1);  // bins up to ~0.6 GeV/c
  const int colors[8] = {kGreen + 2, kRed + 1, kMagenta + 1, kOrange + 7, kYellow + 2, kAzure + 2, kCyan + 1, kGray + 2};
  TCanvas *c = NewPage(out, "lowpt_diagnostics", 2, 2);
  for (size_t g = 0; g < lowpt::kGroups.size(); ++g)
  {
    c->cd(g + 1);
    gPad->SetRightMargin(0.36);
    gPad->SetTopMargin(0.14);
    auto *stack = new THStack(Form("stack_%s", lowpt::kGroups[g].name), "");
    stack->SetBit(kCanDelete);
    auto *leg = new TLegend(0.645, 0.15, 0.995, 0.90);
    leg->SetBit(kCanDelete);
    leg->SetTextSize(0.033);
    leg->SetHeader(Form("%s: category fraction", lowpt::kGroups[g].label));
    std::vector<double> totals(nShow, 0);
    for (int b = 0; b < nShow; ++b)
      for (long n : R.counts[g][b]) totals[b] += n;
    for (size_t k = 0; k < lowpt::kCategories.size(); ++k)
    {
      if (k == 4 && R.cuts.nmvtxMin <= 0) continue;  // nMVTX category impossible without an MVTX cut
      auto *h = new TH1D(Form("lowpt_%s_%zu", lowpt::kGroups[g].name, k), "", nShow, 0, nShow);
      h->SetDirectory(nullptr);
      h->SetBit(kCanDelete);
      for (int b = 0; b < nShow; ++b)
      {
        h->SetBinContent(b + 1, totals[b] > 0 ? R.counts[g][b][k] / totals[b] : 0.);
        h->GetXaxis()->SetBinLabel(b + 1, Form("%.2f-%.2f", R.cuts.ptEdges[b], R.cuts.ptEdges[b + 1]));
      }
      h->SetFillColor(colors[k]);
      h->SetLineColor(colors[k]);
      stack->Add(h);
      leg->AddEntry(h, lowpt::kCategories[k].c_str(), "f");
      out.Keep(h);
    }
    stack->Draw("HIST");
    stack->GetXaxis()->SetTitle("p_{T,true} bin (GeV/c)");
    stack->GetYaxis()->SetTitle("fraction of generated particles");
    stack->GetXaxis()->SetLabelSize(0.040);
    stack->SetMaximum(1.0);
    TLatex t;
    t.SetTextSize(0.030);
    t.SetTextAlign(21);
    for (int b = 0; b < nShow; ++b) t.DrawLatex(b + 0.5, 1.01, Form("gen. %.0f", totals[b]));
    leg->Draw();
    DrawHeader(sample, 0.14, 0.935);
  }
  // 4th pad: reconstructed pT of the best matched track for true pT in the lowest bin (pions)
  c->cd(4);
  auto *h = static_cast<TH1D *>(R.recoPtLowestBin[0]->Clone());  // the result is reused in both passes
  h->SetBit(kCanDelete);
  h->SetLineColor(kAzure + 2);
  h->SetLineWidth(2);
  h->Draw("HIST");
  TLine cut(R.cuts.trkPtMin, 0, R.cuts.trkPtMin, h->GetMaximum() * 1.05);
  cut.SetLineStyle(2);
  cut.SetLineColor(kRed + 1);
  cut.Draw();
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.036);
  const double below = h->Integral(0, h->FindBin(R.cuts.trkPtMin - 1e-6));
  t.DrawLatex(0.45, 0.82, Form("matched pions: %.0f; reco p_{T} #leq %.2f: %.0f", h->GetEntries(), R.cuts.trkPtMin, below));
  t.DrawLatex(0.45, 0.76, "dashed: track p_{T} cut");
  DrawHeader(sample);
  out.Keep(h);
  ClosePage(out, c);
}

void WriteNumbers(const Output &out, const Inputs &in)
{
  FILE *f = fopen((out.dir + "/TrackEfficiency_numbers.txt").c_str(), "w");
  fprintf(f, "# species  generated  reconstructed  reconstructed_correct_charge  eff  eff_lo  eff_up  effQ  effQ_lo  effQ_up  empty_ypt  empty_pt  empty_y  nbins_ypt\n");
  for (const auto &s : Species())
  {
    const TEfficiency *e = in.eff.at(Form("eff_%s_ypt", s.tag)), *q = in.eff.at(Form("effQ_%s_ypt", s.tag));
    const Integrated a = Integrate(e), b = Integrate(q);
    fprintf(f, "%-5s %6.0f %6.0f %6.0f  %.4f %.4f %.4f  %.4f %.4f %.4f  %d %d %d %d\n", s.tag, a.total, a.passed, b.passed, a.eff,
            a.low, a.up, b.eff, b.low, b.up, CountUndefinedBins(e), CountUndefinedBins(in.eff.at(Form("eff_%s_pt", s.tag))),
            CountUndefinedBins(in.eff.at(Form("eff_%s_y", s.tag))),
            e->GetTotalHistogram()->GetNbinsX() * e->GetTotalHistogram()->GetNbinsY());
  }
  fclose(f);
}
}  // namespace plottrkeff

void PlotTrkEff(const char *trkEffFile, const char *outDir, const char *sampleLabel = "",
                const char *pdfName = "TrackEfficiency_K1.pdf", const char *skimList = "")
{
  using namespace plottrkeff;
  SetTrkEffStyle();
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);  // plotted objects are written explicitly (Output::Keep), never auto-attached
  const Inputs in = ReadInputs(trkEffFile);

  Output out;
  out.dir = outDir;
  gSystem->mkdir((out.dir + "/figures").c_str(), true);
  out.pdf = out.dir + "/" + pdfName;
  out.objects = TFile::Open((out.dir + "/TrackEfficiency_plots.root").c_str(), "RECREATE");
  const std::string sample = sampleLabel;

  lowpt::Result lowPt;
  const bool haveSkims = skimList && *skimList;
  if (haveSkims)
  {
    lowPt = lowpt::Run(trkEffFile, skimList);
    FILE *txt = fopen((out.dir + "/LowPtDiagnostics.txt").c_str(), "w");
    lowpt::WriteTable(lowPt, txt);
    fclose(txt);
  }
  for (const bool pagesPass : {true, false})
  {
    out.pagesPass = pagesPass;
    out.page = 0;
    if (pagesPass)
    {
      TCanvas c0("c0", "", 1600, 900);
      c0.Print((out.pdf + "[").c_str());  // open the multipage PDF
    }
    SummaryPage(out, in, sample);                                               // page 1
    MapPages(out, in, sample, "eff", "#varepsilon (any charge)");              // pages 2-4
    ComparisonPage(out, in, sample, "pt");                                      // page 5
    ComparisonPage(out, in, sample, "y");                                       // page 6
    ChargePage(out, in, sample);                                                // page 7
    CountsPage(out, in, sample);                                                // page 8
    MapPages(out, in, sample, "effQ", "#varepsilon_{Q} (correct charge)");     // pages 9-11
    if (haveSkims) LowPtPage(out, lowPt, sample);                               // page 12
    if (pagesPass)
    {
      TCanvas c0("c0", "", 1600, 900);
      c0.Print((out.pdf + "]").c_str());  // close the multipage PDF
    }
  }
  WriteNumbers(out, in);
  out.objects->Close();
  printf("PlotTrkEff: %d pages -> %s\n", out.page, out.pdf.c_str());
}
