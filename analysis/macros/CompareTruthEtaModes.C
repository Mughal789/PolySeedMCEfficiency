// =============================================================================
// CompareTruthEtaModes.C - the two tracking-efficiency definitions side by side
//   A "inclusive"  : use_truth_eta_cut = 0 - all generated primaries in the
//                    (y, pT) bin; the reconstructed |eta| < 1.1 acceptance is
//                    part of the efficiency (FINAL definition)
//   B "restricted" : use_truth_eta_cut = 1 - generated primaries with
//                    |eta_true| < 1.1 only; efficiency conditional on the
//                    particle being inside the eta acceptance
// Read-only; every number comes from the stored raw counts / TEfficiency.
// Pages: 1 integrated D, N, eps (all in-range bins and the species analysis
//          windows); 2 eps vs true pT; 3 eps vs true y; 4 acceptance boundary:
//          fraction of generated particles inside |eta_true| < 1.1 (D_B/D_A)
//          and eps_B - eps_A vs true y.
// No smoothing or interpolation; empty bins have no point.
//
// Usage: root -l -b -q 'CompareTruthEtaModes.C+("A.root","B.root","OUTDIR","sample label")'
//   -> OUTDIR/TruthEtaModes.pdf, OUTDIR/TruthEtaModes.txt, OUTDIR/figures/
// =============================================================================
#include "TrkEffPlotStyle.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace cmpeta
{
FILE *gTxt = nullptr;
void Txt(const std::string &s)
{
  printf("%s\n", s.c_str());
  if (gTxt) fprintf(gTxt, "%s\n", s.c_str());
}
struct Sp
{
  const char *tag, *label;
  int color;
  double yw, plo, phi;  // analysis window
};
const Sp kSp[6] = {{"pip", "#pi^{+}", kAzure + 2, 1.0, 0.1, 99}, {"pim", "#pi^{-}", kAzure + 2, 1.0, 0.1, 99},
                   {"kp", "K^{+}", kRed + 1, 0.7, 0.1, 0.6},      {"km", "K^{-}", kRed + 1, 0.7, 0.1, 0.6},
                   {"p", "p", kGreen + 2, 0.6, 0.4, 1.0},          {"pbar", "#bar{p}", kGreen + 2, 0.6, 0.4, 1.0}};
struct Cnt
{
  double d = 0, n = 0, q = 0;
};
Cnt Count(TFile &f, const char *sp, double ylo, double yhi, double plo, double phi)
{
  Cnt c;
  const TH2D *D = f.Get<TH2D>(Form("hDen_%s_ypt", sp)), *N = f.Get<TH2D>(Form("hNum_%s_ypt", sp)), *Q = f.Get<TH2D>(Form("hNumQ_%s_ypt", sp));
  for (int i = 1; i <= D->GetNbinsX(); ++i)
    for (int j = 1; j <= D->GetNbinsY(); ++j)
    {
      const double y = D->GetXaxis()->GetBinCenter(i), p = D->GetYaxis()->GetBinCenter(j);
      if (y > ylo && y < yhi && p > plo && p < phi)
      {
        c.d += D->GetBinContent(i, j);
        c.n += N->GetBinContent(i, j);
        c.q += Q->GetBinContent(i, j);
      }
    }
  return c;
}
std::string Eff(const Cnt &c)
{
  if (c.d <= 0) return "undefined";
  const double lo = TEfficiency::ClopperPearson(c.d, c.n, 0.683, false), up = TEfficiency::ClopperPearson(c.d, c.n, 0.683, true);
  return Form("%.3f [%.3f,%.3f]", c.n / c.d, lo, up);
}
}  // namespace cmpeta

void CompareTruthEtaModes(const char *pathA, const char *pathB, const char *outDir, const char *sampleLabel)
{
  using namespace cmpeta;
  using namespace trkeffplot;
  const std::string dir = outDir;
  gSystem->mkdir((dir + "/figures").c_str(), true);
  gTxt = fopen((dir + "/TruthEtaModes.txt").c_str(), "w");
  TH1::AddDirectory(false);
  std::unique_ptr<TFile> A(TFile::Open(pathA, "READ")), B(TFile::Open(pathB, "READ"));
  SetTrkEffStyle();
  const std::string selA = A->Get<TNamed>("selection")->GetTitle(), selB = B->Get<TNamed>("selection")->GetTitle();
  Txt(std::string("A (inclusive, use_truth_eta_cut = 0): ") + pathA + "\n    " + selA);
  Txt(std::string("B (restricted, use_truth_eta_cut = 1): ") + pathB + "\n    " + selB);
  std::vector<TCanvas *> pages;

  // ---- page 1: integrated table ----------------------------------------------------------
  Txt("\n== integrated, Clopper-Pearson 68.3%");
  Txt("   species  region                     Gen_A   Reco_A    eps_A                  Gen_B   Reco_B    eps_B                  Gen_B/Gen_A");
  std::vector<std::string> lines;
  for (const Sp &sp : kSp)
    for (int w = 0; w < 2; ++w)
    {
      const double ylo = w ? -sp.yw : -99, yhi = w ? sp.yw : 99, plo = w ? sp.plo : -1, phi = w ? sp.phi : 99;
      const Cnt a = Count(*A, sp.tag, ylo, yhi, plo, phi), b = Count(*B, sp.tag, ylo, yhi, plo, phi);
      const std::string region = w ? Form("window |y|<%.1f, %.1f<=pT<%s", sp.yw, sp.plo, sp.phi > 50 ? "5" : Form("%.1f", sp.phi))
                                   : "all in-range bins (|y|<1.2)";
      const std::string l = Form("   %-6s   %-28s %7.0f  %7.0f  %-22s %7.0f  %7.0f  %-22s %.3f", sp.tag, region.c_str(), a.d, a.n,
                                 Eff(a).c_str(), b.d, b.n, Eff(b).c_str(), a.d > 0 ? b.d / a.d : 0.);
      Txt(l);
      lines.push_back(l);
    }
  {
    auto *c = new TCanvas("page1_integrated", "integrated", 1400, 850);
    pages.push_back(c);
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.027);
    double y = 0.95;
    t.DrawLatex(0.03, y, Form("#bf{#it{sPHENIX}} Simulation   %s   -   two tracking-efficiency definitions", sampleLabel));
    y -= 0.045;
    t.SetTextSize(0.020);
    t.DrawLatex(0.03, y, ("A (inclusive, use_truth_eta_cut = 0):  " + selA).c_str());
    y -= 0.032;
    t.DrawLatex(0.03, y, ("B (restricted, use_truth_eta_cut = 1): " + selB).c_str());
    y -= 0.045;
    t.SetTextFont(82);
    t.SetTextSize(0.0165);
    t.DrawLatex(0.03, y, "species  region                       Gen_A   Reco_A    eps_A [CP 68%]         Gen_B   Reco_B    eps_B [CP 68%]        Gen_B/Gen_A");
    y -= 0.032;
    for (size_t i = 0; i < lines.size(); ++i)
    {
      t.DrawLatex(0.03, y, lines[i].substr(3).c_str());
      y -= (i % 2) ? 0.042 : 0.030;
    }
    t.SetTextFont(42);
    t.SetTextSize(0.019);
    t.DrawLatex(0.03, 0.10, "A: generated = all primaries in the (y, p_{T}) bin; particles leaving the |#eta| < 1.1 acceptance count as inefficiency.");
    t.DrawLatex(0.03, 0.07, "B: generated = primaries with |#eta_{true}| < 1.1 only; efficiency conditional on the particle being inside the #eta acceptance.");
    t.DrawLatex(0.03, 0.04, "Both: same events, track cuts (incl. reconstructed |#eta| < 1.1) and binning; reconstructed #subset generated. #varepsilon_{A} = #varepsilon_{B} #times Gen_{B}/Gen_{A} + (reconstructed with |#eta_{true}| #geq 1.1)/Gen_{A}.");
  }

  // ---- pages 2 and 3: eps vs pT and vs y ---------------------------------------------------
  for (int isPt = 1; isPt >= 0; --isPt)
  {
    auto *c = new TCanvas(isPt ? "page2_pt" : "page3_y", "", 1700, 950);
    pages.push_back(c);
    c->Divide(3, 2);
    for (int s = 0; s < 6; ++s)
    {
      c->cd((s % 2) * 3 + s / 2 + 1);
      const TEfficiency *ea = A->Get<TEfficiency>(Form("eff_%s_%s", kSp[s].tag, isPt ? "pt" : "y"));
      const TEfficiency *eb = B->Get<TEfficiency>(Form("eff_%s_%s", kSp[s].tag, isPt ? "pt" : "y"));
      std::vector<Curve> curves = {{ea, std::string(kSp[s].label) + " A incl.", "A", kGray + 2, 24, 1, ""},
                                   {eb, std::string(kSp[s].label) + " B |#eta_{true}|<1.1", "B", kSp[s].color, 20, 1, ""}};
      DrawEfficiencyCurves(curves, isPt, sampleLabel, true);
      DrawHeader(sampleLabel);
    }
  }

  // ---- page 4: acceptance boundary ----------------------------------------------------------------
  Txt("\n== acceptance boundary vs true y (pT summed over the in-range bins): Gen_B/Gen_A and eps_B - eps_A");
  {
    auto *c = new TCanvas("page4_boundary", "boundary", 1700, 950);
    pages.push_back(c);
    c->Divide(3, 2);
    for (int s = 0; s < 6; ++s)
    {
      c->cd((s % 2) * 3 + s / 2 + 1);
      gPad->SetRightMargin(0.04);
      const TH1D *dA = A->Get<TH1D>(Form("hDen_%s_y", kSp[s].tag)), *dB = B->Get<TH1D>(Form("hDen_%s_y", kSp[s].tag));
      const TH1D *nA = A->Get<TH1D>(Form("hNum_%s_y", kSp[s].tag)), *nB = B->Get<TH1D>(Form("hNum_%s_y", kSp[s].tag));
      const double lo = dA->GetXaxis()->GetXmin(), hi = dA->GetXaxis()->GetXmax();
      auto *f = new TH1D(Form("fb_%s", kSp[s].tag), ";y_{true};", 10, lo, hi);
      f->SetBit(kCanDelete);
      f->SetMinimum(-0.2);
      f->SetMaximum(1.6);  // free band above 1.2 for the legend
      f->Draw("AXIS");
      TLine l;
      l.SetLineStyle(2);
      l.SetLineColor(kGray + 1);
      l.DrawLine(lo, 0, hi, 0);
      l.DrawLine(lo, 1, hi, 1);
      auto *gFrac = new TGraphAsymmErrors();
      auto *gDiff = new TGraphErrors();
      gFrac->SetBit(kCanDelete);
      gDiff->SetBit(kCanDelete);
      std::string row = Form("   %-5s", kSp[s].tag);
      for (int b = 1; b <= dA->GetNbinsX(); ++b)
      {
        const double DA = dA->GetBinContent(b), DB = dB->GetBinContent(b), x = dA->GetXaxis()->GetBinCenter(b), hw = 0.5 * dA->GetXaxis()->GetBinWidth(b);
        if (DA <= 0) continue;
        const double fr = DB / DA, flo = TEfficiency::ClopperPearson(DA, DB, 0.683, false), fup = TEfficiency::ClopperPearson(DA, DB, 0.683, true);
        gFrac->AddPoint(x, fr);
        gFrac->SetPointError(gFrac->GetN() - 1, hw, hw, fr - flo, fup - fr);
        std::string d = "    -   ";
        if (DB > 0)
        {
          // B is a subset of A (correlated): show the difference with the error of the larger-error efficiency (A)
          const double eA = nA->GetBinContent(b) / DA, eB = nB->GetBinContent(b) / DB;
          gDiff->AddPoint(x, eB - eA);
          gDiff->SetPointError(gDiff->GetN() - 1, hw, std::sqrt(std::max(eA * (1 - eA), 1.0 / DA) / DA));
          d = Form("%+.3f", eB - eA);
        }
        row += Form(" [%+.1f] %.3f %s |", dA->GetXaxis()->GetBinLowEdge(b), fr, d.c_str());
      }
      Txt(row);
      gFrac->SetMarkerStyle(21);
      gFrac->SetMarkerColor(kBlack);
      gFrac->SetLineColor(kBlack);
      gFrac->SetLineWidth(2);
      gDiff->SetMarkerStyle(20);
      gDiff->SetMarkerColor(kSp[s].color);
      gDiff->SetLineColor(kSp[s].color);
      gDiff->SetLineWidth(2);
      gFrac->Draw("P SAME");
      gDiff->Draw("P SAME");
      auto *leg = new TLegend(0.18, 0.74, 0.95, 0.90);
      leg->SetFillStyle(0);
      leg->SetBit(kCanDelete);
      leg->SetTextSize(0.042);
      leg->SetHeader(kSp[s].label);
      leg->AddEntry(gFrac, "fraction of generated inside |#eta_{true}|<1.1", "pe");
      leg->AddEntry(gDiff, "#varepsilon_{B} - #varepsilon_{A}", "pe");
      leg->Draw();
      DrawHeader(sampleLabel);
    }
  }

  const std::string pdf = dir + "/TruthEtaModes.pdf";
  pages.front()->Print((pdf + "[").c_str());
  for (auto *c : pages) c->Print(pdf.c_str(), Form("Title:%s", c->GetName()));
  pages.front()->Print((pdf + "]").c_str());
  for (auto *c : pages)
  {
    c->SaveAs((dir + "/figures/" + c->GetName() + ".pdf").c_str());
    c->SaveAs((dir + "/figures/" + c->GetName() + ".png").c_str());
  }
  Txt("wrote " + pdf);
  if (gTxt) fclose(gTxt);
}
