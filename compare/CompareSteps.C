// =============================================================================
// CompareSteps.C
//
// Compares several production steps using their qa/EvalQA_hists.root
// (written by PlotEvalQA.C).
//
// Usage (compare/compare.sh calls it):
//   root -l -b -q 'CompareSteps.C+("stepA/qa/EvalQA_hists.root,stepB/qa/EvalQA_hists.root",
//                                           "labelA,labelB","compare/outdir")'
// Writes into outDir:
//   CompareSteps.pdf, png/      overlays of every 1D histogram and efficiency
//   CompareSteps_table.txt      side-by-side table of the key numbers
// Steps produced with the same seeds contain the same events, so raw counts
// are compared directly (no normalisation).
// =============================================================================
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::vector<std::string> Split(const std::string &s)
{
  std::vector<std::string> v;
  std::stringstream ss(s);
  std::string x;
  while (std::getline(ss, x, ',')) v.push_back(x);
  return v;
}
const int kCol[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1};

double Mode(TH1 *h) { return h->GetXaxis()->GetBinCenter(h->GetMaximumBin()); }
double FracIn(TH1 *h, double lo, double hi)
{
  const double all = h->Integral(0, h->GetNbinsX() + 1);
  return all > 0 ? h->Integral(h->FindBin(lo + 1e-9), h->FindBin(hi - 1e-9)) / all : 0;
}
}  // namespace

void CompareSteps(const char *files, const char *labels, const char *outDir)
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TH1::AddDirectory(false);
  gSystem->mkdir(Form("%s/png", outDir), true);
  const auto fn = Split(files), lb = Split(labels);
  if (fn.size() != lb.size() || fn.empty())
  {
    std::cerr << "FATAL: need the same number of files and labels" << std::endl;
    return;
  }
  std::vector<TFile *> F;
  for (const auto &f : fn)
  {
    F.push_back(TFile::Open(f.c_str(), "READ"));
    if (!F.back() || F.back()->IsZombie())
    {
      std::cerr << "FATAL: cannot open " << f << std::endl;
      return;
    }
  }
  const size_t n = F.size();
  auto get = [&](size_t i, const char *name) { return dynamic_cast<TH1 *>(F[i]->Get(name)); };

  // ---------------------------------------------------------------------------
  // Table of key numbers
  // ---------------------------------------------------------------------------
  std::ofstream t(Form("%s/CompareSteps_table.txt", outDir));
  auto row = [&](const std::string &name, std::function<std::string(size_t)> fnc) {
    t << Form("%-46s", name.c_str());
    for (size_t i = 0; i < n; ++i) t << Form(" %18s", fnc(i).c_str());
    t << "\n";
  };
  t << Form("%-46s", "quantity");
  for (size_t i = 0; i < n; ++i) t << Form(" %18s", lb[i].c_str());
  t << "\n" << std::string(46 + 19 * n, '-') << "\n";
  auto evts = [&](size_t i) { return get(i, "H_n_tpc")->GetEntries(); };
  row("events", [&](size_t i) { return Form("%.0f", evts(i)); });
  row("[A] events with >=1 vertex (fraction)", [&](size_t i) {
    TH1 *h = get(i, "A_nvtx");
    return Form("%.3f", 1. - h->GetBinContent(1) / h->GetEntries());
  });
  row("[A] vertex dz to truth: mean / RMS (um)", [&](size_t i) {
    TH1 *h = get(i, "A_dz");
    return Form("%.0f / %.0f", 1e4 * h->GetMean(), 1e4 * h->GetRMS());
  });
  row("[A] vertex x, y mean (cm)", [&](size_t i) { return Form("%.3f, %.3f", get(i, "A_vx")->GetMean(), get(i, "A_vy")->GetMean()); });
  row("[A] tracks per vertex (mean)", [&](size_t i) { return Form("%.2f", get(i, "A_vtxntrk")->GetMean()); });
  row("[C] tracks per event", [&](size_t i) { return Form("%.2f", get(i, "C_pt")->GetEntries() / evts(i)); });
  row("[C] tracks with silicon seed (fraction)", [&](size_t i) {
    TH1 *h = get(i, "C_silseed");
    return Form("%.3f", h->GetBinContent(2) / h->GetEntries());
  });
  row("[C] good tracks: primary / secondary", [&](size_t i) {
    TH1 *h = get(i, "C_truthcomp");
    return Form("%.0f / %.0f", h->GetBinContent(2), h->GetBinContent(3));
  });
  row("[D] tracks with vertexID < 0 (fraction)", [&](size_t i) {
    TH1 *h = get(i, "D_vertexID");
    return Form("%.3f", (h->GetBinContent(1) + h->GetBinContent(3)) / h->GetEntries());
  });
  row("[D] track vertex == true vertex (fraction)", [&](size_t i) {
    TH1 *h = get(i, "D_vtx_is_truth");
    return Form("%.3f", h->GetBinContent(2) / h->GetEntries());
  });
  for (const char *nm : {"E_pcaz_tpconly_etaneg", "E_pcaz_tpconly_etapos", "E_pcaz_mvtx"})
    row(std::string("[E] ") + (nm + 7) + ": mode/mean (cm)", [&](size_t i) {
      TH1 *h = get(i, nm);
      return h->GetEntries() > 0 ? Form("%.2f / %.2f", Mode(h), h->GetMean()) : std::string("-");
    });
  for (const char *nm : {"E_pcaz_tpconly_etaneg", "E_pcaz_tpconly_etapos", "E_pcaz_mvtx"})
    row(std::string("[E] ") + (nm + 7) + ": |dz|<1 cm", [&](size_t i) {
      TH1 *h = get(i, nm);
      return h->GetEntries() > 0 ? Form("%.3f (N=%.0f)", FracIn(h, -1, 1), h->GetEntries()) : std::string("-");
    });
  for (const char *nm : {"E_tanl_tpconly_etaneg", "E_tanl_tpconly_etapos", "E_tanl_mvtx"})
    row(std::string("[E] ") + (nm + 7) + ": mean", [&](size_t i) {
      TH1 *h = get(i, nm);
      return h->GetEntries() > 0 ? Form("%.3f", h->GetMean()) : std::string("-");
    });
  {
    const char *cat[5] = {"no track", "fails good", "good TPC-only", "good 1 MVTX", "good >=2 MVTX"};
    for (int k = 1; k <= 5; ++k)
      row(std::string("[F] primaries: ") + cat[k - 1], [&](size_t i) {
        auto *h = dynamic_cast<TH2 *>(F[i]->Get("F_fate_vs_pt"));
        double s = 0, c = 0;
        for (int b = 1; b <= h->GetNbinsX(); ++b)
          for (int j = 1; j <= 5; ++j)
          {
            s += h->GetBinContent(b, j);
            if (j == k) c += h->GetBinContent(b, j);
          }
        return Form("%.3f", s > 0 ? c / s : 0.);
      });
  }
  row("[F] crossed >=2 MVTX layers -> reco >=2 MVTX", [&](size_t i) {
    TH1 *h = get(i, "F_prim_truemvtx_good");
    const double all = h->GetEntries();
    double m2 = 0;
    for (int b = h->FindBin(2.); b <= h->GetNbinsX(); ++b) m2 += h->GetBinContent(b);
    return Form("%.3f (N=%.0f)", all > 0 ? m2 / all : 0., all);
  });
  row("[H] mean n_TPC", [&](size_t i) { return Form("%.3f", get(i, "H_n_tpc")->GetMean()); });
  row("[H] mean n_MVTX2", [&](size_t i) { return Form("%.3f", get(i, "H_n_mvtx2")->GetMean()); });
  row("[H] mean N (true)", [&](size_t i) { return Form("%.3f", get(i, "H_N_true")->GetMean()); });
  t.close();
  std::cout << std::ifstream(Form("%s/CompareSteps_table.txt", outDir)).rdbuf() << std::endl;

  // ---------------------------------------------------------------------------
  // Overlays: every 1D histogram and every efficiency of the first file
  // ---------------------------------------------------------------------------
  std::vector<std::string> h1, eff;
  for (TObject *o : *F[0]->GetListOfKeys())
  {
    TKey *k = static_cast<TKey *>(o);
    TClass *cl = TClass::GetClass(k->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom(TH2::Class())) continue;
    if (cl->InheritsFrom(TH1::Class())) h1.push_back(k->GetName());
    if (cl->InheritsFrom(TEfficiency::Class())) eff.push_back(k->GetName());
  }
  std::vector<TCanvas *> pages;
  auto newPage = [&](int idx) {
    auto *c = new TCanvas(Form("cmp%02d", idx), "", 1100, 850);
    c->Divide(3, 2, 0.004, 0.004);
    pages.push_back(c);
    return c;
  };
  auto legend = [&](std::vector<TObject *> objs) {
    auto *lg = new TLegend(0.50, 0.72, 0.93, 0.90);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetTextSize(0.035);
    for (size_t i = 0; i < objs.size(); ++i)
      if (objs[i]) lg->AddEntry(objs[i], lb[i].c_str(), "l");
    lg->Draw();
  };
  int pad = 6, idx = 0;
  TCanvas *c = nullptr;
  for (const auto &name : h1)
  {
    if (pad == 6)
    {
      c = newPage(idx++);
      pad = 0;
    }
    c->cd(++pad);
    gPad->SetTopMargin(0.09);
    double mx = 0;
    std::vector<TObject *> objs(n, nullptr);
    for (size_t i = 0; i < n; ++i)
      if (TH1 *h = get(i, name.c_str())) mx = std::max(mx, h->GetMaximum());
    for (size_t i = 0; i < n; ++i)
    {
      TH1 *h = get(i, name.c_str());
      if (!h) continue;
      h->SetLineColor(kCol[i % 6]);
      h->SetMarkerColor(kCol[i % 6]);
      h->SetLineWidth(2);
      h->SetLineStyle(i == 0 ? 1 : 1 + int(i));
      h->SetMaximum(mx * 1.35);
      h->SetMinimum(0);
      const bool prof = h->InheritsFrom(TProfile::Class());
      h->Draw(i == 0 ? (prof ? "e" : "hist") : (prof ? "e same" : "hist same"));
      objs[i] = h;
    }
    legend(objs);
  }
  for (const auto &name : eff)
  {
    if (pad == 6)
    {
      c = newPage(idx++);
      pad = 0;
    }
    c->cd(++pad);
    gPad->SetTopMargin(0.09);
    std::vector<TObject *> objs(n, nullptr);
    for (size_t i = 0; i < n; ++i)
    {
      auto *e = dynamic_cast<TEfficiency *>(F[i]->Get(name.c_str()));
      if (!e) continue;
      e->SetLineColor(kCol[i % 6]);
      e->SetMarkerColor(kCol[i % 6]);
      e->SetMarkerStyle(20 + int(i));
      e->Draw(i == 0 ? "AP" : "P same");
      gPad->Update();
      if (i == 0 && e->GetPaintedGraph())
      {
        e->GetPaintedGraph()->SetMinimum(0);
        e->GetPaintedGraph()->SetMaximum(1.3);
      }
      objs[i] = e;
    }
    legend(objs);
  }
  const std::string pdf = std::string(outDir) + "/CompareSteps.pdf";
  for (size_t i = 0; i < pages.size(); ++i)
  {
    pages[i]->Print((pdf + (i == 0 ? "(" : (i + 1 == pages.size() ? ")" : ""))).c_str());
    pages[i]->Print(Form("%s/png/%02zu.png", outDir, i + 1));
  }
  std::cout << "wrote " << pdf << " (" << pages.size() << " pages) and " << outDir << "/CompareSteps_table.txt" << std::endl;
}
