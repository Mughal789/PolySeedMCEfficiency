// =============================================================================
// CompareVtxEff.C
//
// Compares any number of samples analysed with analysis/run_analysis.sh.
// Each sample is given by its analysis directory SAMPLEDIR/analysis, which contains
//   z<all|8>/VtxEff_TPC.root, z<all|8>/breakdown/{VtxEffBreakdown_hists.root,VtxEffBreakdown_summary.txt},
//   qa/EvalQA_hists.root, audit/ResidPhi.txt
// Samples are drawn in the given order (first page: epsilon_vtx per sample, so an
// ordered list of steps gives a "staircase").
// Samples produced with the same seeds contain the same events.
// Errors: Clopper-Pearson 68.3%. pT axes are logarithmic; pT spectra also log y.
//
// Usage: root -l -b -q 'CompareVtxEff.C+("dirA/analysis,dirB/analysis","label A,label B","outDir")'
// Writes into outDir: CompareVtxEff.pdf, png/NN_name.png, CompareVtxEff_table.txt
// =============================================================================
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Sample
{
  std::string dir, label;
  int col, marker;
};

const int kSampleColor[] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kGray + 2, kCyan + 2, kViolet + 1};
const int kSampleMarker[] = {20, 21, 24, 26, 25, 22, 32, 23};

std::vector<std::string> Split(const std::string &s)
{
  std::vector<std::string> v;
  std::stringstream ss(s);
  for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
  return v;
}

std::ofstream gTable;
void Out(const std::string &s)
{
  std::cout << s << std::endl;
  gTable << s << std::endl;
}

// "(k/n)" from the VtxEffBreakdown summary line starting with `start`
bool ReadCount(const std::string &summary, const std::string &start, double &k, double &n)
{
  std::ifstream in(summary);
  std::regex re(R"(\((\d+)/(\d+)\))");
  for (std::string l; std::getline(in, l);)
  {
    if (l.find(start) != 0) continue;
    std::smatch m;
    if (!std::regex_search(l, m, re)) continue;
    k = std::stod(m[1]);
    n = std::stod(m[2]);
    return true;
  }
  std::cerr << "WARNING: '" << start << "' not found in " << summary << std::endl;
  return false;
}

// no-vertex reasons inside the central window: R2..R5 counts
std::vector<double> ReadReasons(const std::string &summary)
{
  std::ifstream in(summary);
  std::vector<double> r;
  bool inBlock = false;
  std::regex re(R"(R[2-5] .*\s(\d+) /\s*(\d+) =)");
  for (std::string l; std::getline(in, l);)
  {
    if (l.find("events WITHOUT a vertex inside") != std::string::npos) { inBlock = true; continue; }
    if (!inBlock) continue;
    std::smatch m;
    if (!std::regex_search(l, m, re)) break;
    r.push_back(std::stod(m[1]));
  }
  return r;
}

double CPlo(double n, double k) { return TEfficiency::ClopperPearson(n, k, 0.683, false); }
double CPup(double n, double k) { return TEfficiency::ClopperPearson(n, k, 0.683, true); }

TEfficiency *GetEff(const std::string &file, const char *name)
{
  auto *f = TFile::Open(file.c_str(), "READ");
  if (!f || f->IsZombie()) { std::cerr << "WARNING: cannot open " << file << std::endl; return nullptr; }
  auto *e = (TEfficiency *) f->Get(name);
  if (!e) { std::cerr << "WARNING: " << name << " not in " << file << std::endl; return nullptr; }
  e = (TEfficiency *) e->Clone(Form("%s_%p", name, (void *) f));
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(0.683);
  return e;
}

TEfficiency *EffFromHists(const std::string &file, const char *pass, const char *tot)
{
  auto *f = TFile::Open(file.c_str(), "READ");
  if (!f || f->IsZombie()) { std::cerr << "WARNING: cannot open " << file << std::endl; return nullptr; }
  auto *p = (TH1 *) f->Get(pass);
  auto *t = (TH1 *) f->Get(tot);
  if (!p || !t) { std::cerr << "WARNING: " << pass << "/" << tot << " not in " << file << std::endl; return nullptr; }
  auto *e = new TEfficiency(*p, *t);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(0.683);
  return e;
}

TLegend *Leg(double x1, double y1, double x2, double y2)
{
  auto *l = new TLegend(x1, y1, x2, y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.033);
  return l;
}

void Header(const char *text)
{
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.032);
  t.DrawLatex(0.12, 0.955, text);
}

template <class T>
void Style(T *o, const Sample &s)
{
  o->SetLineColor(s.col);
  o->SetMarkerColor(s.col);
  o->SetMarkerStyle(s.marker);
  o->SetMarkerSize(1.3);
  o->SetLineWidth(2);
}

void EffOverlay(const std::vector<Sample> &S, const std::vector<TEfficiency *> &E, const char *frame, double xmin, double xmax,
                double ymin, double ymax, bool logx, const char *header, bool legTop)
{
  gPad->SetLogx(logx);
  auto *fr = gPad->DrawFrame(xmin, ymin, xmax, ymax, frame);
  fr->GetXaxis()->SetMoreLogLabels(logx);
  fr->GetXaxis()->SetNoExponent(logx);
  auto *l = legTop ? Leg(0.14, 0.89 - 0.055 * S.size(), 0.6, 0.89) : Leg(0.45, 0.15, 0.89, 0.15 + 0.055 * S.size());
  for (size_t i = 0; i < S.size(); ++i)
  {
    if (!E[i]) continue;
    Style(E[i], S[i]);
    E[i]->Draw("P same");
    l->AddEntry(E[i], S[i].label.c_str(), "pl");
  }
  l->Draw();
  Header(header);
}
}  // namespace

void CompareVtxEff(const char *dirs, const char *labels, const char *outDir)
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetTitleSize(0.045, "XY");
  gStyle->SetLabelSize(0.04, "XY");
  TH1::AddDirectory(false);

  const auto D = Split(dirs), L = Split(labels);
  if (D.empty() || D.size() != L.size() || D.size() > 8)
  {
    std::cerr << "FATAL: need 1-8 analysis directories and the same number of labels" << std::endl;
    return;
  }
  std::vector<Sample> S;
  for (size_t i = 0; i < D.size(); ++i) S.push_back({D[i], L[i], kSampleColor[i], kSampleMarker[i]});
  gSystem->mkdir(Form("%s/png", outDir), true);
  gTable.open(Form("%s/CompareVtxEff_table.txt", outDir));
  std::vector<std::pair<TCanvas *, std::string>> pages;
  auto NewPage = [&](const char *name) {
    auto *c = new TCanvas(Form("c_%s", name), name, 900, 700);
    pages.push_back({c, name});
    return c;
  };

  // 1. epsilon_vtx per sample (all z and |z_true| < 8 cm)
  {
    auto *c = NewPage("eps_per_sample");
    c->SetBottomMargin(0.16);
    auto *fr = gPad->DrawFrame(-0.5, 0, S.size() - 0.5, 1.0, ";;#varepsilon_{vtx}");
    fr->GetXaxis()->SetNdivisions(0);
    TLatex lab;
    lab.SetTextAlign(23);
    lab.SetTextSize(S.size() > 5 ? 0.024 : 0.03);
    auto *l = Leg(0.14, 0.78, 0.6, 0.9);
    Out("epsilon_vtx (k/n), every collision in the denominator");
    Out(Form("%-40s %-24s %-24s", "sample", "all z", "|z_true| < 8 cm"));
    for (size_t i = 0; i < S.size(); ++i)
    {
      lab.DrawLatex(i, -0.03, S[i].label.c_str());
      std::string row = Form("%-40s", S[i].label.c_str());
      for (int zc = 0; zc < 2; ++zc)
      {
        double k = 0, n = 0;
        const std::string summ = S[i].dir + (zc ? "/z8" : "/zall") + "/breakdown/VtxEffBreakdown_summary.txt";
        if (!ReadCount(summ, zc ? "  |z_true| < 8 cm " : "  all collisions", k, n) || n <= 0)
        {
          row += Form(" %-24s", "-");
          continue;
        }
        row += Form(" %-24s", Form("%.3f +%.3f -%.3f (%.0f/%.0f)", k / n, CPup(n, k) - k / n, k / n - CPlo(n, k), k, n));
        auto *g = new TGraphAsymmErrors(1);
        g->SetPoint(0, i, k / n);
        g->SetPointError(0, 0, 0, k / n - CPlo(n, k), CPup(n, k) - k / n);
        g->SetMarkerStyle(zc ? 20 : 21);
        g->SetMarkerSize(1.6);
        g->SetMarkerColor(zc ? kBlue + 1 : kRed + 1);
        g->SetLineColor(zc ? kBlue + 1 : kRed + 1);
        g->SetLineWidth(2);
        g->Draw("P same");
        if (i == 0) l->AddEntry(g, zc ? "|z_{true}| < 8 cm" : "all z", "p");
      }
      Out(row);
    }
    l->Draw();
    Header("Vertex efficiency per sample; every collision in the denominator");
  }

  // 2. MVTX truth-reco cluster distance vs phi (log y)
  {
    auto *c = NewPage("mvtx_residual_vs_phi");
    c->SetLogy();
    gPad->DrawFrame(-180, 2, 180, 3000, ";cluster #phi (deg);median |truth - reco| MVTX cluster distance (#mum)");
    auto *l = Leg(0.14, 0.62 - 0.05 * S.size(), 0.7, 0.62);
    std::regex re(R"(phi\s+(-?\d+)\.\.\s*(-?\d+)\s*:\s*median\s+(\d+))");
    Out("\nMVTX truth-reco distance (median per 20 deg phi bin, um): min - max");
    for (const auto &s : S)
    {
      std::ifstream in(s.dir + "/audit/ResidPhi.txt");
      auto *g = new TGraphAsymmErrors();
      double mn = 1e9, mx = 0;
      for (std::string line; std::getline(in, line);)
      {
        std::smatch m;
        if (!std::regex_search(line, m, re)) continue;
        const double v = std::stod(m[3]);
        g->SetPoint(g->GetN(), 0.5 * (std::stod(m[1]) + std::stod(m[2])), v);
        g->SetPointError(g->GetN() - 1, 10, 10, 0, 0);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      if (g->GetN() == 0) { std::cerr << "WARNING: no ResidPhi.txt for " << s.label << std::endl; continue; }
      Style(g, s);
      g->Draw("P same");
      l->AddEntry(g, s.label.c_str(), "p");
      Out(Form("  %-40s %4.0f - %4.0f", s.label.c_str(), mn, mx));
    }
    l->Draw();
    Header("MVTX layers 0-2: truth vs nearest reco cluster (spatial match < 3 mm)");
  }

  // 3-7. efficiencies
  std::vector<TEfficiency *> eN, ez, em, et, emx;
  for (const auto &s : S)
  {
    eN.push_back(EffFromHists(s.dir + "/z8/VtxEff_TPC.root", "hNW_vtx_zc", "hNW_all_zc"));
    ez.push_back(GetEff(s.dir + "/zall/breakdown/VtxEffBreakdown_hists.root", "eff_z"));
    em.push_back(GetEff(s.dir + "/z8/breakdown/VtxEffBreakdown_hists.root", "eff_nMVTX3_zc"));
    et.push_back(GetEff(s.dir + "/zall/breakdown/VtxEffBreakdown_hists.root", "eff_trackable_pt"));
    emx.push_back(GetEff(s.dir + "/zall/breakdown/VtxEffBreakdown_hists.root", "eff_mvtxmatch_pt"));
  }
  NewPage("eff_vs_N_z8");
  EffOverlay(S, eN, ";N (charged primaries, p_{T} > 0.1, |#eta| < 1.1);#varepsilon_{vtx}", -0.5, 40.5, 0, 1.1, false,
             "Vertex efficiency vs true multiplicity, |z_{true}| < 8 cm", false);
  NewPage("eff_vs_ztrue");
  EffOverlay(S, ez, ";z_{true} (cm);#varepsilon_{vtx}", -40, 40, 0, 1.4, false, "Vertex efficiency vs true z (all events)", true);
  NewPage("eff_vs_nMVTX3_z8");
  EffOverlay(S, em, ";n_{MVTX3} (good tracks with #geq 3 MVTX clusters);#varepsilon_{vtx}", -0.5, 40.5, 0, 1.1, false,
             "Vertex efficiency vs n_{MVTX3}, |z_{true}| < 8 cm", false);
  NewPage("trackeff_vs_pt");
  EffOverlay(S, et, ";p_{T,true} (GeV/c);good-track efficiency", 0.1, 5, 0, 1.1, true,
             "Good-track efficiency, trackable charged primaries (#geq 18 TPC layers)", false);
  NewPage("mvtxmatch_vs_pt");
  EffOverlay(S, emx, ";p_{T,true} (GeV/c);fraction with #geq 2 reco MVTX clusters", 0.1, 5, 0, 1.1, true,
             "MVTX association of reconstructed primaries that crossed #geq 2 MVTX layers", false);

  // 8-10. spectra (per event, log x / log y) and vertex z residual (shape)
  auto Spectrum = [&](const char *name, const char *hname, const char *xtitle, const char *header, bool logxy) {
    auto *c = NewPage(name);
    c->SetLogy(logxy);
    c->SetLogx(logxy);
    auto *l = logxy ? Leg(0.5, 0.9 - 0.055 * S.size(), 0.94, 0.9) : Leg(0.14, 0.9 - 0.055 * S.size(), 0.6, 0.9);
    std::vector<TH1 *> hs;
    double ymax = 0;
    for (const auto &s : S)
    {
      auto *f = TFile::Open((s.dir + "/qa/EvalQA_hists.root").c_str(), "READ");
      auto *h = f ? (TH1 *) f->Get(hname) : nullptr;
      if (!h) { std::cerr << "WARNING: " << hname << " missing for " << s.label << std::endl; hs.push_back(nullptr); continue; }
      h = (TH1 *) h->Clone(Form("%s_%zu", hname, hs.size()));
      if (h->GetNbinsX() >= 100) h->Rebin(2);
      auto *ne = (TH1 *) f->Get("A_ztrue_all");  // one entry per event
      if (logxy)
        h->Scale(1.0 / (ne ? ne->GetEntries() : 1), "width");
      else if (h->Integral() > 0)
        h->Scale(1.0 / h->Integral());
      h->SetLineColor(s.col);
      h->SetLineWidth(2);
      ymax = std::max(ymax, h->GetMaximum());
      hs.push_back(h);
    }
    bool first = true;
    Out(Form("\n%s", header));
    for (size_t i = 0; i < hs.size(); ++i)
    {
      if (!hs[i]) continue;
      if (logxy)
      {
        hs[i]->GetXaxis()->SetRangeUser(0.1, 5);
        hs[i]->GetXaxis()->SetMoreLogLabels();
        hs[i]->GetXaxis()->SetNoExponent();
        hs[i]->SetMinimum(1e-2);
        hs[i]->SetMaximum(ymax * 3);
      }
      else
        hs[i]->SetMaximum(ymax * 1.7);
      hs[i]->GetXaxis()->SetTitle(xtitle);
      hs[i]->GetYaxis()->SetTitle(logxy ? "per event per GeV/c" : "fraction of vertices");
      hs[i]->Draw(first ? "hist" : "hist same");
      first = false;
      l->AddEntry(hs[i], logxy ? S[i].label.c_str() : Form("%s, RMS %.4f cm", S[i].label.c_str(), hs[i]->GetRMS()), "l");
      Out(Form("  %-40s mean %.4f  RMS %.4f", S[i].label.c_str(), hs[i]->GetMean(), hs[i]->GetRMS()));
    }
    l->Draw();
    Header(header);
  };
  Spectrum("track_pt", "C_pt", "track p_{T} (GeV/c)", "Reconstructed track p_{T} (all tracks)", true);
  Spectrum("primary_pt", "F_prim_pt_good", "p_{T,true} (GeV/c)", "Charged primaries with a good track", true);
  Spectrum("vertex_dz", "A_dz", "z_{vtx} - z_{true} (cm)", "Vertex z residual (shape)", false);

  // 11. why no vertex, |z_true| < 8 cm
  {
    auto *c = NewPage("novertex_reasons_z8");
    c->SetBottomMargin(0.16);
    const char *rl[4] = {"#splitline{no good}{TPC track}", "#splitline{good tracks, none}{with #geq 3 MVTX}", "#splitline{one track}{with #geq 3 MVTX}",
                         "#splitline{#geq 2 such tracks,}{no vertex (finder)}"};
    auto *fr = gPad->DrawFrame(-0.5, 0, 3.5, 0.6, ";;no-vertex events / all events, |z_{true}| < 8 cm");
    fr->GetXaxis()->SetNdivisions(0);
    TLatex lab;
    lab.SetTextAlign(23);
    lab.SetTextSize(0.03);
    for (int i = 0; i < 4; ++i) lab.DrawLatex(i, -0.018, rl[i]);
    auto *l = Leg(0.5, 0.9 - 0.06 * S.size(), 0.94, 0.9);
    Out("\nNo-vertex reasons inside |z_true| < 8 cm (fraction of all events in the window): no TPC track / none MVTX3 / one MVTX3 / finder");
    const double w = 0.6 / S.size();
    for (size_t j = 0; j < S.size(); ++j)
    {
      const std::string summ = S[j].dir + "/z8/breakdown/VtxEffBreakdown_summary.txt";
      double kz = 0, nz = 0;
      ReadCount(summ, "  |z_true| < 8 cm ", kz, nz);
      const auto r = ReadReasons(summ);
      if (r.size() != 4 || nz <= 0) { std::cerr << "WARNING: reasons not read from " << summ << std::endl; continue; }
      auto *g = new TGraphAsymmErrors();
      std::string row = Form("  %-40s n=%5.0f :", S[j].label.c_str(), nz);
      for (int i = 0; i < 4; ++i)
      {
        g->SetPoint(i, i - 0.3 + w * (j + 0.5), r[i] / nz);
        g->SetPointError(i, 0, 0, r[i] / nz - CPlo(nz, r[i]), CPup(nz, r[i]) - r[i] / nz);
        row += Form(" %.3f", r[i] / nz);
      }
      Style(g, S[j]);
      g->SetMarkerSize(1.5);
      g->Draw("P same");
      l->AddEntry(g, S[j].label.c_str(), "p");
      Out(row);
    }
    l->Draw();
    Header("Events without a vertex inside |z_{true}| < 8 cm, by first failing condition");
  }

  const std::string pdf = std::string(outDir) + "/CompareVtxEff.pdf";
  for (size_t i = 0; i < pages.size(); ++i)
  {
    pages[i].first->Print((pdf + (i == 0 ? "(" : (i + 1 == pages.size() ? ")" : ""))).c_str());
    pages[i].first->Print(Form("%s/png/%02zu_%s.png", outDir, i + 1, pages[i].second.c_str()));
  }
  gTable.close();
  std::cout << "wrote " << pdf << " (" << pages.size() << " pages), " << outDir << "/png/, " << outDir << "/CompareVtxEff_table.txt" << std::endl;
}
