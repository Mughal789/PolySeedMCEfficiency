// =============================================================================
// PlotEvalQA.C
//
// QA of one SvtxEvaluator output file (outputeval_*.root), for comparing
// production modes (release, evaluator fixes, clock fix, ...).
//
// Usage (from analysis/macros/; analysis/run_analysis.sh calls it with a file list):
//   root -l -b -q 'PlotEvalQA.C+("SAMPLE/evaluator_files.list","SAMPLE/qa","label")'
//
// Writes into outDir:
//   EvalQA.pdf          all pages
//   EvalQA_hists.root   every histogram, fixed names (input for later comparisons)
//   EvalQA_summary.txt  key numbers
//
// Sections
//   A  events and reconstructed vertices (ntp_vertex, ntp_gpoint)
//   B  ntp_vertex placeholder rows (evaluator fix E1)
//   C  reconstructed tracks (ntp_track): kinematics, clusters, silicon seed
//   D  track vertex association and DCA (evaluator fixes E2, E4)
//   E  TPC timing diagnostics: pcaz - gvz and tan(lambda) ratio (clock fix)
//   F  charged primaries (ntp_gtrack): fate, efficiency, silicon matching
//   G  clusters (ntp_cluster): per layer, TPC z residual to truth
//   H  event multiplicities and vertex efficiency
// Good track (as VtxEffAnalysis): pt > 0.1, |eta| < 1.1, ntpc >= 18.
// N: charged G4 primaries with pt > 0.1, |eta| < 1.1.
// =============================================================================
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TParticlePDG.h>
#include <TPad.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TChain.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::string gLabel;
std::ofstream gSum;
std::vector<TObject *> gStore;  // everything written to EvalQA_hists.root

bool Fin(float x) { return std::isfinite(x); }

// Binds float branches by name; missing branches are reported once.
// Event key: the evaluator "event" counter restarts at 0 in every file, so the
// unique key is (file index in the chain, event).
struct Reader
{
  TChain *t = nullptr;
  std::map<std::string, float> v;
  Reader(const std::vector<std::string> &files, const char *name, const std::vector<std::string> &br)
  {
    t = new TChain(name);
    for (const auto &f : files) t->Add(f.c_str());
    if (t->GetNtrees() == 0 || t->GetEntries() == 0)
    {
      std::cout << "WARNING: tree " << name << " missing or empty" << std::endl;
    }
    t->SetBranchStatus("*", 0);
    for (const auto &b : br)
    {
      v[b] = std::numeric_limits<float>::quiet_NaN();
      if (!t->GetBranch(b.c_str()))
      {
        std::cout << "WARNING: " << name << "." << b << " missing" << std::endl;
        continue;
      }
      t->SetBranchStatus(b.c_str(), 1);
      t->SetBranchAddress(b.c_str(), &v[b]);
    }
  }
  float operator[](const std::string &b) { return v[b]; }
  long long key() { return 1000000LL * t->GetTreeNumber() + std::llround(v["event"]); }
  Long64_t N() const { return t ? t->GetEntries() : 0; }
};

int Charge(float f)
{
  if (!Fin(f)) return 0;
  const long long pdg = std::llround(f);
  if (std::llabs(pdg) >= 1000000000LL) return int(((std::llabs(pdg) + 5000) / 10000) % 1000);
  static std::map<long long, int> cache;
  auto it = cache.find(pdg);
  if (it != cache.end()) return it->second;
  TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(int(pdg));
  const int q = p ? int(std::lround(p->Charge() / 3.)) : 0;
  cache[pdg] = q;
  return q;
}

TH1D *H1(const char *n, const char *t, int nb, double lo, double hi)
{
  auto *h = new TH1D(n, t, nb, lo, hi);
  h->SetDirectory(nullptr);
  gStore.push_back(h);
  return h;
}
TH2D *H2(const char *n, const char *t, int nx, double xl, double xh, int ny, double yl, double yh)
{
  auto *h = new TH2D(n, t, nx, xl, xh, ny, yl, yh);
  h->SetDirectory(nullptr);
  gStore.push_back(h);
  return h;
}

void Out(const std::string &s)
{
  std::cout << s << std::endl;
  gSum << s << "\n";
}

std::vector<TCanvas *> gPages;
TCanvas *Page(const char *name, const char *title, int nx = 2, int ny = 2)
{
  auto *c = new TCanvas(name, title, 1100, 850);
  c->SetTitle(title);
  c->Divide(nx, ny, 0.004, 0.004);
  // leave the top 3.5% of the page for the page label
  for (int i = 1; i <= nx * ny; ++i)
  {
    auto *p = (TPad *) c->GetPad(i);
    p->SetPad(p->GetXlowNDC(), p->GetYlowNDC() * 0.965, p->GetXlowNDC() + p->GetWNDC(), (p->GetYlowNDC() + p->GetHNDC()) * 0.965);
  }
  gPages.push_back(c);
  return c;
}
void PageLabel(TCanvas *c)
{
  c->cd(0);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.018);
  l.SetTextColor(kGray + 2);
  l.DrawLatex(0.01, 0.985, Form("%s  |  %s", c->GetTitle(), gLabel.c_str()));
}
void Draw(TCanvas *c, int pad, TH1 *h, const char *opt = "hist", bool logy = false)
{
  c->cd(pad);
  gPad->SetLogy(logy);
  gPad->SetTopMargin(0.09);
  h->SetLineWidth(2);
  h->Draw(opt);
}
// Several histograms (same binning) overlaid with a legend.
void Overlay(TCanvas *c, int pad, std::vector<TH1 *> hs, std::vector<std::string> labs, bool norm, bool logy = false)
{
  const int cols[] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kBlack};
  c->cd(pad);
  gPad->SetLogy(logy);
  gPad->SetTopMargin(0.09);
  double mx = 0;
  std::vector<TH1 *> d;
  for (size_t i = 0; i < hs.size(); ++i)
  {
    TH1 *h = norm ? (TH1 *) hs[i]->Clone(Form("%s_n", hs[i]->GetName())) : hs[i];
    if (norm && h->Integral() > 0) h->Scale(1. / h->Integral());
    h->SetLineColor(cols[i % 6]);
    h->SetLineWidth(2);
    mx = std::max(mx, h->GetMaximum());
    d.push_back(h);
  }
  auto *lg = new TLegend(0.55, 0.72, 0.93, 0.90);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.035);
  for (size_t i = 0; i < d.size(); ++i)
  {
    d[i]->SetMaximum(logy ? mx * 20 : mx * 1.35);
    if (logy) d[i]->SetMinimum(norm ? 1e-4 : 0.5);
    d[i]->Draw(i ? "hist same" : "hist");
    lg->AddEntry(d[i], Form("%s (%.0f)", labs[i].c_str(), hs[i]->GetEntries()), "l");
  }
  lg->Draw();
}
void Text(TCanvas *c, int pad, const std::vector<std::string> &lines, double size = 0.042)
{
  c->cd(pad);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(size);
  l.SetTextFont(42);
  for (size_t i = 0; i < lines.size(); ++i) l.DrawLatex(0.04, 0.93 - i * 1.5 * size, lines[i].c_str());
}
std::string Frac(double k, double n) { return n > 0 ? Form("%.0f / %.0f = %.3f", k, n, k / n) : "-"; }
double Mode(TH1 *h) { return h->GetXaxis()->GetBinCenter(h->GetMaximumBin()); }
double FracIn(TH1 *h, double lo, double hi)
{
  const double all = h->Integral(0, h->GetNbinsX() + 1);
  return all > 0 ? h->Integral(h->FindBin(lo + 1e-9), h->FindBin(hi - 1e-9)) / all : 0;
}
}  // namespace

void PlotEvalQA(const char *evalFile, const char *outDir, const char *label = "")
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.045);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.2f");
  TH1::AddDirectory(false);
  gLabel = label;
  gSystem->mkdir(outDir, true);
  gSum.open(Form("%s/EvalQA_summary.txt", outDir));

  // evalFile: one file, a comma-separated list, or a text file *.list (one path per line)
  std::vector<std::string> files;
  {
    std::string in = evalFile;
    if (in.size() > 5 && in.substr(in.size() - 5) == ".list")
    {
      std::ifstream l(in);
      std::string x;
      while (std::getline(l, x))
        if (!x.empty()) files.push_back(x);
    }
    else
    {
      std::stringstream ss(in);
      std::string x;
      while (std::getline(ss, x, ','))
        if (!x.empty()) files.push_back(x);
    }
  }
  for (const auto &fname : files)
  {
    TFile test(fname.c_str(), "READ");
    if (test.IsZombie())
    {
      std::cerr << "FATAL: cannot open " << fname << std::endl;
      return;
    }
  }
  Out(Form("input files : %zu", files.size()));
  Out(std::string("Evaluator QA: ") + evalFile);
  Out(std::string("Label       : ") + label);

  // ===========================================================================
  // Truth vertex per event (ntp_gpoint) and event list (ntp_info)
  // ===========================================================================
  std::map<long long, std::array<float, 3>> gv;
  {
    Reader r(files, "ntp_gpoint", {"event", "gvx", "gvy", "gvz"});
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      gv[r.key()] = {r["gvx"], r["gvy"], r["gvz"]};
    }
  }
  std::set<long long> events;
  {
    Reader r(files, "ntp_info", {"event"});
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      events.insert(r.key());
    }
  }
  const double nEv = events.size();
  Out(Form("events (ntp_info): %.0f ; ntp_gpoint rows: %zu", nEv, gv.size()));

  // ===========================================================================
  // A + B: vertices and placeholder rows
  // ===========================================================================
  auto *hNvtx = H1("A_nvtx", "Reconstructed vertices per event;vertices (finite vx);events", 6, -0.5, 5.5);
  auto *hVtxRows = H1("A_vtxrows", "ntp_vertex rows per event;rows;events", 6, -0.5, 5.5);
  auto *hVx = H1("A_vx", "Vertex x;x (cm);vertices", 100, -0.5, 0.5);
  auto *hVy = H1("A_vy", "Vertex y;y (cm);vertices", 100, -0.5, 0.5);
  auto *hVz = H1("A_vz", "Vertex z;z (cm);vertices", 80, -40, 40);
  auto *hDx = H1("A_dx", "Vertex x - true x;#Deltax (cm);vertices", 100, -0.1, 0.1);
  auto *hDy = H1("A_dy", "Vertex y - true y;#Deltay (cm);vertices", 100, -0.1, 0.1);
  auto *hDz = H1("A_dz", "Vertex z - true z;#Deltaz (cm);vertices", 100, -0.2, 0.2);
  auto *hVtxNtrk = H1("A_vtxntrk", "Tracks per vertex;tracks;vertices", 41, -0.5, 40.5);
  auto *hChi2 = H1("A_chi2ndf", "Vertex #chi^{2}/ndf (NaN/0 ndf -> -1);#chi^{2}/ndf;vertices", 60, -1.5, 10);
  auto *hVtxCross = H1("A_vtxcrossing", "Vertex crossing;crossing;vertices", 21, -10.5, 10.5);
  auto *hZtAll = H1("A_ztrue_all", "True z, all events;z_{true} (cm);events", 16, -40, 40);
  auto *hZtVtx = H1("A_ztrue_vtx", "True z, vertexed events;z_{true} (cm);events", 16, -40, 40);
  // placeholder rows (vx NaN): column contents
  auto *hPhGvx = H1("B_ph_gvx", "Placeholder rows: value in column gvx;gvx column;rows", 60, -3, 3);
  auto *hPhNtrk = H1("B_ph_ntracks", "Placeholder rows: value in column ntracks;ntracks column;rows", 60, -3, 3);
  auto *hPhVtxId = H1("B_ph_vertexID", "Placeholder rows: value in column vertexID;vertexID column;rows", 60, -3, 3);
  auto *hPhGembed = H1("B_ph_gembed", "Placeholder rows: value in column gembed;gembed column;rows", 60, -3, 3);
  std::map<long long, int> nVtxEv, nRowsEv;
  double phRows = 0, phFiniteCols = 0, phBad = 0, phInVtxEvt = 0, vtxRowsN = 0;
  {
    const std::vector<std::string> cols = {"event", "vertexID", "vx", "vy", "vz", "ntracks", "chi2", "ndof", "crossing", "gvx", "gvy",
                                           "gvz", "gvt", "gembed", "gntracks", "gntracksmaps", "gnembed", "nfromtruth"};
    Reader r(files, "ntp_vertex", cols);
    // pass 1: finite vertices per event
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      const long long e = r.key();
      nRowsEv[e]++;
      if (Fin(r["vx"]) && Fin(r["vy"]) && Fin(r["vz"])) nVtxEv[e]++;
    }
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      const long long e = r.key();
      const bool fin = Fin(r["vx"]) && Fin(r["vy"]) && Fin(r["vz"]);
      if (fin)
      {
        vtxRowsN++;
        hVx->Fill(r["vx"]);
        hVy->Fill(r["vy"]);
        hVz->Fill(r["vz"]);
        hVtxNtrk->Fill(r["ntracks"]);
        hChi2->Fill(r["ndof"] > 0 ? r["chi2"] / r["ndof"] : -1.);
        hVtxCross->Fill(r["crossing"]);
        if (gv.count(e))
        {
          hDx->Fill(r["vx"] - gv[e][0]);
          hDy->Fill(r["vy"] - gv[e][1]);
          hDz->Fill(r["vz"] - gv[e][2]);
        }
      }
      else
      {
        phRows++;
        if (nVtxEv[e] > 0) phInVtxEvt++;
        hPhGvx->Fill(Fin(r["gvx"]) ? r["gvx"] : -2.9);  // NaN shown at -2.9
        hPhNtrk->Fill(Fin(r["ntracks"]) ? r["ntracks"] : -2.9);
        hPhVtxId->Fill(Fin(r["vertexID"]) ? r["vertexID"] : -2.9);
        hPhGembed->Fill(Fin(r["gembed"]) ? r["gembed"] : -2.9);
        int nf = 0;
        bool bad = false;
        for (const auto &cname : cols)
        {
          if (cname == "event") continue;
          const float x = r[cname];
          if (Fin(x)) nf++;
          if (Fin(x) && (std::fabs(x) > 1e6 || (x != 0 && std::fabs(x) < 1e-30))) bad = true;  // garbage
        }
        phFiniteCols += nf;
        phBad += bad;
      }
    }
  }
  double nVtxd = 0;
  for (long long e : events)
  {
    hNvtx->Fill(nVtxEv[e]);
    hVtxRows->Fill(nRowsEv[e]);
    nVtxd += nVtxEv[e] > 0;
    if (gv.count(e))
    {
      hZtAll->Fill(gv[e][2]);
      if (nVtxEv[e] > 0) hZtVtx->Fill(gv[e][2]);
    }
  }

  // ===========================================================================
  // C + D + E: tracks
  // ===========================================================================
  auto *hPt = H1("C_pt", "Track p_{T};p_{T} (GeV/c);tracks", 100, 0, 5);
  auto *hEta = H1("C_eta", "Track #eta;#eta;tracks", 60, -1.5, 1.5);
  auto *hPhi = H1("C_phi", "Track #phi;#phi;tracks", 64, -3.2, 3.2);
  auto *hNtpc = H1("C_ntpc", "TPC clusters (seed);ntpc;tracks", 61, -0.5, 60.5);
  auto *hNmaps = H1("C_nmaps", "MVTX clusters (seed);nmaps;tracks", 8, -0.5, 7.5);
  auto *hNintt = H1("C_nintt", "INTT clusters (seed);nintt;tracks", 8, -0.5, 7.5);
  auto *hNmms = H1("C_nmms", "TPOT clusters;nmms;tracks", 4, -0.5, 3.5);
  auto *hQual = H1("C_quality", "Track quality (#chi^{2}/ndf);quality;tracks", 100, 0, 100);
  auto *hTrkCross = H1("C_crossing", "Track crossing;crossing;tracks", 21, -10.5, 10.5);
  auto *hSil = H1("C_silseed", "Silicon seed (siqr finite);0 = no, 1 = yes;tracks", 2, -0.5, 1.5);
  auto *hComp = H1("C_truthcomp", "Good tracks: truth match;1 unmatched  2 primary  3 secondary;good tracks", 3, 0.5, 3.5);
  auto *hPur = H1("C_purity", "Good tracks: nfromtruth/nhits;purity;good tracks", 22, 0, 1.1);
  auto *hVid = H1("D_vertexID", "Track vertexID (NaN -> -3);vertexID;tracks", 14, -3.5, 10.5);
  auto *hVxTruth = H1("D_vtx_is_truth", "Track vertex position == true vertex;0 = no, 1 = yes;tracks", 2, -0.5, 1.5);
  auto *hDxyV = H1("D_dca3dxy_vtxd", "dca3dxy, tracks with vertexID >= 0;dca3dxy (cm);tracks", 100, -1, 1);
  auto *hDxyN = H1("D_dca3dxy_novtx", "dca3dxy, tracks with vertexID < 0;dca3dxy (cm);tracks", 100, -1, 1);
  auto *hDzV = H1("D_dca3dz_vtxd", "dca3dz, tracks with vertexID >= 0;dca3dz (cm);tracks", 100, -20, 20);
  auto *hDzN = H1("D_dca3dz_novtx", "dca3dz, tracks with vertexID < 0;dca3dz (cm);tracks", 100, -20, 20);
  auto *hDca2d = H1("D_dca2d_finite", "dca2d finite;0 = NaN, 1 = finite;tracks", 2, -0.5, 1.5);
  // timing diagnostics: primary-matched, ntpc > 25, pt > 0.4, |eta| < 1
  auto *hPczTn = H1("E_pcaz_tpconly_etaneg", "TPC-only primaries, #eta<0: pcaz - gvz;pcaz - gvz (cm);tracks", 100, -20, 20);
  auto *hPczTp = H1("E_pcaz_tpconly_etapos", "TPC-only primaries, #eta>0: pcaz - gvz;pcaz - gvz (cm);tracks", 100, -20, 20);
  auto *hPczS = H1("E_pcaz_mvtx", "MVTX-matched primaries: pcaz - gvz;pcaz - gvz (cm);tracks", 100, -20, 20);
  auto *hTlTn = H1("E_tanl_tpconly_etaneg", "TPC-only, #eta<0: tan#lambda reco/true;(p_{z}/p_{T})/(gp_{z}/gp_{T});tracks", 80, 0.8, 1.2);
  auto *hTlTp = H1("E_tanl_tpconly_etapos", "TPC-only, #eta>0: tan#lambda reco/true;(p_{z}/p_{T})/(gp_{z}/gp_{T});tracks", 80, 0.8, 1.2);
  auto *hTlS = H1("E_tanl_mvtx", "MVTX-matched: tan#lambda reco/true;(p_{z}/p_{T})/(gp_{z}/gp_{T});tracks", 80, 0.8, 1.2);
  std::map<long long, int> nGood, nM2, nGoodPrim;
  double nTrk = 0, nGoodAll = 0, nSil = 0, nVidNeg = 0, nVxTruth = 0, nVxTruthNovtx = 0, nNovtxTrk = 0, nDca2dFin = 0;
  {
    Reader r(files, "ntp_track", {"event", "pt", "eta", "phi", "pz", "ntpc", "nmaps", "nintt", "nmms", "quality", "crossing", "siqr",
                               "vertexID", "vx", "vy", "vz", "dca2d", "dca3dxy", "dca3dz", "pcaz", "gtrackID", "gprimary", "gpt",
                               "gpz", "gvx", "gvy", "gvz", "nfromtruth", "nhits"});
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      const long long e = r.key();
      nTrk++;
      hPt->Fill(r["pt"]);
      hEta->Fill(r["eta"]);
      hPhi->Fill(r["phi"]);
      hNtpc->Fill(r["ntpc"]);
      const int nm = Fin(r["nmaps"]) ? std::lround(r["nmaps"]) : 0;
      const int ni = Fin(r["nintt"]) ? std::lround(r["nintt"]) : 0;
      hNmaps->Fill(nm);
      hNintt->Fill(ni);
      hNmms->Fill(Fin(r["nmms"]) ? r["nmms"] : 0);
      hQual->Fill(r["quality"]);
      hTrkCross->Fill(r["crossing"]);
      const bool sil = Fin(r["siqr"]);
      hSil->Fill(sil);
      nSil += sil;
      const float vid = r["vertexID"];
      hVid->Fill(Fin(vid) ? vid : -3.);
      const bool novtx = !(Fin(vid) && vid >= 0);
      nVidNeg += novtx;
      nNovtxTrk += novtx;
      // track "vertex" position identical to the true collision vertex
      const bool isTruth = gv.count(e) && Fin(r["vz"]) && r["vx"] == gv[e][0] && r["vy"] == gv[e][1] && r["vz"] == gv[e][2];
      hVxTruth->Fill(isTruth);
      nVxTruth += isTruth;
      if (novtx) nVxTruthNovtx += isTruth;
      (novtx ? hDxyN : hDxyV)->Fill(r["dca3dxy"]);
      (novtx ? hDzN : hDzV)->Fill(r["dca3dz"]);
      hDca2d->Fill(Fin(r["dca2d"]));
      nDca2dFin += Fin(r["dca2d"]);

      const bool good = r["pt"] > 0.1 && std::fabs(r["eta"]) < 1.1 && r["ntpc"] >= 18;
      const bool matched = Fin(r["gtrackID"]) && std::lround(r["gtrackID"]) != 0;
      const bool prim = matched && Fin(r["gprimary"]) && std::lround(r["gprimary"]) == 1;
      if (good)
      {
        nGoodAll++;
        nGood[e]++;
        if (nm >= 2) nM2[e]++;
        if (prim) nGoodPrim[e]++;
        hComp->Fill(!matched ? 1 : (prim ? 2 : 3));
        if (r["nhits"] > 0) hPur->Fill(r["nfromtruth"] / r["nhits"]);
      }
      if (prim && r["ntpc"] > 25 && r["pt"] > 0.4 && std::fabs(r["eta"]) < 1.0 && Fin(r["pcaz"]) && Fin(r["gvz"]))
      {
        const double dz = r["pcaz"] - r["gvz"];
        const double tl = (r["pt"] > 0 && r["gpt"] > 0 && r["gpz"] != 0) ? (r["pz"] / r["pt"]) / (r["gpz"] / r["gpt"]) : -1;
        const bool fwdTl = std::fabs(r["gpz"] / r["gpt"]) > 0.2;  // tan(lambda) ratio only where pz is not ~0
        if (nm >= 2)
        {
          hPczS->Fill(dz);
          if (fwdTl) hTlS->Fill(tl);
        }
        else if (nm == 0 && ni == 0)
        {
          (r["eta"] < 0 ? hPczTn : hPczTp)->Fill(dz);
          if (fwdTl) (r["eta"] < 0 ? hTlTn : hTlTp)->Fill(tl);
        }
      }
    }
  }

  // ===========================================================================
  // F: charged primaries in acceptance (ntp_gtrack)
  // ===========================================================================
  const double ptb[] = {0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0};
  auto *hFate = H2("F_fate_vs_pt", "Charged primaries: fate vs p_{T};p_{T,true} (GeV/c);1 no track  2 fails good  3 good TPC-only  4 good 1 MVTX  5 good #geq2 MVTX",
                   6, 0, 6, 5, 0.5, 5.5);
  for (int b = 0; b < 6; ++b) hFate->GetXaxis()->SetBinLabel(b + 1, Form("%.1f-%.1f", ptb[b], ptb[b + 1]));
  auto *hPrimPt = H1("F_prim_pt", "Charged primaries;p_{T,true} (GeV/c);particles", 25, 0, 5);
  auto *hPrimPtGood = H1("F_prim_pt_good", "Charged primaries with a good track;p_{T,true} (GeV/c);particles", 25, 0, 5);
  auto *hPrimPtM2 = H1("F_prim_pt_mvtx2", "Charged primaries with a good #geq2-MVTX track;p_{T,true} (GeV/c);particles", 25, 0, 5);
  auto *hPrimEta = H1("F_prim_eta", "Charged primaries;#eta_{true};particles", 22, -1.1, 1.1);
  auto *hPrimEtaGood = H1("F_prim_eta_good", "Charged primaries with a good track;#eta_{true};particles", 22, -1.1, 1.1);
  auto *hPrimMvtxT = H1("F_prim_truemvtx_good", "Good-reco primaries that crossed #geq2 MVTX layers: reco MVTX clusters;nmaps (reco);particles", 8, -0.5, 7.5);
  std::map<long long, int> Ntrue;
  double nPrim = 0, nPrimTrueMvtx = 0, nPrimTrueMvtxM2 = 0;
  {
    Reader r(files, "ntp_gtrack", {"event", "gprimary", "gflavor", "gpt", "geta", "trackID", "pt", "eta", "ntpc", "nmaps", "gnlmaps"});
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      if (!(Fin(r["gprimary"]) && std::lround(r["gprimary"]) == 1)) continue;
      if (!(r["gpt"] > 0.1 && std::fabs(r["geta"]) < 1.1)) continue;
      if (Charge(r["gflavor"]) == 0) continue;
      nPrim++;
      Ntrue[r.key()]++;
      const bool reco = Fin(r["trackID"]);
      const bool good = reco && r["pt"] > 0.1 && std::fabs(r["eta"]) < 1.1 && r["ntpc"] >= 18;
      const int nm = reco && Fin(r["nmaps"]) ? std::lround(r["nmaps"]) : 0;
      const int cat = !reco ? 1 : (!good ? 2 : (nm == 0 ? 3 : (nm == 1 ? 4 : 5)));
      int b = 0;
      while (b < 5 && r["gpt"] >= ptb[b + 1]) ++b;
      hFate->Fill(b + 0.5, cat);
      hPrimPt->Fill(r["gpt"]);
      hPrimEta->Fill(r["geta"]);
      if (good)
      {
        hPrimPtGood->Fill(r["gpt"]);
        hPrimEtaGood->Fill(r["geta"]);
        if (nm >= 2) hPrimPtM2->Fill(r["gpt"]);
        if (Fin(r["gnlmaps"]) && r["gnlmaps"] >= 2)
        {
          nPrimTrueMvtx++;
          nPrimTrueMvtxM2 += nm >= 2;
          hPrimMvtxT->Fill(nm);
        }
      }
    }
  }

  // ===========================================================================
  // G: clusters
  // ===========================================================================
  auto *hClLayer = H1("G_clusters_per_layer", "Clusters per layer (per event);layer;clusters / event", 57, -0.5, 56.5);
  auto *hClDzN = H1("G_tpc_cluster_dz_south", "TPC clusters, z<0: z - gz;z - gz (cm);clusters", 100, -10, 10);
  auto *hClDzP = H1("G_tpc_cluster_dz_north", "TPC clusters, z>0: z - gz;z - gz (cm);clusters", 100, -10, 10);
  auto *pClDz = new TProfile("G_tpc_cluster_dz_vs_gz", "TPC clusters: mean z - gz vs true z;g_{z} (cm);#LTz - gz#GT (cm)", 42, -105, 105);
  pClDz->SetDirectory(nullptr);
  gStore.push_back(pClDz);
  double nClTpc = 0, nClTpcTruth = 0;
  {
    Reader r(files, "ntp_cluster", {"event", "layer", "z", "gz"});
    for (Long64_t i = 0; i < r.N(); ++i)
    {
      r.t->GetEntry(i);
      const int L = std::lround(r["layer"]);
      hClLayer->Fill(L);
      if (L >= 7 && L <= 54)
      {
        nClTpc++;
        if (Fin(r["gz"]))
        {
          nClTpcTruth++;
          const double d = r["z"] - r["gz"];
          (r["gz"] < 0 ? hClDzN : hClDzP)->Fill(d);
          if (std::fabs(d) < 20) pClDz->Fill(r["gz"], d);
        }
      }
    }
  }
  if (nEv > 0) hClLayer->Scale(1. / nEv);

  // ===========================================================================
  // H: event multiplicities and vertex efficiency
  // ===========================================================================
  auto *hnT = H1("H_n_tpc", "Good tracks per event (n_{TPC});n_{TPC};events", 101, -0.5, 100.5);
  auto *hnM = H1("H_n_mvtx2", "Good #geq2-MVTX tracks per event (n_{MVTX2});n_{MVTX2};events", 21, -0.5, 20.5);
  auto *hNt = H1("H_N_true", "Charged primaries per event (N);N;events", 101, -0.5, 100.5);
  auto *hnTv = H1("H_n_tpc_vtx", "n_{TPC}, vertexed events;n_{TPC};events", 101, -0.5, 100.5);
  auto *hnMv = H1("H_n_mvtx2_vtx", "n_{MVTX2}, vertexed events;n_{MVTX2};events", 21, -0.5, 20.5);
  auto *hNtv = H1("H_N_true_vtx", "N, vertexed events;N;events", 101, -0.5, 100.5);
  // wide bins for small samples
  const double wb[] = {-0.5, 0.5, 1.5, 2.5, 4.5, 7.5, 12.5, 20.5, 40.5};
  auto *hnTw = new TH1D("H_n_tpc_wide", ";n_{TPC};events", 8, wb);
  auto *hnTwv = new TH1D("H_n_tpc_wide_vtx", ";n_{TPC};events", 8, wb);
  auto *hnMw = new TH1D("H_n_mvtx2_wide", ";n_{MVTX2};events", 8, wb);
  auto *hnMwv = new TH1D("H_n_mvtx2_wide_vtx", ";n_{MVTX2};events", 8, wb);
  auto *hNw = new TH1D("H_N_wide", ";N;events", 8, wb);
  auto *hNwv = new TH1D("H_N_wide_vtx", ";N;events", 8, wb);
  for (TH1 *h : {(TH1 *) hnTw, (TH1 *) hnTwv, (TH1 *) hnMw, (TH1 *) hnMwv, (TH1 *) hNw, (TH1 *) hNwv})
  {
    h->SetDirectory(nullptr);
    gStore.push_back(h);
  }
  double sumT = 0, sumM = 0, sumN = 0;
  for (long long e : events)
  {
    const bool v = nVtxEv[e] > 0;
    const int a = nGood[e], m = nM2[e], N = Ntrue[e];
    sumT += a;
    sumM += m;
    sumN += N;
    hnT->Fill(a);
    hnM->Fill(m);
    hNt->Fill(N);
    hnTw->Fill(a);
    hnMw->Fill(m);
    hNw->Fill(N);
    if (v)
    {
      hnTv->Fill(a);
      hnMv->Fill(m);
      hNtv->Fill(N);
      hnTwv->Fill(a);
      hnMwv->Fill(m);
      hNwv->Fill(N);
    }
  }

  // ===========================================================================
  // Summary
  // ===========================================================================
  Out("\n[A] Events and vertices");
  Out("  events with >=1 finite vertex : " + Frac(nVtxd, nEv));
  Out(Form("  finite vertex rows            : %.0f ; mean tracks/vertex %.2f", vtxRowsN, hVtxNtrk->GetMean()));
  Out(Form("  vertex - truth  mean (RMS) um : dx %.0f (%.0f)  dy %.0f (%.0f)  dz %.0f (%.0f)", 1e4 * hDx->GetMean(),
           1e4 * hDx->GetRMS(), 1e4 * hDy->GetMean(), 1e4 * hDy->GetRMS(), 1e4 * hDz->GetMean(), 1e4 * hDz->GetRMS()));
  Out(Form("  vertex x, y mean (cm)         : %.4f, %.4f   (simulated beam spot -0.07, 0.141)", hVx->GetMean(), hVy->GetMean()));
  Out(Form("  vertices with chi2/ndf unset  : %.0f of %.0f", hChi2->GetBinContent(hChi2->FindBin(-1.)), hChi2->GetEntries()));
  Out("\n[B] ntp_vertex placeholder rows (vx NaN)  [evaluator fix E1]");
  Out(Form("  placeholder rows              : %.0f (in events that also have a vertex: %.0f)", phRows, phInVtxEvt));
  Out(Form("  finite entries per row (mean) : %.2f ; rows with garbage values: %.0f", phRows > 0 ? phFiniteCols / phRows : 0., phBad));
  Out(Form("  column gvx holds 1 (=gembed)  : %.0f rows", hPhGvx->GetBinContent(hPhGvx->FindBin(1.0))));
  Out("\n[C] Tracks");
  Out(Form("  tracks %.0f ; good %.0f ; per event %.2f", nTrk, nGoodAll, nEv > 0 ? nTrk / nEv : 0.));
  Out("  tracks with a silicon seed    : " + Frac(nSil, nTrk));
  Out(Form("  good tracks: unmatched %.0f, primary %.0f, secondary %.0f", hComp->GetBinContent(1), hComp->GetBinContent(2),
           hComp->GetBinContent(3)));
  Out("\n[D] Track-vertex association and DCA  [evaluator fixes E2, E4]");
  Out("  tracks with vertexID < 0      : " + Frac(nVidNeg, nTrk));
  Out("  track vertex == true vertex   : " + Frac(nVxTruth, nTrk) + "  (of tracks with vertexID<0: " + Frac(nVxTruthNovtx, nNovtxTrk) + ")");
  Out("  dca2d finite                  : " + Frac(nDca2dFin, nTrk));
  Out("\n[E] TPC timing (primary-matched, ntpc>25, pT>0.4, |eta|<1)  [clock fix]");
  auto tim = [&](const char *l, TH1 *dz, TH1 *tl) {
    Out(Form("  %-16s N=%5.0f  pcaz-gvz mode %6.2f mean %6.2f cm, |dz|<1 cm %.3f ; tanl ratio mean %.3f (N=%.0f)", l,
             dz->GetEntries(), Mode(dz), dz->GetMean(), FracIn(dz, -1, 1), tl->GetMean(), tl->GetEntries()));
  };
  tim("TPC-only eta<0", hPczTn, hTlTn);
  tim("TPC-only eta>0", hPczTp, hTlTp);
  tim("MVTX-matched", hPczS, hTlS);
  Out("\n[F] Charged primaries pT>0.1, |eta|<1.1");
  Out(Form("  primaries %.0f ; per event %.2f", nPrim, nEv > 0 ? nPrim / nEv : 0.));
  {
    double c[6] = {};
    for (int b = 1; b <= 6; ++b)
      for (int k = 1; k <= 5; ++k) c[k] += hFate->GetBinContent(b, k);
    Out(Form("  no track %.3f | fails good %.3f | good TPC-only %.3f | good 1 MVTX %.3f | good >=2 MVTX %.3f", c[1] / nPrim,
             c[2] / nPrim, c[3] / nPrim, c[4] / nPrim, c[5] / nPrim));
  }
  Out("  good primaries crossing >=2 MVTX layers, reco >=2 MVTX: " + Frac(nPrimTrueMvtxM2, nPrimTrueMvtx));
  Out("\n[G] Clusters");
  Out("  TPC clusters with finite truth gz: " + Frac(nClTpcTruth, nClTpc));
  Out(Form("  TPC cluster z-gz mean: south %.3f cm, north %.3f cm", hClDzN->GetMean(), hClDzP->GetMean()));
  Out("\n[H] Event multiplicity");
  Out(Form("  mean n_TPC %.3f ; mean n_MVTX2 %.3f ; mean N %.3f", sumT / nEv, sumM / nEv, sumN / nEv));
  Out("  vertex fraction               : " + Frac(nVtxd, nEv));

  // ===========================================================================
  // Pages
  // ===========================================================================
  auto eff = [](TH1 *pass, TH1 *tot, const char *name, const char *title) {
    auto *e = new TEfficiency(*pass, *tot);
    e->SetName(name);
    e->SetTitle(title);
    e->SetStatisticOption(TEfficiency::kFCP);
    e->SetConfidenceLevel(0.683);
    e->SetMarkerStyle(20);
    gStore.push_back(e);
    return e;
  };
  {
    auto *c = Page("pA1", "A. Vertices", 3, 2);
    Draw(c, 1, hNvtx, "hist text0");
    Draw(c, 2, hVtxNtrk);
    Draw(c, 3, hChi2);
    Draw(c, 4, hVx);
    Draw(c, 5, hVy);
    Draw(c, 6, hVz);
    PageLabel(c);
  }
  {
    auto *c = Page("pA2", "A. Vertex resolution and z acceptance", 3, 2);
    Draw(c, 1, hDx);
    Draw(c, 2, hDy);
    Draw(c, 3, hDz);
    Draw(c, 4, hVtxCross);
    Overlay(c, 5, {hZtAll, hZtVtx}, {"all events", "vertexed"}, false);
    c->cd(6);
    eff(hZtVtx, hZtAll, "eff_ztrue", "Vertex efficiency vs true z;z_{true} (cm);#varepsilon_{vtx}")->Draw("AP");
    PageLabel(c);
  }
  {
    auto *c = Page("pB", "B. ntp_vertex placeholder rows (E1)", 3, 2);
    Draw(c, 1, hVtxRows, "hist text0");
    Draw(c, 2, hPhVtxId);
    Draw(c, 3, hPhNtrk);
    Draw(c, 4, hPhGvx);
    Draw(c, 5, hPhGembed);
    Text(c, 6,
         {"Placeholder rows (vx = NaN):", Form("rows: %.0f", phRows), Form("in events with a vertex: %.0f", phInVtxEvt),
          Form("finite entries/row: %.2f", phRows > 0 ? phFiniteCols / phRows : 0.), Form("rows with garbage: %.0f", phBad),
          "NaN entries are drawn at -2.9", "Correct: vertexID, ntracks NaN;", "gvx = true x, gembed = 1"},
         0.05);
    PageLabel(c);
  }
  {
    auto *c = Page("pC1", "C. Tracks: kinematics", 3, 2);
    Draw(c, 1, hPt, "hist", true);
    Draw(c, 2, hEta);
    Draw(c, 3, hPhi);
    Draw(c, 4, hQual, "hist", true);
    Draw(c, 5, hTrkCross, "hist text0");
    Draw(c, 6, hSil, "hist text0");
    PageLabel(c);
  }
  {
    auto *c = Page("pC2", "C. Tracks: clusters and truth", 3, 2);
    Draw(c, 1, hNtpc);
    Draw(c, 2, hNmaps, "hist text0");
    Draw(c, 3, hNintt, "hist text0");
    Draw(c, 4, hNmms, "hist text0");
    Draw(c, 5, hComp, "hist text0");
    Draw(c, 6, hPur);
    PageLabel(c);
  }
  {
    auto *c = Page("pD", "D. Track-vertex association and DCA (E2, E4)", 3, 2);
    Draw(c, 1, hVid, "hist text0");
    Draw(c, 2, hVxTruth, "hist text0");
    Draw(c, 3, hDca2d, "hist text0");
    Overlay(c, 4, {hDxyV, hDxyN}, {"vertexID #geq 0", "vertexID < 0"}, true, true);
    Overlay(c, 5, {hDzV, hDzN}, {"vertexID #geq 0", "vertexID < 0"}, true, true);
    Text(c, 6,
         {"vertexID < 0: " + Frac(nVidNeg, nTrk), "vertex == truth: " + Frac(nVxTruth, nTrk), "dca2d finite: " + Frac(nDca2dFin, nTrk),
          "If 'vertex == truth' is non-zero, DCAs", "of those tracks are measured from the", "TRUE vertex (bug E2)."},
         0.05);
    PageLabel(c);
  }
  {
    auto *c = Page("pE", "E. TPC timing: z offset and tan#lambda (clock fix)", 3, 2);
    Draw(c, 1, hPczTn);
    Draw(c, 2, hPczTp);
    Draw(c, 3, hPczS);
    Draw(c, 4, hTlTn);
    Draw(c, 5, hTlTp);
    Draw(c, 6, hTlS);
    PageLabel(c);
  }
  {
    auto *c = Page("pF", "F. Charged primaries: fate and efficiency", 3, 2);
    c->cd(1);
    gPad->SetRightMargin(0.14);
    auto *fr = (TH2D *) hFate->Clone("F_fate_frac");
    for (int b = 1; b <= 6; ++b)
    {
      double s = 0;
      for (int k = 1; k <= 5; ++k) s += fr->GetBinContent(b, k);
      for (int k = 1; k <= 5; ++k) fr->SetBinContent(b, k, s > 0 ? fr->GetBinContent(b, k) / s : 0);
    }
    fr->SetTitle("Fate fraction per p_{T} bin;p_{T,true};1 none 2 fail 3 TPC-only 4 1MVTX 5 #geq2MVTX");
    gStore.push_back(fr);
    fr->Draw("colz text");
    c->cd(2);
    eff(hPrimPtGood, hPrimPt, "eff_good_vs_pt", "Good-track efficiency;p_{T,true} (GeV/c);efficiency")->Draw("AP");
    c->cd(3);
    eff(hPrimPtM2, hPrimPt, "eff_mvtx2_vs_pt", "Good #geq2-MVTX efficiency;p_{T,true} (GeV/c);efficiency")->Draw("AP");
    c->cd(4);
    eff(hPrimEtaGood, hPrimEta, "eff_good_vs_eta", "Good-track efficiency;#eta_{true};efficiency")->Draw("AP");
    Draw(c, 5, hPrimMvtxT, "hist text0");
    Text(c, 6, {Form("primaries: %.0f", nPrim), "crossed #geq2 MVTX layers, good:", "  reco #geq2 MVTX: " + Frac(nPrimTrueMvtxM2, nPrimTrueMvtx)},
         0.05);
    PageLabel(c);
  }
  {
    auto *c = Page("pG", "G. Clusters", 2, 2);
    Draw(c, 1, hClLayer);
    Draw(c, 2, hClDzN);
    Draw(c, 3, hClDzP);
    c->cd(4);
    pClDz->SetMarkerStyle(20);
    pClDz->Draw();
    PageLabel(c);
  }
  {
    auto *c = Page("pH", "H. Multiplicity and vertex efficiency", 3, 2);
    Overlay(c, 1, {hnT, hnTv}, {"all", "vertexed"}, false);
    Overlay(c, 2, {hnM, hnMv}, {"all", "vertexed"}, false);
    Overlay(c, 3, {hNt, hNtv}, {"all", "vertexed"}, false);
    c->cd(4);
    eff(hnTwv, hnTw, "eff_n_tpc", "#varepsilon_{vtx}(n_{TPC});n_{TPC};#varepsilon_{vtx}")->Draw("AP");
    c->cd(5);
    eff(hnMwv, hnMw, "eff_n_mvtx2", "#varepsilon_{vtx}(n_{MVTX2});n_{MVTX2};#varepsilon_{vtx}")->Draw("AP");
    c->cd(6);
    eff(hNwv, hNw, "eff_N", "#varepsilon_{vtx}(N);N;#varepsilon_{vtx}")->Draw("AP");
    PageLabel(c);
  }

  const std::string pdf = std::string(outDir) + "/EvalQA.pdf";
  gSystem->mkdir(Form("%s/png", outDir), true);
  for (size_t i = 0; i < gPages.size(); ++i)
  {
    gPages[i]->Print((pdf + (i == 0 ? "(" : (i + 1 == gPages.size() ? ")" : ""))).c_str(), Form("Title:%s", gPages[i]->GetName()));
    gPages[i]->Print(Form("%s/png/%02zu_%s.png", outDir, i + 1, gPages[i]->GetName()));
  }
  TFile out(Form("%s/EvalQA_hists.root", outDir), "RECREATE");
  for (TObject *o : gStore) o->Write();
  TNamed("label", label).Write();
  TNamed("input", evalFile).Write();
  out.Close();
  Out(Form("\nwrote %s (%zu pages), %s/EvalQA_hists.root, %s/EvalQA_summary.txt", pdf.c_str(), gPages.size(), outDir, outDir));
  gSum.close();
}
