// =============================================================================
// VtxEffBreakdown.C - breakdown of the vertex efficiency of one sample: efficiency vs
// multiplicity and z, why collisions have no vertex, why tracks are missing,
// reconstructed vs true multiplicity
//
// Usage (from analysis/macros/; analysis/run_analysis.sh calls it per z window):
//   root -l -b -q 'VtxEffBreakdown.C+("SAMPLE/z8/VtxEff_TPC.root","SAMPLE/evaluator_files.list","SAMPLE/z8/breakdown","label"[,zCentral[,mvtxMin]])'
//   zCentral (default 10 cm): central |z| window. The analysis uses |z_vtx| < 8 cm:
//   run with 8. The summary then also gives the efficiency
//   with the cut applied to the RECONSTRUCTED vertex (zReco) and the z migration.
//   mvtxMin (default 2): MVTX clusters a good track needs to count as MVTX-matched
//   (n_MVTX2 or n_MVTX3 from the evt tree). 3 = PHSimpleVertexFinder setNmvtxRequired(3).
//
// Inputs
//   VtxEff_TPC.root      VtxEffAnalysis output (per-event tree "evt"; n = good
//                        tracks pt>0.1, |eta|<1.1, ntpc>=18; N = charged G4
//                        primaries pt>0.1, |eta|<1.1; hasVtx = >=1 tracking vertex)
//   evaluator list       the step's evaluator files (ntp_gtrack, ntp_g4cluster)
// Output (outDir): VtxEffBreakdown.pdf, png/, VtxEffBreakdown_summary.txt, VtxEffBreakdown_hists.root
//
// Sections of the summary and pages of the PDF (histogram names keep the prefixes Q1/Q2/Q3/Q5):
// 1.  vertex efficiency and why events have no vertex (Q1_*)
//     no-vertex reasons (exclusive, in this order):
//       R1 |z_true| > 20 cm (outside the MVTX acceptance)
//       R2 no good track at all (n_TPC = 0)
//       R3 good tracks, none with >= 2 MVTX clusters (n_MVTX2 = 0)
//       R4 exactly one MVTX-matched good track (n_MVTX2 = 1)
//       R5 n_MVTX2 >= 2 and still no vertex (vertex finder)
// 2.  events with 8-12 good tracks: vertexed vs not (Q2_*)
// 3.  charged primaries: why no track / why no MVTX match (Q3_*)
//     (truth TPC/MVTX layers per particle from ntp_g4cluster, gprimary == 1)
// 4.  epsilon(n) vs epsilon(N): n versus N, composition of n (Q5_*)
// =============================================================================
#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TParticlePDG.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
std::ofstream gSum;
std::string gLabel;
std::vector<TObject *> gStore;
std::vector<TCanvas *> gPages;

void Out(const std::string &s)
{
  std::cout << s << std::endl;
  gSum << s << "\n";
}
std::string Frac(double k, double n) { return n > 0 ? Form("%5.0f / %5.0f = %.3f", k, n, k / n) : "-"; }
std::string CP(double k, double n)
{
  if (n <= 0) return "-";
  const double e = k / n, lo = TEfficiency::ClopperPearson(n, k, 0.683, false), up = TEfficiency::ClopperPearson(n, k, 0.683, true);
  return Form("%.3f +%.3f -%.3f  (%.0f/%.0f)", e, up - e, e - lo, k, n);
}
int Charge(float f)
{
  if (!std::isfinite(f)) return 0;
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
template <class T>
T *Keep(T *h)
{
  h->SetDirectory(nullptr);
  gStore.push_back(h);
  return h;
}
TEfficiency *Eff(TH1 *pass, TH1 *tot, const char *name, const char *title)
{
  auto *e = new TEfficiency(*pass, *tot);
  e->SetName(name);
  e->SetTitle(title);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(0.683);
  e->SetMarkerStyle(20);
  gStore.push_back(e);
  return e;
}
TCanvas *Page(const char *name, const char *title, int nx = 3, int ny = 2)
{
  auto *c = new TCanvas(name, title, 1100, 850);
  c->SetTitle(title);
  c->Divide(nx, ny, 0.004, 0.004);
  for (int i = 1; i <= nx * ny; ++i)
  {
    auto *p = (TPad *) c->GetPad(i);
    p->SetPad(p->GetXlowNDC(), p->GetYlowNDC() * 0.965, p->GetXlowNDC() + p->GetWNDC(), (p->GetYlowNDC() + p->GetHNDC()) * 0.965);
    p->SetTopMargin(0.10);
  }
  gPages.push_back(c);
  return c;
}
void Label(TCanvas *c)
{
  c->cd(0);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.018);
  l.SetTextColor(kGray + 2);
  l.DrawLatex(0.01, 0.985, Form("%s  |  %s", c->GetTitle(), gLabel.c_str()));
}
void Text(TCanvas *c, int pad, const std::vector<std::string> &lines, double size = 0.045)
{
  c->cd(pad);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(size);
  for (size_t i = 0; i < lines.size(); ++i) l.DrawLatex(0.04, 0.92 - i * 1.5 * size, lines[i].c_str());
}
void DrawEff(TCanvas *c, int pad, TEfficiency *e)
{
  c->cd(pad);
  e->Draw("AP");
  gPad->Update();
  if (e->GetPaintedGraph())
  {
    e->GetPaintedGraph()->SetMinimum(0);
    e->GetPaintedGraph()->SetMaximum(1.05);
  }
}
void DrawH(TCanvas *c, int pad, TH1 *h, const char *opt = "hist text0", int col = kBlue + 1)
{
  c->cd(pad);
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetMinimum(0);
  h->SetMaximum(h->GetMaximum() * 1.3 + 1);
  h->Draw(opt);
}
}  // namespace

void VtxEffBreakdown(const char *vtxeffFile, const char *evalList, const char *outDir, const char *label = "", double zCentral = 10, int mvtxMin = 2)
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  gStyle->SetPaintTextFormat("4.0f");
  gStyle->SetPalette(kBird);
  TH1::AddDirectory(false);
  gLabel = label;
  gSystem->mkdir(Form("%s/png", outDir), true);
  gSum.open(Form("%s/VtxEffBreakdown_summary.txt", outDir));
  Out(std::string("VtxEffBreakdown: ") + label);
  Out(std::string("  VtxEff file : ") + vtxeffFile);
  Out(std::string("  evaluator   : ") + evalList);

  // ===========================================================================
  // Event level (VtxEffAnalysis tree)
  // ===========================================================================
  TFile fv(vtxeffFile, "READ");
  auto *t = dynamic_cast<TTree *>(fv.Get("evt"));
  if (!t)
  {
    std::cerr << "FATAL: no evt tree in " << vtxeffFile << std::endl;
    return;
  }
  int n, N, hasVtx, nMVTX2, nGoodPrim, nGoodSec, nGoodDup, nUniquePrim, vtxNtracks, nAll;
  float zTrue, zReco;
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("N", &N);
  t->SetBranchAddress("hasVtx", &hasVtx);
  t->SetBranchAddress(mvtxMin == 3 ? "nMVTX3" : "nMVTX2", &nMVTX2);  // variable keeps its name; label via ML()
  t->SetBranchAddress("nGoodPrim", &nGoodPrim);
  t->SetBranchAddress("nGoodSec", &nGoodSec);
  t->SetBranchAddress("nGoodDup", &nGoodDup);
  t->SetBranchAddress("nUniquePrim", &nUniquePrim);
  t->SetBranchAddress("vtxNtracks", &vtxNtracks);
  t->SetBranchAddress("nAll", &nAll);
  t->SetBranchAddress("zTrue", &zTrue);
  t->SetBranchAddress("zReco", &zReco);
  if (mvtxMin != 2 && mvtxMin != 3)
  {
    std::cerr << "FATAL: mvtxMin must be 2 or 3" << std::endl;
    return;
  }
  // labels: "MVTX2" -> "MVTX<mvtxMin>", "#geq2 MVTX" -> "#geq<mvtxMin> MVTX"
  auto ML = [mvtxMin](std::string x)
  {
    if (mvtxMin == 2) return x;
    const std::string to = "MVTX" + std::to_string(mvtxMin);
    for (size_t q = 0; (q = x.find("MVTX2", q)) != std::string::npos; q += to.size()) x.replace(q, 5, to);
    const std::string geq = "#geq" + std::to_string(mvtxMin) + " MVTX";
    for (size_t q = 0; (q = x.find("#geq2 MVTX", q)) != std::string::npos; q += geq.size()) x.replace(q, 10, geq);
    return x;
  };

  const double wb[] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 8.5, 12.5, 16.5, 24.5, 40.5};
  const int nwb = 11;
  auto mk = [&](const char *nm, const char *t1) { return Keep(new TH1D(nm, t1, nwb, wb)); };
  TH1D *hnA = mk("Q1_nTPC_all", ";n_{TPC};events"), *hnV = mk("Q1_nTPC_vtx", ";n_{TPC};events");
  TH1D *hmA = mk("Q1_nMVTX2_all", ML(";n_{MVTX2};events").c_str()), *hmV = mk("Q1_nMVTX2_vtx", ML(";n_{MVTX2};events").c_str());
  TH1D *hNA = mk("Q1_N_all", ";N;events"), *hNV = mk("Q1_N_vtx", ";N;events");
  TH1D *hzA = Keep(new TH1D("Q1_z_all", ";z_{true} (cm);events", 16, -40, 40));
  TH1D *hzV = Keep(new TH1D("Q1_z_vtx", ";z_{true} (cm);events", 16, -40, 40));
  TH1D *hmAc = mk("Q1_nMVTX2_all_zc", ML(";n_{MVTX2};events").c_str()), *hmVc = mk("Q1_nMVTX2_vtx_zc", ML(";n_{MVTX2};events").c_str());
  TH1D *hReason = Keep(new TH1D("Q1_novtx_reason", "Events without a vertex: reason;;events", 5, 0.5, 5.5));
  const std::string rl[5] = {"|z|>20 cm", "n_{TPC}=0", ML("n_{MVTX2}=0"), ML("n_{MVTX2}=1"), ML("n_{MVTX2}#geq2")};
  for (int i = 0; i < 5; ++i) hReason->GetXaxis()->SetBinLabel(i + 1, rl[i].c_str());
  TH1D *hReason812 = (TH1D *) hReason->Clone("Q2_novtx_reason_n8to12");
  hReason812->SetTitle("8 #leq n_{TPC} #leq 12, no vertex: reason;;events");
  Keep(hReason812);
  // Q5
  auto *pnN = Keep(new TProfile("Q5_nTPC_vs_N", "Mean n vs true N;N;#LTn#GT", 21, -0.5, 20.5));
  auto *pmN = Keep(new TProfile("Q5_nMVTX2_vs_N", ML("Mean n_{MVTX2} vs N;N;#LTn_{MVTX2}#GT").c_str(), 21, -0.5, 20.5));
  auto *ppN = Keep(new TProfile("Q5_nPrimUnique_vs_N", "Mean unique primary-matched good tracks vs N;N;#LTn_{prim,unique}#GT", 21, -0.5, 20.5));
  auto *hComp = Keep(new TH1D("Q5_composition", "Composition of good tracks;;tracks", 3, 0.5, 3.5));
  hComp->GetXaxis()->SetBinLabel(1, "primary, unique");
  hComp->GetXaxis()->SetBinLabel(2, "duplicate");
  hComp->GetXaxis()->SetBinLabel(3, "secondary");
  auto *hPA = mk("Q5_nPrimU_all", ";n_{prim,unique};events"), *hPV = mk("Q5_nPrimU_vtx", ";n_{prim,unique};events");

  double nEv = 0, nVtx = 0, nVtxZc = 0, nZc = 0, nZ20 = 0, nVtxZ20 = 0;
  double nVtxZcRecoZc = 0, nRecoZc = 0, nRecoZcTrueOut = 0;  // cut on the reconstructed vertex
  double reasonZc[6] = {};                                       // no-vertex reasons inside |z_true| < zCentral
  double g[4][2] = {}, gzc[4][2] = {};  // groups n_MVTX2 = 0, 1, >=2, >=3
  // Q2: 8..12
  double q2[2] = {}, q2m2[2] = {}, q2m0[2] = {}, q2z[2] = {}, q2prim[2] = {}, q2sec[2] = {}, q2dup[2] = {}, q2N[2] = {}, q2nm[2] = {};
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    nEv++;
    nVtx += hasVtx;
    const bool zc = std::fabs(zTrue) < zCentral, z20 = std::fabs(zTrue) < 20;
    const bool recoZc = hasVtx && std::isfinite(zReco) && std::fabs(zReco) < zCentral;
    nZc += zc;
    nVtxZc += zc && hasVtx;
    nVtxZcRecoZc += zc && recoZc;
    nRecoZc += recoZc;
    nRecoZcTrueOut += recoZc && !zc;
    nZ20 += z20;
    nVtxZ20 += z20 && hasVtx;
    hnA->Fill(n);
    hmA->Fill(nMVTX2);
    hNA->Fill(N);
    hzA->Fill(zTrue);
    hPA->Fill(nUniquePrim);
    if (zc) hmAc->Fill(nMVTX2);
    if (hasVtx)
    {
      hnV->Fill(n);
      hmV->Fill(nMVTX2);
      hNV->Fill(N);
      hzV->Fill(zTrue);
      hPV->Fill(nUniquePrim);
      if (zc) hmVc->Fill(nMVTX2);
    }
    const int gi[4] = {nMVTX2 == 0, nMVTX2 == 1, nMVTX2 >= 2, nMVTX2 >= 3};
    for (int k = 0; k < 4; ++k)
      if (gi[k])
      {
        g[k][0]++;
        g[k][1] += hasVtx;
        if (zc)
        {
          gzc[k][0]++;
          gzc[k][1] += hasVtx;
        }
      }
    int reason = 0;
    if (!hasVtx)
    {
      reason = !z20 ? 1 : (n == 0 ? 2 : (nMVTX2 == 0 ? 3 : (nMVTX2 == 1 ? 4 : 5)));
      hReason->Fill(reason);
      if (zc) reasonZc[reason]++;
    }
    pnN->Fill(N, n);
    pmN->Fill(N, nMVTX2);
    ppN->Fill(N, nUniquePrim);
    hComp->Fill(1, nUniquePrim);
    hComp->Fill(2, nGoodDup);
    hComp->Fill(3, nGoodSec);
    if (n >= 8 && n <= 12)
    {
      const int v = hasVtx ? 1 : 0;
      q2[v]++;
      q2m2[v] += (nMVTX2 >= 2);
      q2m0[v] += (nMVTX2 == 0);
      q2nm[v] += nMVTX2;
      q2z[v] += !z20;
      q2prim[v] += nGoodPrim;
      q2sec[v] += nGoodSec;
      q2dup[v] += nGoodDup;
      q2N[v] += N;
      if (!hasVtx) hReason812->Fill(reason);
    }
  }

  Out("\n=== 1. Vertex efficiency (every collision in the denominator) ===");
  Out("  all collisions                 : " + CP(nVtx, nEv));
  Out("  |z_true| < 20 cm               : " + CP(nVtxZ20, nZ20));
  Out(Form("  %-31s: ", Form("|z_true| < %g cm", zCentral)) + CP(nVtxZc, nZc));
  Out(Form("  |z_true| < %g cm AND vertex |z_reco| < %g cm : ", zCentral, zCentral) + CP(nVtxZcRecoZc, nZc));
  Out(Form("  events selected by |z_reco| < %g cm: %.0f, of which |z_true| >= %g cm: %.0f",
           zCentral, nRecoZc, zCentral, nRecoZcTrueOut));
  const std::string gn[4] = {ML("n_MVTX2 = 0 "), ML("n_MVTX2 = 1 "), ML("n_MVTX2 >= 2"), ML("n_MVTX2 >= 3")};
  Out("  group           fraction of events          eps (all z)                     " + std::string(Form("eps (|z|<%g)", zCentral)));
  for (int k = 0; k < 4; ++k)
    Out(Form("  %s   %-26s  %-30s  %s", gn[k].c_str(), Frac(g[k][0], nEv).c_str(), CP(g[k][1], g[k][0]).c_str(), CP(gzc[k][1], gzc[k][0]).c_str()));
  Out("  events WITHOUT a vertex, by reason (exclusive, in this order):");
  const double nNo = nEv - nVtx;
  const char *rlt[5] = {"R1 |z_true| > 20 cm (MVTX acceptance)", "R2 no good track (n_TPC = 0)", "R3 good tracks, none MVTX-matched", "R4 one MVTX-matched good track", "R5 >= 2 MVTX-matched, no vertex"};
  for (int r = 1; r <= 5; ++r) Out(Form("    %-40s %s", rlt[r - 1], Frac(hReason->GetBinContent(r), nNo).c_str()));
  if (zCentral < 20)
  {
    Out(Form("  events WITHOUT a vertex inside |z_true| < %g cm, by reason:", zCentral));
    for (int r = 2; r <= 5; ++r) Out(Form("    %-40s %s", rlt[r - 1], Frac(reasonZc[r], nZc - nVtxZc).c_str()));
  }

  Out("\n=== 2. Events with 8 <= n_TPC <= 12 ===");
  Out(Form("  %-40s %12s %12s", "", "NO vertex", "vertexed"));
  auto q2row = [&](const char *l, double *a, bool mean) {
    Out(Form("  %-40s %12.3f %12.3f", l, mean ? a[0] / q2[0] : a[0] / q2[0], mean ? a[1] / q2[1] : a[1] / q2[1]));
  };
  Out(Form("  %-40s %12.0f %12.0f", "events", q2[0], q2[1]));
  q2row(ML("mean MVTX-matched good tracks (n_MVTX2)").c_str(), q2nm, true);
  q2row(ML("fraction with n_MVTX2 = 0").c_str(), q2m0, false);
  q2row(ML("fraction with n_MVTX2 >= 2").c_str(), q2m2, false);
  q2row("fraction |z_true| > 20 cm", q2z, false);
  q2row("mean primary-matched good tracks", q2prim, true);
  q2row("mean secondary-matched good tracks", q2sec, true);
  q2row("mean duplicate good tracks", q2dup, true);
  q2row("mean true N", q2N, true);
  Out("  no-vertex events by reason:");
  for (int r = 1; r <= 5; ++r) Out(Form("    %-40s %s", rlt[r - 1], Frac(hReason812->GetBinContent(r), q2[0]).c_str()));

  // ===========================================================================
  // Track level: charged primaries (ntp_gtrack) and their truth layers (ntp_g4cluster)
  // ===========================================================================
  std::vector<std::string> files;
  {
    std::ifstream l(evalList);
    std::string x;
    while (std::getline(l, x))
      if (!x.empty()) files.push_back(x);
  }
  // truth layers per primary particle: key (file, event, gtrackID) -> bitmask of layers
  std::map<std::tuple<int, int, int>, std::pair<unsigned, unsigned long long>> layers;  // (mvtx/intt bits 0..6, tpc bits)
  {
    TChain c("ntp_g4cluster");
    for (auto &f : files) c.Add(f.c_str());
    float ev, layer, gtrackID, gprimary;
    c.SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "gtrackID", "gprimary"}) c.SetBranchStatus(b, 1);
    c.SetBranchAddress("event", &ev);
    c.SetBranchAddress("layer", &layer);
    c.SetBranchAddress("gtrackID", &gtrackID);
    c.SetBranchAddress("gprimary", &gprimary);
    for (Long64_t i = 0; i < c.GetEntries(); ++i)
    {
      c.GetEntry(i);
      if (!(std::isfinite(gprimary) && std::lround(gprimary) == 1)) continue;
      const int L = std::lround(layer);
      auto &e = layers[{c.GetTreeNumber(), int(std::lround(ev)), int(std::lround(gtrackID))}];
      if (L >= 0 && L <= 6) e.first |= (1u << L);
      else if (L >= 7 && L <= 54) e.second |= (1ULL << (L - 7));
    }
  }
  auto nbits = [](unsigned long long x) { return __builtin_popcountll(x); };

  const double ptb[] = {0.1, 0.2, 0.3, 0.5, 1.0, 5.0};
  auto *hMissTpc = Keep(new TH1D("Q3_missed_truthTPClayers", "Primaries WITHOUT a track: TPC layers with truth hits;TPC layers hit;particles", 49, -0.5, 48.5));
  auto *hRecoTpc = Keep(new TH1D("Q3_reco_truthTPClayers", "Primaries WITH a good track: TPC layers with truth hits;TPC layers hit;particles", 49, -0.5, 48.5));
  auto *hMissCat = Keep(new TH1D("Q3_missed_category", "Primaries without a track;;particles", 4, 0.5, 4.5));
  const char *mcl[4] = {"0 TPC layers", "1-17 layers", "18-29 layers", "#geq30 layers"};
  for (int i = 0; i < 4; ++i) hMissCat->GetXaxis()->SetBinLabel(i + 1, mcl[i]);
  auto *hPt = Keep(new TH1D("Q3_prim_pt", ";p_{T,true} (GeV/c);particles", 5, ptb));
  auto *hPtReco = Keep(new TH1D("Q3_prim_pt_reco", ";p_{T,true} (GeV/c);particles", 5, ptb));
  auto *hPtTrk = Keep(new TH1D("Q3_prim_pt_trackable", ";p_{T,true} (GeV/c);particles", 5, ptb));
  auto *hPtTrkReco = Keep(new TH1D("Q3_prim_pt_trackable_reco", ";p_{T,true} (GeV/c);particles", 5, ptb));
  // MVTX matching for good-track primaries that truly crossed >= 2 MVTX layers
  auto *hMxPt = Keep(new TH1D("Q3_mvtxcand_pt", ";p_{T,true} (GeV/c);particles", 5, ptb));
  auto *hMxPtM = Keep(new TH1D("Q3_mvtxcand_pt_matched", ";p_{T,true} (GeV/c);particles", 5, ptb));
  auto *hMxEta = Keep(new TH1D("Q3_mvtxcand_eta", ";#eta_{true};particles", 11, -1.1, 1.1));
  auto *hMxEtaM = Keep(new TH1D("Q3_mvtxcand_eta_matched", ";#eta_{true};particles", 11, -1.1, 1.1));
  auto *hMxZ = Keep(new TH1D("Q3_mvtxcand_z", ";z_{vertex,true} (cm);particles", 8, -20, 20));
  auto *hMxZM = Keep(new TH1D("Q3_mvtxcand_z_matched", ";z_{vertex,true} (cm);particles", 8, -20, 20));
  double nPrim = 0, nNoTrk = 0, nNoTrkTrackable = 0, nTrackable = 0, nNoTrkZ = 0, nMissCat[5] = {};
  {
    TChain c("ntp_gtrack");
    for (auto &f : files) c.Add(f.c_str());
    float ev, gtrackID, gprimary, gflavor, gpt, geta, gvz, trackID, pt, eta, ntpc, nmaps;
    c.SetBranchStatus("*", 0);
    std::map<std::string, float *> br = {{"event", &ev}, {"gtrackID", &gtrackID}, {"gprimary", &gprimary}, {"gflavor", &gflavor},
                                         {"gpt", &gpt}, {"geta", &geta}, {"gvz", &gvz}, {"trackID", &trackID}, {"pt", &pt},
                                         {"eta", &eta}, {"ntpc", &ntpc}, {"nmaps", &nmaps}};
    for (auto &b : br)
    {
      c.SetBranchStatus(b.first.c_str(), 1);
      c.SetBranchAddress(b.first.c_str(), b.second);
    }
    for (Long64_t i = 0; i < c.GetEntries(); ++i)
    {
      c.GetEntry(i);
      if (!(std::isfinite(gprimary) && std::lround(gprimary) == 1)) continue;
      if (!(gpt > 0.1 && std::fabs(geta) < 1.1) || Charge(gflavor) == 0) continue;
      nPrim++;
      const auto it = layers.find({c.GetTreeNumber(), int(std::lround(ev)), int(std::lround(gtrackID))});
      const int ltpc = it != layers.end() ? nbits(it->second.second) : 0;
      const int lmvtx = it != layers.end() ? nbits(it->second.first & 0x7u) : 0;
      const bool reco = std::isfinite(trackID);
      const bool good = reco && pt > 0.1 && std::fabs(eta) < 1.1 && ntpc >= 18;
      const int nm = reco && std::isfinite(nmaps) ? std::lround(nmaps) : 0;
      const bool trackable = ltpc >= 18;
      hPt->Fill(gpt);
      if (good) hPtReco->Fill(gpt);
      if (trackable)
      {
        nTrackable++;
        hPtTrk->Fill(gpt);
        if (good) hPtTrkReco->Fill(gpt);
      }
      if (!reco)
      {
        nNoTrk++;
        nNoTrkZ += std::fabs(gvz) > 20;
        hMissTpc->Fill(ltpc);
        const int cat = ltpc == 0 ? 1 : (ltpc < 18 ? 2 : (ltpc < 30 ? 3 : 4));
        hMissCat->Fill(cat);
        nMissCat[cat]++;
        if (trackable) nNoTrkTrackable++;
      }
      if (good)
      {
        hRecoTpc->Fill(ltpc);
        if (lmvtx >= 2 && std::fabs(gvz) < 20)
        {
          hMxPt->Fill(gpt);
          hMxEta->Fill(geta);
          hMxZ->Fill(gvz);
          if (nm >= 2)
          {
            hMxPtM->Fill(gpt);
            hMxEtaM->Fill(geta);
            hMxZM->Fill(gvz);
          }
        }
      }
    }
  }
  Out("\n=== 3. Charged primaries, pT > 0.1 GeV/c, |eta| < 1.1 ===");
  Out(Form("  primaries %.0f ; without any reconstructed track: %s", nPrim, Frac(nNoTrk, nPrim).c_str()));
  Out("  primaries without a track, by TPC layers with truth hits:");
  for (int k = 1; k <= 4; ++k) Out(Form("    %-16s %s", mcl[k - 1], Frac(nMissCat[k], nNoTrk).c_str()));
  Out("    (of these, produced at |z_true| > 20 cm: " + Frac(nNoTrkZ, nNoTrk) + ")");
  Out("  'trackable' primaries (>= 18 TPC layers with truth hits): " + Frac(nTrackable, nPrim));
  Out("  good-track efficiency, all primaries      : " + CP(hPtReco->GetEntries(), nPrim));
  Out("  good-track efficiency, trackable primaries: " + CP(hPtTrkReco->GetEntries(), nTrackable));
  Out("  MVTX matching (good track, truly crossed >= 2 MVTX layers, |z|<20): " + CP(hMxPtM->GetEntries(), hMxPt->GetEntries()));

  Out("\n=== 4. epsilon(n) vs epsilon(N) ===");
  const double tot = hComp->Integral();
  Out(Form("  good tracks: unique primary %.3f, duplicate %.3f, secondary %.3f  (fractions of %.0f)", hComp->GetBinContent(1) / tot,
           hComp->GetBinContent(2) / tot, hComp->GetBinContent(3) / tot, tot));
  Out(ML("   N   <n_TPC>   <n_MVTX2>  <n_prim,unique>").c_str());
  for (int b = 1; b <= 13; ++b)
    if (pnN->GetBinEntries(b) > 0)
      Out(Form("  %2d   %6.2f    %6.2f     %6.2f", b - 1, pnN->GetBinContent(b), pmN->GetBinContent(b), ppN->GetBinContent(b)));

  // ===========================================================================
  // Pages
  // ===========================================================================
  {
    auto *c = Page("eff_vs_multiplicity_z", "Vertex efficiency vs multiplicity and z");
    DrawEff(c, 1, Eff(hnV, hnA, "eff_nTPC", "#varepsilon_{vtx} vs n_{TPC};n_{TPC} (good tracks);#varepsilon_{vtx}"));
    DrawEff(c, 2, Eff(hmV, hmA, ML("eff_nMVTX2").c_str(), ML("#varepsilon_{vtx} vs n_{MVTX2};n_{MVTX2} (good tracks #geq2 MVTX);#varepsilon_{vtx}").c_str()));
    DrawEff(c, 3, Eff(hNV, hNA, "eff_N", "#varepsilon_{vtx} vs N;N (charged primaries);#varepsilon_{vtx}"));
    DrawEff(c, 4, Eff(hzV, hzA, "eff_z", "#varepsilon_{vtx} vs z_{true};z_{true} (cm);#varepsilon_{vtx}"));
    DrawEff(c, 5, Eff(hmVc, hmAc, ML("eff_nMVTX2_zc").c_str(), Form(ML("#varepsilon_{vtx} vs n_{MVTX2}, |z_{true}|<%g cm;n_{MVTX2};#varepsilon_{vtx}").c_str(), zCentral)));
    Text(c, 6,
         {"#varepsilon_{vtx} (every collision counted):", "all: " + CP(nVtx, nEv), "|z|<20: " + CP(nVtxZ20, nZ20), Form("|z|<%g: ", zCentral) + CP(nVtxZc, nZc),
          ML("n_{MVTX2}=0: ").c_str() + CP(g[0][1], g[0][0]), ML("n_{MVTX2}=1: ").c_str() + CP(g[1][1], g[1][0]), ML("n_{MVTX2}#geq2: ").c_str() + CP(g[2][1], g[2][0])},
         0.042);
    Label(c);
  }
  {
    auto *c = Page("no_vertex_reasons", "Why collisions have no vertex", 2, 2);
    DrawH(c, 1, hReason);
    DrawH(c, 2, hReason812);
    Text(c, 3, {"R1 |z_{true}| > 20 cm: outside MVTX acceptance", "R2 no good track", "R3 good tracks, none with #geq2 MVTX clusters",
                "R4 exactly one MVTX-matched good track", "R5 #geq2 MVTX-matched good tracks, no vertex", "(exclusive, assigned in this order)"},
         0.05);
    Text(c, 4,
         {"8 #leq n_{TPC} #leq 12:", Form("events: %.0f no vertex, %.0f vertexed", q2[0], q2[1]),
          Form(ML("mean n_{MVTX2}: %.2f vs %.2f").c_str(), q2nm[0] / std::max(1., q2[0]), q2nm[1] / std::max(1., q2[1])),
          Form(ML("n_{MVTX2}=0: %.2f vs %.2f").c_str(), q2m0[0] / std::max(1., q2[0]), q2m0[1] / std::max(1., q2[1])),
          Form("|z|>20 cm: %.2f vs %.2f", q2z[0] / std::max(1., q2[0]), q2z[1] / std::max(1., q2[1]))},
         0.05);
    Label(c);
  }
  {
    auto *c = Page("missing_tracks", "Why tracks are missing / not MVTX-matched");
    DrawH(c, 1, hMissCat);
    c->cd(2);
    hRecoTpc->SetLineColor(kBlue + 1);
    hMissTpc->SetLineColor(kRed + 1);
    hRecoTpc->SetLineWidth(2);
    hMissTpc->SetLineWidth(2);
    hRecoTpc->SetMaximum(std::max(hRecoTpc->GetMaximum(), hMissTpc->GetMaximum()) * 1.3);
    hRecoTpc->SetTitle("TPC layers with truth hits;TPC layers hit;particles");
    hRecoTpc->Draw("hist");
    hMissTpc->Draw("hist same");
    auto *lg = new TLegend(0.15, 0.72, 0.6, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hRecoTpc, "with a good track", "l");
    lg->AddEntry(hMissTpc, "without a track", "l");
    lg->Draw();
    DrawEff(c, 3, Eff(hPtTrkReco, hPtTrk, "eff_trackable_pt", "Good-track eff., trackable primaries (#geq18 TPC layers);p_{T,true};efficiency"));
    DrawEff(c, 4, Eff(hMxPtM, hMxPt, "eff_mvtxmatch_pt", "MVTX match vs p_{T} (crossed #geq2 MVTX layers);p_{T,true};matched"));
    DrawEff(c, 5, Eff(hMxEtaM, hMxEta, "eff_mvtxmatch_eta", "MVTX match vs #eta;#eta_{true};matched"));
    DrawEff(c, 6, Eff(hMxZM, hMxZ, "eff_mvtxmatch_z", "MVTX match vs vertex z;z_{true} (cm);matched"));
    Label(c);
  }
  {
    auto *c = Page("reco_vs_true_multiplicity", "#varepsilon(n) vs #varepsilon(N): what n contains");
    c->cd(1);
    pnN->SetMarkerStyle(20);
    pnN->SetMarkerColor(kBlack);
    pnN->SetMinimum(0);
    pnN->Draw();
    pmN->SetMarkerStyle(21);
    pmN->SetMarkerColor(kRed + 1);
    pmN->SetLineColor(kRed + 1);
    pmN->Draw("same");
    ppN->SetMarkerStyle(22);
    ppN->SetMarkerColor(kBlue + 1);
    ppN->SetLineColor(kBlue + 1);
    ppN->Draw("same");
    auto *lg = new TLegend(0.15, 0.70, 0.62, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(pnN, "n_{TPC}", "p");
    lg->AddEntry(ppN, "unique primary-matched", "p");
    lg->AddEntry(pmN, ML("n_{MVTX2}").c_str(), "p");
    lg->Draw();
    gStyle->SetPaintTextFormat("4.0f");
    DrawH(c, 2, hComp);
    DrawEff(c, 3, Eff(hPV, hPA, "eff_nPrimUnique", "#varepsilon_{vtx} vs unique primary-matched good tracks;n_{prim,unique};#varepsilon_{vtx}"));
    c->cd(4);
    for (auto *o : gStore)
      if (std::string(o->GetName()) == "eff_N") ((TEfficiency *) o)->Draw("AP");
    for (auto *o : gStore)
      if (std::string(o->GetName()) == "eff_nTPC")
      {
        auto *e = (TEfficiency *) o;
        e->SetMarkerColor(kRed + 1);
        e->SetLineColor(kRed + 1);
        e->SetMarkerStyle(24);
        e->Draw("P same");
      }
    Text(c, 5, {"#varepsilon(N) (black) vs #varepsilon(n_{TPC}) (red, open)", "x axes differ: n #approx 2N"}, 0.05);
    Label(c);
  }
  const std::string pdf = std::string(outDir) + "/VtxEffBreakdown.pdf";
  for (size_t i = 0; i < gPages.size(); ++i)
  {
    gPages[i]->Print((pdf + (i == 0 ? "(" : (i + 1 == gPages.size() ? ")" : ""))).c_str());
    gPages[i]->Print(Form("%s/png/%02zu_%s.png", outDir, i + 1, gPages[i]->GetName()));
  }
  TFile out(Form("%s/VtxEffBreakdown_hists.root", outDir), "RECREATE");
  for (TObject *o : gStore) o->Write();
  out.Close();
  Out(Form("\nwrote %s (%zu pages), VtxEffBreakdown_summary.txt, VtxEffBreakdown_hists.root", pdf.c_str(), gPages.size()));
  gSum.close();
}
