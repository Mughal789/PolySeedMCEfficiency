// Read-only: pT, eta, phi of reconstructed vs lost charged primaries (SvtxEvaluator trees).
// Truth: charged G4 primaries, pT > 0.1, |eta| < 1.1 (ntp_gtrack).
// Reconstructed = has a TPC-good track (pt>0.1, |eta|<1.1, ntpc>=18; best by nfromtruth).
// MVTX-associated = that track has nmaps >= 3 (locked definition = vertex finder).
// Categories by truth TPC layers crossed (gnltpc): 0 | 1-17 | >=18 ("trackable").
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TParticlePDG.h>
#include <TStyle.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace
{
  bool fin(float x) { return std::isfinite(x); }
  int I(float x) { return fin(x) ? (int) std::lround(x) : 0; }
  bool charged(float f)
  {
    if (!fin(f)) return false;
    auto *p = TDatabasePDG::Instance()->GetParticle((int) f);
    return p && std::fabs(p->Charge()) > 0.1;
  }
}  // namespace

void LostTracks(const char *listFile, const char *outPdf, const char *label)
{
  gStyle->SetOptStat(0);
  std::vector<std::string> files;
  {
    std::ifstream in(listFile);
    for (std::string s; std::getline(in, s);)
      if (!s.empty()) files.push_back(s);
  }
  const double ptb[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 3.0};
  const int npt = 10;
  // denominators/numerators: [sel] sel 0 = all primaries, 1 = trackable (gnltpc>=18), 2 = |gvz|<8 trackable
  struct H { TH1D *pt, *eta, *phi, *phis, *z; };
  auto mk = [&](const char *n) {
    H h;
    h.pt = new TH1D(Form("pt_%s", n), ";p_{T} (GeV/c)", npt, ptb);
    h.eta = new TH1D(Form("eta_%s", n), ";#eta", 11, -1.1, 1.1);
    h.phi = new TH1D(Form("phi_%s", n), ";#phi (rad)", 36, -M_PI, M_PI);
    h.phis = new TH1D(Form("phis_%s", n), ";#phi mod 30#circ (deg)", 15, 0, 30);
    h.z = new TH1D(Form("z_%s", n), ";z_{vtx,true} (cm)", 16, -40, 40);
    return h;
  };
  auto fill = [](H &h, float pt, float eta, float phi, float z) {
    h.pt->Fill(std::min(pt, 2.99f));
    h.eta->Fill(eta);
    h.phi->Fill(phi);
    double d = std::fmod(phi * 180 / M_PI + 360., 30.);
    h.phis->Fill(d);
    h.z->Fill(z);
  };
  H hAll = mk("all"), hRec = mk("rec"), hTr = mk("trk"), hTrRec = mk("trkrec"), hTrM3 = mk("trkm3"), hRecM = mk("recm"), hRecM3 = mk("recm3");
  // hRecM: reconstructed AND crossed >=2 MVTX layers, |zvtx|<10 (MVTX reachable) ; hRecM3: of those, nmaps>=3
  double cat[3] = {}, catRec[3] = {};
  // lost trackable primaries: list summary
  double nLostTr = 0, lostLowPt = 0, lostEdgeEta = 0;
  for (auto &fn : files)
  {
    TFile f(fn.c_str());
    auto *tt = (TTree *) f.Get("ntp_track");
    auto *tg = (TTree *) f.Get("ntp_gtrack");
    auto *tc = (TTree *) f.Get("ntp_g4cluster");
    std::map<std::pair<int, int>, int> mvtxLayers;  // truth MVTX layers crossed
    {
      float ev, layer, gid;
      tc->SetBranchAddress("event", &ev);
      tc->SetBranchAddress("layer", &layer);
      tc->SetBranchAddress("gtrackID", &gid);
      for (Long64_t i = 0; i < tc->GetEntries(); ++i)
      {
        tc->GetEntry(i);
        const int l = I(layer);
        if (l >= 0 && l <= 2) mvtxLayers[{I(ev), I(gid)}] |= 1 << l;
      }
    }
    // best TPC-good track per (event, gtrackID)
    std::map<std::pair<int, int>, std::pair<float, int>> best;  // -> (nfromtruth, nmaps)
    {
      float ev, pt, eta, ntpc, nmaps, gid, nfrom;
      tt->SetBranchAddress("event", &ev);
      tt->SetBranchAddress("pt", &pt);
      tt->SetBranchAddress("eta", &eta);
      tt->SetBranchAddress("ntpc", &ntpc);
      tt->SetBranchAddress("nmaps", &nmaps);
      tt->SetBranchAddress("gtrackID", &gid);
      tt->SetBranchAddress("nfromtruth", &nfrom);
      for (Long64_t i = 0; i < tt->GetEntries(); ++i)
      {
        tt->GetEntry(i);
        if (!(fin(pt) && pt > 0.1 && fin(eta) && std::fabs(eta) < 1.1 && I(ntpc) >= 18 && fin(gid))) continue;
        auto k = std::make_pair(I(ev), I(gid));
        if (!best.count(k) || nfrom > best[k].first) best[k] = {nfrom, I(nmaps)};
      }
    }
    float ev, gid, gprim, gflav, gpt, geta, gphi, gvz, gnltpc;
    tg->SetBranchAddress("event", &ev);
    tg->SetBranchAddress("gtrackID", &gid);
    tg->SetBranchAddress("gprimary", &gprim);
    tg->SetBranchAddress("gflavor", &gflav);
    tg->SetBranchAddress("gpt", &gpt);
    tg->SetBranchAddress("geta", &geta);
    tg->SetBranchAddress("gphi", &gphi);
    tg->SetBranchAddress("gvz", &gvz);
    tg->SetBranchAddress("gnltpc", &gnltpc);
    for (Long64_t i = 0; i < tg->GetEntries(); ++i)
    {
      tg->GetEntry(i);
      if (I(gprim) != 1 || !charged(gflav) || !(gpt > 0.1) || !(std::fabs(geta) < 1.1)) continue;
      auto k = std::make_pair(I(ev), I(gid));
      auto b = best.find(k);
      const bool rec = b != best.end();
      const int nl = I(gnltpc);
      const int c = nl == 0 ? 0 : (nl < 18 ? 1 : 2);
      cat[c]++;
      catRec[c] += rec;
      fill(hAll, gpt, geta, gphi, gvz);
      if (rec) fill(hRec, gpt, geta, gphi, gvz);
      if (c == 2)
      {
        fill(hTr, gpt, geta, gphi, gvz);
        if (rec) fill(hTrRec, gpt, geta, gphi, gvz);
        else
        {
          nLostTr++;
          lostLowPt += gpt < 0.2;
          lostEdgeEta += std::fabs(geta) > 0.9;
        }
      }
      auto m = mvtxLayers.find(k);
      const bool reach = m != mvtxLayers.end() && __builtin_popcount(m->second) >= 2 && std::fabs(gvz) < 10;
      if (rec && reach)
      {
        fill(hRecM, gpt, geta, gphi, gvz);
        if (b->second.second >= 3) fill(hRecM3, gpt, geta, gphi, gvz);
      }
    }
  }
  printf("\n##### %s #####\n", label);
  const char *cn[3] = {"0 TPC layers (never reach TPC)", "1-17 TPC layers", ">=18 TPC layers (trackable)"};
  double tot = cat[0] + cat[1] + cat[2];
  printf("charged primaries pT>0.1 |eta|<1.1: %.0f\n", tot);
  for (int c = 0; c < 3; ++c) printf("   %-32s %5.0f (%.3f of all)  reconstructed %.3f  lost %4.0f\n", cn[c], cat[c], cat[c] / tot, catRec[c] / cat[c], cat[c] - catRec[c]);
  printf("   lost trackable: %.0f ; with pT<0.2: %.0f (%.2f) ; with |eta|>0.9: %.0f (%.2f)\n", nLostTr, lostLowPt, lostLowPt / nLostTr, lostEdgeEta, lostEdgeEta / nLostTr);
  auto effTable = [&](const char *title, TH1D *num, TH1D *den) {
    printf("  %s\n    ", title);
    for (int i = 1; i <= den->GetNbinsX(); ++i)
      printf("[%g,%g) %.2f(%.0f)  ", den->GetXaxis()->GetBinLowEdge(i), den->GetXaxis()->GetBinUpEdge(i), den->GetBinContent(i) ? num->GetBinContent(i) / den->GetBinContent(i) : NAN, den->GetBinContent(i));
    printf("\n");
  };
  printf("\nTPC-good reconstruction efficiency of TRACKABLE primaries (>=18 TPC layers):\n");
  effTable("vs pT", hTrRec.pt, hTr.pt);
  effTable("vs eta", hTrRec.eta, hTr.eta);
  effTable("vs phi mod 30 deg", hTrRec.phis, hTr.phis);
  effTable("vs z_vtx", hTrRec.z, hTr.z);
  printf("\nMVTX>=3 association of reconstructed primaries that crossed >=2 MVTX layers, |z_vtx|<10:\n");
  effTable("vs pT", hRecM3.pt, hRecM.pt);
  effTable("vs eta", hRecM3.eta, hRecM.eta);
  effTable("vs phi", hRecM3.phi, hRecM.phi);
  // PDF
  TCanvas c("c", "", 1500, 1000);
  c.Print(Form("%s[", outPdf));
  // build per-variable pairs
  auto pick = [](H &h, int iv) { TH1D *a[5] = {h.pt, h.eta, h.phi, h.phis, h.z}; return a[iv]; };
  for (int pg = 0; pg < 2; ++pg)
  {
    c.Clear();
    c.Divide(3, 2);
    for (int iv = 0; iv < 5; ++iv)
    {
      c.cd(iv + 1);
      if (iv == 0) gPad->SetLogx();
      auto *lg = new TLegend(0.15, 0.13, 0.88, 0.30);
      std::vector<std::pair<TH1D *, TH1D *>> v;
      std::vector<const char *> lab;
      if (pg == 0)
      {
        v = {{pick(hTrRec, iv), pick(hTr, iv)}, {pick(hRec, iv), pick(hAll, iv)}};
        lab = {"TPC-good, trackable primaries (#geq18 TPC layers)", "TPC-good, all charged primaries"};
      }
      else
      {
        v = {{pick(hRecM3, iv), pick(hRecM, iv)}};
        lab = {"MVTX#geq3 of reco primaries crossing #geq2 MVTX layers, |z|<10"};
      }
      for (size_t k = 0; k < v.size(); ++k)
      {
        auto *e = new TEfficiency(*v[k].first, *v[k].second);
        e->SetTitle(Form("%s;%s;efficiency", label, v[k].second->GetXaxis()->GetTitle()));
        e->SetMarkerStyle(20 + k);
        e->SetMarkerColor(k == 0 ? kBlue + 1 : kRed + 1);
        e->SetLineColor(k == 0 ? kBlue + 1 : kRed + 1);
        e->Draw(k == 0 ? "AP" : "P same");
        gPad->Update();
        if (k == 0 && e->GetPaintedGraph()) { e->GetPaintedGraph()->SetMinimum(0); e->GetPaintedGraph()->SetMaximum(1.05); gPad->Update(); }
        lg->AddEntry(e, lab[k], "p");
      }
      lg->SetTextSize(0.035);
      lg->Draw();
    }
    c.Print(outPdf);
  }
  c.Print(Form("%s]", outPdf));
  printf("wrote %s\n", outPdf);
}
