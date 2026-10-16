// Read-only audit of MVTX association of TPC-good tracks, from SvtxEvaluator trees.
// Usage: root -l -b -q 'MvtxAudit.C+("file.list","label")'
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <tuple>
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
  int bin4(int n) { return n >= 3 ? 3 : n; }
  struct Stat
  {
    double n = 0, sum = 0, sum2 = 0;
    void add(double x) { n++; sum += x; sum2 += x * x; }
    double mean() const { return n ? sum / n : NAN; }
  };
}  // namespace

void MvtxAudit(const char *listFile, const char *label)
{
  std::vector<std::string> files;
  {
    std::ifstream in(listFile);
    for (std::string s; std::getline(in, s);)
      if (!s.empty()) files.push_back(s);
  }
  // per-class MVTX multiplicity, class: 0 unique primary, 1 additional (duplicate), 2 secondary, 3 unmatched
  const char *cls[4] = {"unique primary", "additional match", "secondary", "unmatched"};
  double cnt[4][4] = {};
  // pass/fail MVTX>=2 kinematics & purity, per class
  Stat ptS[4][2], etaS[4][2], purS[4][2], tpcpurS[4][2], mvtxpurS[4][2], ntpcS[4][2];
  double ptBins[5] = {0.1, 0.3, 0.5, 1.0, 1e9}, etaBins[4] = {0, 0.5, 0.8, 1.1};
  double ptCnt[4][2] = {}, etaCnt[3][2] = {};  // unique primaries only: [bin][pass]
  // events
  double nEv = 0, sumNtpc = 0, sumNm2 = 0;
  double sumNtpcZ[3] = {}, sumNm2Z[3] = {}, nEvZ[3] = {};  // |zvtx| <10, 10-20, >20
  // stage chain for unique primaries (TPC-good, matched to G4 primary)
  // A: crossed >=2 MVTX layers (truth g4clusters in >=2 layers)
  // B: >=2 of those layers have a reconstructed cluster (nreco>0)
  // C: track has >=2 MVTX clusters
  // D: track has >=2 MVTX clusters from the right particle (ntrumaps)
  double st[3][5] = {};  // per |z| bin: [0]=all, A,B,C,D
  // truth-level: charged primaries pT>0.1 |eta|<1.1
  double tr[3][5] = {};  // per |z| bin: all, >=2 MVTX layers crossed, TPC-good track, ..with nmaps>=2, ..ntrumaps>=2
  double wrongM[2] = {};
  const double dmax = 0.01;  // cm, truth-reco MVTX cluster spatial match
  double nTruthMvtx = 0, nTruthMvtxReco = 0;
  std::vector<double> hDmin;  // MVTX clusters on unique primaries: [0] total nmaps, [1] nwrongmaps

  for (size_t fi = 0; fi < files.size(); ++fi)
  {
    TFile f(files[fi].c_str());
    auto *ti = (TTree *) f.Get("ntp_info");
    auto *tt = (TTree *) f.Get("ntp_track");
    auto *tg = (TTree *) f.Get("ntp_gtrack");
    auto *tc = (TTree *) f.Get("ntp_g4cluster");
    auto *tv = (TTree *) f.Get("ntp_gpoint");
    if (!ti || !tt || !tg || !tc || !tv) { printf("missing tree in %s\n", files[fi].c_str()); continue; }
    // true primary vertex z per event
    std::map<int, float> gvz;
    {
      float ev, z;
      tv->SetBranchAddress("event", &ev);
      tv->SetBranchAddress("gvz", &z);
      for (Long64_t i = 0; i < tv->GetEntries(); ++i) { tv->GetEntry(i); if (!gvz.count(I(ev))) gvz[I(ev)] = z; }
    }
    auto zb = [&](int ev) { float z = gvz.count(ev) ? std::fabs(gvz[ev]) : 99; return z < 10 ? 0 : (z < 20 ? 1 : 2); };
    std::set<int> events;
    {
      float ev;
      ti->SetBranchAddress("event", &ev);
      for (Long64_t i = 0; i < ti->GetEntries(); ++i) { ti->GetEntry(i); events.insert(I(ev)); }
    }
    // reco MVTX clusters per (event, layer) from ntp_cluster (no truth link for MVTX in the evaluator)
    std::map<std::pair<int, int>, std::vector<std::array<float, 3>>> rc;
    {
      auto *tcl = (TTree *) f.Get("ntp_cluster");
      float ev, layer, x, y, z;
      tcl->SetBranchAddress("event", &ev);
      tcl->SetBranchAddress("layer", &layer);
      tcl->SetBranchAddress("x", &x);
      tcl->SetBranchAddress("y", &y);
      tcl->SetBranchAddress("z", &z);
      for (Long64_t i = 0; i < tcl->GetEntries(); ++i)
      {
        tcl->GetEntry(i);
        const int l = I(layer);
        if (l < 0 || l > 2) continue;
        rc[{I(ev), l}].push_back({x, y, z});
      }
    }
    // truth MVTX clusters per (event, gtrackID): layers crossed and layers with a reco cluster within dmax (spatial match)
    std::map<std::pair<int, int>, std::pair<int, int>> g4m;  // -> (bitmask layers crossed, bitmask layers reconstructed)
    {
      float ev, layer, gid, gx, gy, gz;
      tc->SetBranchAddress("event", &ev);
      tc->SetBranchAddress("layer", &layer);
      tc->SetBranchAddress("gtrackID", &gid);
      tc->SetBranchAddress("gx", &gx);
      tc->SetBranchAddress("gy", &gy);
      tc->SetBranchAddress("gz", &gz);
      for (Long64_t i = 0; i < tc->GetEntries(); ++i)
      {
        tc->GetEntry(i);
        int l = I(layer);
        if (l < 0 || l > 2) continue;
        auto &m = g4m[{I(ev), I(gid)}];
        m.first |= 1 << l;
        nTruthMvtx++;
        double dmin = 1e9;
        for (auto &r : rc[{I(ev), l}]) dmin = std::min(dmin, (double) std::hypot((double) r[0] - gx, (double) r[1] - gy, (double) r[2] - gz));
        hDmin.push_back(dmin);
        if (dmin < dmax) { m.second |= 1 << l; nTruthMvtxReco++; }
      }
    }
    auto nl = [](int mask) { return __builtin_popcount(mask); };
    // reco tracks
    float ev, pt, eta, ntpc, nmaps, gid, gprim, nfrom, nhits, ntrumaps, nwrongmaps, ntrutpc, gpt, geta;
    tt->SetBranchAddress("event", &ev);
    tt->SetBranchAddress("pt", &pt);
    tt->SetBranchAddress("eta", &eta);
    tt->SetBranchAddress("ntpc", &ntpc);
    tt->SetBranchAddress("nmaps", &nmaps);
    tt->SetBranchAddress("gtrackID", &gid);
    tt->SetBranchAddress("gprimary", &gprim);
    tt->SetBranchAddress("nfromtruth", &nfrom);
    tt->SetBranchAddress("nhits", &nhits);
    tt->SetBranchAddress("ntrumaps", &ntrumaps);
    tt->SetBranchAddress("nwrongmaps", &nwrongmaps);
    tt->SetBranchAddress("ntrutpc", &ntrutpc);
    // first pass: best track per (event,gid) among good tracks = highest nfromtruth
    struct T { int ev, gid, prim, ntpc, nm, ntm, nwm; float pt, eta, nfrom, nhits, ntrutpc; bool matched; };
    std::vector<T> good;
    for (Long64_t i = 0; i < tt->GetEntries(); ++i)
    {
      tt->GetEntry(i);
      if (!(fin(pt) && pt > 0.1 && fin(eta) && std::fabs(eta) < 1.1 && I(ntpc) >= 18)) continue;
      good.push_back({I(ev), fin(gid) ? I(gid) : -1, I(gprim), I(ntpc), I(nmaps), I(ntrumaps), I(nwrongmaps), pt, eta, nfrom, nhits, ntrutpc, fin(gid)});
    }
    std::map<std::pair<int, int>, size_t> best;
    for (size_t i = 0; i < good.size(); ++i)
    {
      if (!good[i].matched) continue;
      auto k = std::make_pair(good[i].ev, good[i].gid);
      if (!best.count(k) || good[i].nfrom > good[best[k]].nfrom) best[k] = i;
    }
    std::map<int, int> nT, nM2;
    for (size_t i = 0; i < good.size(); ++i)
    {
      const T &t = good[i];
      nT[t.ev]++;
      if (t.nm >= 2) nM2[t.ev]++;
      int c = !t.matched ? 3 : (best[{t.ev, t.gid}] != i ? 1 : (t.prim == 1 ? 0 : 2));
      cnt[c][bin4(t.nm)]++;
      const int p = t.nm >= 2;
      ptS[c][p].add(t.pt);
      etaS[c][p].add(std::fabs(t.eta));
      ntpcS[c][p].add(t.ntpc);
      if (t.matched && t.nhits > 0) purS[c][p].add(t.nfrom / t.nhits);
      if (t.matched && t.ntpc > 0 && fin(t.ntrutpc)) tpcpurS[c][p].add(t.ntrutpc / t.ntpc);
      if (t.matched && t.nm > 0) mvtxpurS[c][p].add((double) t.ntm / t.nm);
      if (c == 0)
      {
        for (int b = 0; b < 4; ++b) if (t.pt >= ptBins[b] && t.pt < ptBins[b + 1]) ptCnt[b][p]++;
        for (int b = 0; b < 3; ++b) if (std::fabs(t.eta) >= etaBins[b] && std::fabs(t.eta) < etaBins[b + 1]) etaCnt[b][p]++;
        wrongM[0] += t.nm;
        wrongM[1] += t.nwm;
        const int z = zb(t.ev);
        auto it = g4m.find({t.ev, t.gid});
        const int crossed = it == g4m.end() ? 0 : nl(it->second.first), recod = it == g4m.end() ? 0 : nl(it->second.second);
        st[z][0]++;
        if (crossed >= 2)
        {
          st[z][1]++;
          if (recod >= 2)
          {
            st[z][2]++;
            if (t.nm >= 2)
            {
              st[z][3]++;
              if (t.ntm >= 2) st[z][4]++;
            }
          }
        }
      }
    }
    for (int e : events)
    {
      nEv++;
      sumNtpc += nT[e];
      sumNm2 += nM2[e];
      const int z = zb(e);
      nEvZ[z]++;
      sumNtpcZ[z] += nT[e];
      sumNm2Z[z] += nM2[e];
    }
    // truth: charged primaries in acceptance
    {
      float ev2, gid2, gprim2, gflav, gpt2, geta2, trk;
      tg->SetBranchAddress("event", &ev2);
      tg->SetBranchAddress("gtrackID", &gid2);
      tg->SetBranchAddress("gprimary", &gprim2);
      tg->SetBranchAddress("gflavor", &gflav);
      tg->SetBranchAddress("gpt", &gpt2);
      tg->SetBranchAddress("geta", &geta2);
      tg->SetBranchAddress("trackID", &trk);
      // best good track per (event,gid)
      for (Long64_t i = 0; i < tg->GetEntries(); ++i)
      {
        tg->GetEntry(i);
        if (I(gprim2) != 1 || !charged(gflav) || !(gpt2 > 0.1) || !(std::fabs(geta2) < 1.1)) continue;
        const int e = I(ev2), g = I(gid2), z = zb(e);
        tr[z][0]++;
        auto it = g4m.find({e, g});
        if (it == g4m.end() || nl(it->second.first) < 2) continue;
        tr[z][1]++;
        auto b = best.find({e, g});
        if (b == best.end()) continue;
        tr[z][2]++;
        if (good[b->second].nm >= 2) tr[z][3]++;
        if (good[b->second].ntm >= 2) tr[z][4]++;
      }
    }
  }

  printf("\n################ %s : %zu files, %.0f events ################\n", label, files.size(), nEv);
  printf("mean good-track multiplicity per event: n_TPC %.3f   n_MVTX2 %.3f   (ratio %.3f)\n", sumNtpc / nEv, sumNm2 / nEv, sumNm2 / sumNtpc);
  const char *zl[3] = {"|zvtx|<10", "10-20", ">20"};
  for (int z = 0; z < 3; ++z)
    printf("   %-10s %4.0f ev: n_TPC %.2f  n_MVTX2 %.2f  ratio %.3f\n", zl[z], nEvZ[z], sumNtpcZ[z] / nEvZ[z], sumNm2Z[z] / nEvZ[z], sumNm2Z[z] / std::max(1., sumNtpcZ[z]));
  printf("\n[2] MVTX clusters on TPC-good tracks, by truth class:\n   %-18s %7s  %6s %6s %6s %6s   frac(>=2)\n", "class", "tracks", "0", "1", "2", "3+");
  double tot[4] = {};
  for (int c = 0; c < 4; ++c)
  {
    double s = 0;
    for (int b = 0; b < 4; ++b) { s += cnt[c][b]; tot[b] += cnt[c][b]; }
    printf("   %-18s %7.0f  %6.3f %6.3f %6.3f %6.3f   %.3f\n", cls[c], s, cnt[c][0] / s, cnt[c][1] / s, cnt[c][2] / s, cnt[c][3] / s, (cnt[c][2] + cnt[c][3]) / s);
  }
  double s = tot[0] + tot[1] + tot[2] + tot[3];
  printf("   %-18s %7.0f  %6.3f %6.3f %6.3f %6.3f   %.3f\n", "ALL", s, tot[0] / s, tot[1] / s, tot[2] / s, tot[3] / s, (tot[2] + tot[3]) / s);
  printf("   composition of tracks passing MVTX>=2: ");
  for (int c = 0; c < 4; ++c) printf("%s %.3f  ", cls[c], (cnt[c][2] + cnt[c][3]) / (tot[2] + tot[3]));
  printf("\n   composition of tracks failing MVTX>=2: ");
  for (int c = 0; c < 4; ++c) printf("%s %.3f  ", cls[c], (cnt[c][0] + cnt[c][1]) / (tot[0] + tot[1]));
  printf("\n\n[3] pass (MVTX>=2) vs fail: mean pT, mean|eta|, mean ntpc, purity nfromtruth/nhits, TPC purity, MVTX purity\n");
  for (int c = 0; c < 3; ++c)
    for (int p = 1; p >= 0; --p)
      printf("   %-18s %-5s N=%6.0f  pT %.3f  |eta| %.3f  ntpc %.1f  pur %.3f  tpcpur %.3f  mvtxpur %s\n", cls[c], p ? "PASS" : "fail", ptS[c][p].n, ptS[c][p].mean(), etaS[c][p].mean(), ntpcS[c][p].mean(),
             purS[c][p].mean(), tpcpurS[c][p].mean(), mvtxpurS[c][p].n ? Form("%.3f", mvtxpurS[c][p].mean()) : "  -  ");
  printf("   unique primaries, pass fraction vs pT: ");
  for (int b = 0; b < 4; ++b) printf("[%.1f,%s) %.3f (N=%.0f)  ", ptBins[b], b == 3 ? "inf" : Form("%.1f", ptBins[b + 1]), ptCnt[b][1] / (ptCnt[b][0] + ptCnt[b][1]), ptCnt[b][0] + ptCnt[b][1]);
  printf("\n   unique primaries, pass fraction vs |eta|: ");
  for (int b = 0; b < 3; ++b) printf("[%.1f,%.1f) %.3f (N=%.0f)  ", etaBins[b], etaBins[b + 1], etaCnt[b][1] / (etaCnt[b][0] + etaCnt[b][1]), etaCnt[b][0] + etaCnt[b][1]);
  printf("\n   MVTX clusters on unique primaries: %.0f, of which wrong-particle %.0f (%.3f)\n", wrongM[0], wrongM[1], wrongM[1] / std::max(1., wrongM[0]));
  std::sort(hDmin.begin(), hDmin.end());
  auto q = [&](double fr) { return hDmin.empty() ? NAN : hDmin[(size_t) (fr * (hDmin.size() - 1))] * 1e4; };
  printf("\n[4a] MVTX truth clusters (all particles): %.0f, with a reco cluster within %.0f um: %.0f (%.3f); dmin quantiles 50%% %.1f um, 90%% %.1f um, 95%% %.1f um\n", nTruthMvtx, dmax * 1e4, nTruthMvtxReco, nTruthMvtxReco / std::max(1., nTruthMvtx), q(0.5), q(0.9), q(0.95));
  printf("\n[4] stage chain, unique primary TPC-good tracks (cumulative):\n   %-10s %7s  %-22s %-22s %-22s %-22s\n", "|zvtx|", "tracks", "A crossed>=2 MVTX lay", "B >=2 lay reco clus", "C track nmaps>=2", "D ntrumaps>=2");
  for (int z = 0; z < 3; ++z)
    printf("   %-10s %7.0f  %6.0f (%.3f)         %6.0f (%.3f|A)      %6.0f (%.3f|B)      %6.0f (%.3f|C)\n", zl[z], st[z][0], st[z][1], st[z][1] / st[z][0], st[z][2], st[z][2] / std::max(1., st[z][1]), st[z][3],
           st[z][3] / std::max(1., st[z][2]), st[z][4], st[z][4] / std::max(1., st[z][3]));
  printf("\n[5] truth charged primaries (pT>0.1, |eta|<1.1), cumulative:\n   %-10s %7s  %-20s %-20s %-20s %-20s\n", "|zvtx|", "prim", ">=2 MVTX lay crossed", "TPC-good track", "track nmaps>=2", "ntrumaps>=2");
  double T[5] = {};
  for (int z = 0; z < 3; ++z)
  {
    printf("   %-10s %7.0f  %6.0f (%.3f)       %6.0f (%.3f)       %6.0f (%.3f)       %6.0f (%.3f)\n", zl[z], tr[z][0], tr[z][1], tr[z][1] / tr[z][0], tr[z][2], tr[z][2] / std::max(1., tr[z][1]), tr[z][3],
           tr[z][3] / std::max(1., tr[z][2]), tr[z][4], tr[z][4] / std::max(1., tr[z][3]));
    for (int k = 0; k < 5; ++k) T[k] += tr[z][k];
  }
  printf("   %-10s %7.0f  %6.0f (%.3f)       %6.0f (%.3f)       %6.0f (%.3f)       %6.0f (%.3f)\n", "all", T[0], T[1], T[1] / T[0], T[2], T[2] / T[1], T[3], T[3] / T[2], T[4], T[4] / T[3]);
}
