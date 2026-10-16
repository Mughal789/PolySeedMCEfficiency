#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>
// median truth-reco MVTX cluster distance vs phi (20 deg bins), layers 0-2 combined
void ResidPhi(const char *list)
{
  std::map<int, std::vector<double>> d;
  std::ifstream in(list);
  for (std::string fn; std::getline(in, fn);)
  {
    TFile f(fn.c_str());
    auto *tc = (TTree *) f.Get("ntp_g4cluster");
    auto *tr = (TTree *) f.Get("ntp_cluster");
    std::map<std::pair<int, int>, std::vector<std::array<float, 3>>> rc;
    float ev, layer, x, y, z;
    tr->SetBranchAddress("event", &ev); tr->SetBranchAddress("layer", &layer);
    tr->SetBranchAddress("x", &x); tr->SetBranchAddress("y", &y); tr->SetBranchAddress("z", &z);
    for (Long64_t i = 0; i < tr->GetEntries(); ++i) { tr->GetEntry(i); int l = std::lround(layer); if (l >= 0 && l <= 2) rc[{(int) std::lround(ev), l}].push_back({x, y, z}); }
    tc->SetBranchAddress("event", &ev); tc->SetBranchAddress("layer", &layer);
    tc->SetBranchAddress("gx", &x); tc->SetBranchAddress("gy", &y); tc->SetBranchAddress("gz", &z);
    for (Long64_t i = 0; i < tc->GetEntries(); ++i)
    {
      tc->GetEntry(i);
      int l = std::lround(layer);
      if (l < 0 || l > 2) continue;
      double best = 1e9;
      for (auto &r : rc[{(int) std::lround(ev), l}]) best = std::min(best, (double) std::hypot((double) r[0] - x, (double) r[1] - y, (double) r[2] - z));
      if (best > 0.3) continue;
      int b = (int) std::floor((std::atan2(y, x) * 180 / M_PI + 180) / 20);
      d[b].push_back(best * 1e4);
    }
  }
  for (auto &[b, v] : d)
  {
    std::sort(v.begin(), v.end());
    printf("phi %5d..%4d : median %4.0f um  (N=%zu)\n", -180 + 20 * b, -160 + 20 * b, v[v.size() / 2], v.size());
  }
}
