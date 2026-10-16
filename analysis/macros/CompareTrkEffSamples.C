// =============================================================================
// CompareTrkEffSamples.C - consistency of two independent samples analysed with
// the same TrkEffAnalysis configuration (e.g. K1 1000 events vs the 10k sample)
//
// For every listed output file (same name in both directories) and species:
// integrated efficiency over the in-range bins (or the species window for the
// final maps), binomial errors, difference and pull. Also the event-level
// quantities (selected events / all events) and the pion charge difference.
// The samples are independent (disjoint seeds), so the pulls should scatter
// around 0 with unit width.
//
// Usage: root -l -b -q 'CompareTrkEffSamples.C+("DIR_A","DIR_B","label A","label B","out.txt")'
// =============================================================================
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace cmpsmp
{
struct E
{
  double n = 0, d = 0;
  double e() const { return d > 0 ? n / d : NAN; }
  double s() const { return d > 0 ? std::sqrt(std::max(e() * (1 - e()), 1.0 / d / d) / d) : NAN; }
};
E Integ(TFile &f, const char *sp, double ylo, double yhi, double plo, double phi)
{
  E r;
  auto *d = f.Get<TH2D>(Form("hDen_%s_ypt", sp));
  auto *n = f.Get<TH2D>(Form("hNum_%s_ypt", sp));
  if (!d || !n) return r;
  for (int i = 1; i <= d->GetNbinsX(); ++i)
    for (int j = 1; j <= d->GetNbinsY(); ++j)
    {
      const double y = d->GetXaxis()->GetBinCenter(i), p = d->GetYaxis()->GetBinCenter(j);
      if (y > ylo && y < yhi && p > plo && p < phi)
      {
        r.d += d->GetBinContent(i, j);
        r.n += n->GetBinContent(i, j);
      }
    }
  return r;
}
}  // namespace cmpsmp

void CompareTrkEffSamples(const char *dirA, const char *dirB, const char *labA, const char *labB, const char *outPath)
{
  using namespace cmpsmp;
  FILE *out = fopen(outPath, "w");
  auto P = [&](const std::string &s) {
    printf("%s\n", s.c_str());
    fprintf(out, "%s\n", s.c_str());
  };
  const char *sp[6] = {"pip", "pim", "kp", "km", "p", "pbar"};
  struct File
  {
    const char *name, *label;
    bool window;  // true: species analysis windows (final maps), false: all in-range bins
  };
  const std::vector<File> files = {{"TrkEff_baseline.root", "Step 2 baseline (true-z, TPC 18, MVTX 0)", false},
                                   {"TrkEff_recoz.root", "Step 2 reco-z", false},
                                   {"variants/TrkEff_mvtx2_ztrue_incl.root", "4A true-z MVTX 2", false},
                                   {"variants/TrkEff_mvtx3_ztrue_incl.root", "4A true-z MVTX 3", false},
                                   {"variants/TrkEff_mvtx2_zreco_incl.root", "4A reco-z MVTX 2", false},
                                   {"variants/TrkEff_mvtx3_zreco_incl.root", "4A reco-z MVTX 3", false},
                                   {"final_mvtx3/TrkEff_final_mvtx3.root", "FINAL (species windows)", true},
                                   {"final_mvtx3/TrkEff_final_mvtx2ref.root", "FINAL with MVTX 2 (windows)", true}};
  P(Form("A = %s (%s)\nB = %s (%s)", labA, dirA, labB, dirB));
  int npull = 0;
  double sumPull = 0, sumPull2 = 0;
  for (const auto &f : files)
  {
    std::unique_ptr<TFile> a(TFile::Open((std::string(dirA) + "/" + f.name).c_str())),
        b(TFile::Open((std::string(dirB) + "/" + f.name).c_str()));
    if (!a || !b || a->IsZombie() || b->IsZombie())
    {
      P(Form("-- %s: missing in one sample, skipped", f.name));
      continue;
    }
    auto *sa = a->Get<TH1D>("hTrkEffStats"), *sb = b->Get<TH1D>("hTrkEffStats");
    P("");
    P(Form("== %s   [selected events: A %.0f/%.0f, B %.0f/%.0f]", f.label, sa->GetBinContent(3), sa->GetBinContent(2),
           sb->GetBinContent(3), sb->GetBinContent(2)));
    P("   species   eps A (D)              eps B (D)              B - A        pull");
    E pa[6], pb[6];
    for (int s = 0; s < 6; ++s)
    {
      double ylo = -99, yhi = 99, plo = -1, phi = 99;
      if (f.window)
      {
        const double yw[6] = {1.0, 1.0, 0.7, 0.7, 0.6, 0.6}, pl[6] = {0.1, 0.1, 0.1, 0.1, 0.4, 0.4}, ph[6] = {99, 99, 0.6, 0.6, 1.0, 1.0};
        ylo = -yw[s];
        yhi = yw[s];
        plo = pl[s];
        phi = ph[s];
      }
      pa[s] = Integ(*a, sp[s], ylo, yhi, plo, phi);
      pb[s] = Integ(*b, sp[s], ylo, yhi, plo, phi);
      const double diff = pb[s].e() - pa[s].e(), sig = std::hypot(pa[s].s(), pb[s].s()), pull = diff / sig;
      P(Form("   %-6s  %.4f+-%.4f (%6.0f)   %.4f+-%.4f (%6.0f)   %+.4f   %+5.2f", sp[s], pa[s].e(), pa[s].s(), pa[s].d, pb[s].e(),
             pb[s].s(), pb[s].d, diff, pull));
      if (std::isfinite(pull))
      {
        ++npull;
        sumPull += pull;
        sumPull2 += pull * pull;
      }
    }
    for (const E *x : {pa, pb})
    {
      const double d = x[1].e() - x[0].e(), s = std::hypot(x[0].s(), x[1].s());
      P(Form("   %s: eps(pi-) - eps(pi+) = %+.4f +- %.4f  (%.1f sigma)", x == pa ? "A" : "B", d, s, d / s));
    }
  }
  const double mean = sumPull / npull, rms = std::sqrt(sumPull2 / npull - mean * mean);
  P("");
  P(Form("pulls B vs A: n = %d, mean %+.2f, rms %.2f (independent samples: expect ~0 and ~1; entries are correlated "
         "across selections of the same sample)",
         npull, mean, rms));
  fclose(out);
}
