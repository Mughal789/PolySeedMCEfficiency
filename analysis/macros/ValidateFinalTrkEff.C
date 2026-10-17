// =============================================================================
// ValidateFinalTrkEff.C - checks of the FINAL tracking-efficiency selection
// (analysis/TRKEFF.md, "Final selection") on the outputs written by run_trkeff_final.sh
//
//  A. one-cut steps on the fine grid, reco vertex, no truth-eta cut:
//       step0  TPC >= 18, MVTX >= 2
//       step1  + TPC >= 20
//       step2  + INTT >= 1
//       step3  + |DCA_xy|, |DCA_z| <= 1.5 cm   (= FINAL)
//     D identical in every bin and species; N and NQ never increase from one
//     step to the next; NQ <= N <= D; step3 identical to the FINAL output.
//  B. validation table (pions; truth |eta| < 1.1, standard binning, in-range
//     bins): true-z / reco-z x MVTX2 / MVTX3 from the validated Step-4A files,
//     FINAL with the same denominator, FINAL with the final denominator in the
//     pion window (|y| < 1.0); plus all species in their analysis windows.
//  C. conditioning check: the MVTX2 -> MVTX3 survival N3/N2 agrees between
//     true-z and reco-z within an absolute tolerance of 0.005 (per charge); the
//     significance (binomial) is printed as information. A statistical test is
//     not used: the true effect (about -0.002, seen in 10k and 1M) is resolved
//     at > 5 sigma with 1M events.
//
// Usage: root -l -b -q 'ValidateFinalTrkEff.C+("FINALDIR","VARIANTSDIR","log.txt")'
//   FINALDIR    output directory of run_trkeff_final.sh
//   VARIANTSDIR Step-4A variants directory of the same sample
// =============================================================================
#include "../src/EffTypes.h"

#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace valfin
{
FILE *gLog = nullptr;
int gN = 0, gFail = 0;
void Out(const std::string &s)
{
  printf("%s\n", s.c_str());
  if (gLog) fprintf(gLog, "%s\n", s.c_str());
}
void Check(bool ok, const std::string &what)
{
  ++gN;
  if (!ok) ++gFail;
  Out(std::string(ok ? "PASS  " : "FAIL  ") + what);
}
const std::vector<const char *> &kSp = efftypes::TrackTags();

TH2D *Get(TFile &f, const char *kind, const char *sp)
{
  auto *h = f.Get<TH2D>(Form("h%s_%s_ypt", kind, sp));
  if (!h) Out(Form("missing h%s_%s_ypt in %s", kind, sp, f.GetName()));
  return h;
}
// in-range sum over y in [ylo, yhi) and pT in [plo, phi)
double Sum(const TH2D *h, double ylo, double yhi, double plo, double phi)
{
  double s = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i)
    for (int j = 1; j <= h->GetNbinsY(); ++j)
    {
      const double y = h->GetXaxis()->GetBinCenter(i), p = h->GetYaxis()->GetBinCenter(j);
      if (y > ylo && y < yhi && p > plo && p < phi) s += h->GetBinContent(i, j);
    }
  return s;
}
struct Eff
{
  double n, d;
  double e() const { return d > 0 ? n / d : NAN; }
  double s() const { return d > 0 ? std::sqrt(e() * (1 - e()) / d) : NAN; }
};
Eff Integrated(TFile &f, const char *sp, double ylo = -99, double yhi = 99, double plo = -1, double phi = 99)
{
  return {Sum(Get(f, "Num", sp), ylo, yhi, plo, phi), Sum(Get(f, "Den", sp), ylo, yhi, plo, phi)};
}
}  // namespace valfin

void ValidateFinalTrkEff(const char *finalDir, const char *variantsDir, const char *logPath)
{
  using namespace valfin;
  gLog = fopen(logPath, "w");
  const std::string F = finalDir, V = variantsDir;
  Out("=== ValidateFinalTrkEff ===");
  Out("final dir   : " + F);
  Out("variants dir: " + V);

  // ---- A. one-cut steps --------------------------------------------------------
  std::vector<std::unique_ptr<TFile>> step;
  for (int s = 0; s <= 3; ++s) step.emplace_back(TFile::Open((F + Form("/TrkEff_step%d.root", s)).c_str()));
  std::unique_ptr<TFile> fin(TFile::Open((F + "/TrkEff_final.root").c_str()));
  bool filesOk = fin && !fin->IsZombie();
  for (auto &f : step) filesOk = filesOk && f && !f->IsZombie();
  Check(filesOk, "A0 step0..step3 and FINAL outputs present");
  if (!filesOk)
  {
    Out(Form("%d checks, %d failed", gN, gFail));
    return;
  }
  for (int s = 0; s <= 3; ++s) Out(Form("  step%d: %s", s, step[s]->Get<TNamed>("selection")->GetTitle()));
  long badD = 0, badMono = 0, badOrder = 0, badFinal = 0, cells = 0;
  for (const char *sp : kSp)
  {
    for (int s = 1; s <= 3; ++s)
    {
      TH2D *d0 = Get(*step[s - 1], "Den", sp), *d1 = Get(*step[s], "Den", sp);
      TH2D *n0 = Get(*step[s - 1], "Num", sp), *n1 = Get(*step[s], "Num", sp);
      TH2D *q0 = Get(*step[s - 1], "NumQ", sp), *q1 = Get(*step[s], "NumQ", sp);
      for (int c = 0; c < d0->GetNcells(); ++c)
      {
        ++cells;
        badD += d0->GetBinContent(c) != d1->GetBinContent(c);
        badMono += n1->GetBinContent(c) > n0->GetBinContent(c) || q1->GetBinContent(c) > q0->GetBinContent(c);
        badOrder += !(q1->GetBinContent(c) <= n1->GetBinContent(c) && n1->GetBinContent(c) <= d1->GetBinContent(c));
      }
    }
    for (const char *k : {"Den", "Num", "NumQ"})
    {
      TH2D *a = Get(*step[3], k, sp), *b = Get(*fin, k, sp);
      for (int c = 0; c < a->GetNcells(); ++c) badFinal += a->GetBinContent(c) != b->GetBinContent(c);
    }
  }
  Check(badD == 0, Form("A1 denominator identical in all steps (%ld cells compared, %ld differ)", cells, badD));
  Check(badMono == 0, Form("A2 N and NQ never increase when a cut is added (%ld violations)", badMono));
  Check(badOrder == 0, Form("A3 NQ <= N <= D in every cell of every step (%ld violations)", badOrder));
  Check(badFinal == 0, Form("A4 step3 identical to the FINAL output (%ld cells differ)", badFinal));
  Out("  integrated, |y| < 1.0 (in range):        step0     step1     step2     step3=FINAL");
  for (const char *sp : kSp)
  {
    std::string line = Form("    %-5s", sp);
    for (int s = 0; s <= 3; ++s) line += Form("   %.4f", Integrated(*step[s], sp, -1.0, 1.0).e());
    Out(line);
  }

  // ---- B. validation table -------------------------------------------------------
  auto open = [&](const std::string &p) { return std::unique_ptr<TFile>(TFile::Open(p.c_str())); };
  auto tz2 = open(V + "/TrkEff_mvtx2_ztrue_incl.root"), tz3 = open(V + "/TrkEff_mvtx3_ztrue_incl.root");
  auto rz2 = open(V + "/TrkEff_mvtx2_zreco_incl.root"), rz3 = open(V + "/TrkEff_mvtx3_zreco_incl.root");
  auto same = open(F + "/TrkEff_final_samedenom.root");
  const bool tabOk = tz2 && tz3 && rz2 && rz3 && same;
  Check(tabOk, "B0 Step-4A variants and the same-denominator FINAL output present");
  if (tabOk)
  {
    Out("");
    Out("  VALIDATION TABLE (pions)                                          pi+              pi-         pi- - pi+");
    struct Row
    {
      const char *label;
      TFile *f;
      double ylo, yhi;
    };
    const Row rows[] = {{"true-z, MVTX >= 2 (TPC 18, truth |eta|<1.1)", tz2.get(), -99, 99},
                        {"true-z, MVTX >= 3 (TPC 18, truth |eta|<1.1)", tz3.get(), -99, 99},
                        {"reco-z, MVTX >= 2 (TPC 18, truth |eta|<1.1)", rz2.get(), -99, 99},
                        {"reco-z, MVTX >= 3 (TPC 18, truth |eta|<1.1)", rz3.get(), -99, 99},
                        {"FINAL reco-z, same denominator (truth |eta|<1.1)", same.get(), -99, 99},
                        {"FINAL reco-z, this output's denominator (|y|<1.0)", fin.get(), -1.0, 1.0}};
    for (const auto &r : rows)
    {
      const Eff p = Integrated(*r.f, "pip", r.ylo, r.yhi), m = Integrated(*r.f, "pim", r.ylo, r.yhi);
      Out(Form("    %-58s %.3f+-%.3f   %.3f+-%.3f   %+.3f+-%.3f", r.label, p.e(), p.s(), m.e(), m.s(), m.e() - p.e(),
               std::hypot(p.s(), m.s())));
    }
    Out("  what each comparison measures:");
    Out("    true-z MVTX2 -> MVTX3      : isolated track-level cost of the third MVTX cluster (no vertex condition)");
    Out("    true-z -> reco-z           : effect of conditioning on the MVTX>=3-based final vertex");
    Out("    reco-z MVTX2 -> MVTX3      : the same track-level cost inside the vertexed sample");
    Out("    reco-z MVTX2 -> FINAL same : extra numerator cuts TPC 20, INTT 1, DCA 1.5 cm");
    Out(std::string("    FINAL same -> this output   : phase space (") + fin->Get<TNamed>("selection")->GetTitle() + "; |y| < 1); not a track-quality effect");
    Out("");
    Out("  FINAL, all species in their analysis windows (true y, true pT):");
    struct Win
    {
      const char *sp;
      double ylo, yhi, plo, phi;
    };
    std::vector<Win> wins;  // EffTypes.h windows: |y| < absMax, ptMin <= pT < ptAnaMax (else the last map edge)
    for (const auto *e : efftypes::MeasuredSpecies(false))
      wins.push_back({e->tag, -e->absMax, e->absMax, e->ptMin, e->ptAnaMax < 50 ? e->ptAnaMax : e->pt.back()});
    for (const auto &w : wins)
    {
      const Eff e = Integrated(*fin, w.sp, w.ylo, w.yhi, w.plo, w.phi);
      Out(Form("    %-5s |y| < %.1f, %.1f <= pT < %4.1f GeV/c : D %6.0f  N %6.0f  eps %.3f +- %.3f", w.sp, w.yhi, w.plo,
               w.phi > 50 ? 5.0 : w.phi, e.d, e.n, e.e(), e.s()));
    }

    // ---- C. conditioning check ---------------------------------------------------
    Out("");
    for (const char *sp : {"pip", "pim"})
    {
      const Eff t2 = Integrated(*tz2, sp), t3 = Integrated(*tz3, sp), r2 = Integrated(*rz2, sp), r3 = Integrated(*rz3, sp);
      const Eff st{t3.n, t2.n}, sr{r3.n, r2.n};  // survival N3/N2
      const double diff = st.e() - sr.e(), sig = std::hypot(st.s(), sr.s());
      constexpr double kCondTol = 0.005;  // absolute tolerance on the survival difference
      Check(std::fabs(diff) < kCondTol,
            Form("C  %s MVTX2->MVTX3 survival true-z %.4f vs reco-z %.4f (diff %+.4f, |diff| < %.3f; %.1f sigma, info)", sp,
                 st.e(), sr.e(), diff, kCondTol, sig > 0 ? diff / sig : 0.));
    }
  }
  Out(Form("%d checks, %d failed", gN, gFail));
  Out(gFail == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  if (gLog) fclose(gLog);
}
