// =============================================================================
// CompareTrkEffVariants.C - comparison plots of tracking-efficiency variants
//
// Input : the outputs of analysis/run_trkeff_variants.sh in VARIANTDIR,
//         TrkEff_mvtx<M>_<ztrue|zreco>_<incl|track>.root (opened read-only).
// Output (OUTDIR):
//   <pdfName>                      page 1: integrated efficiencies of all 12 variants
//                                  pages 2-4: pions, kaons, (anti)protons - per particle
//                                  (positive top row, negative bottom row) three comparisons:
//                                    MVTX >= 0 / 2 / 3          (|z_true| < 8, inclusive)
//                                    inclusive / trackable       (|z_true| < 8, MVTX >= 0)
//                                    true-z / reco-z selection   (MVTX >= 0, inclusive)
//   figures/*.{pdf,png}            every panel as a single figure
//   TrkEffVariants_plots.root      every drawn graph (name g_<variant>_<efficiency>)
//   TrkEffVariants_numbers.txt     integrated D, N, NQ, eps, epsQ for every variant/species
// Statistics as in PlotTrkEff.C: TEfficiency, Clopper-Pearson 68.3%, no point for
// empty bins; every curve shows passed/total and its number of empty bins.
//
// Usage:
//   root -l -b -q 'CompareTrkEffVariants.C+("VARIANTDIR","OUTDIR","sample label","TrkEffVariants.pdf")'
// =============================================================================
#include "TrkEffPlotStyle.h"

#include <TFile.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace trkeffplot;

namespace cmpvar
{
// One efficiency variant = one TrkEffAnalysis output.
struct Variant
{
  std::string tag;    // e.g. "mvtx2_ztrue_incl"
  std::string label;  // legend label within a comparison
  std::string shortLabel;
  int color, marker;
};

// The three comparisons. The first variant of each is the common reference
// (MVTX >= 0, |z_true| < 8 cm, inclusive = the validated baseline).
struct Comparison
{
  const char *name, *fixed;  // fixed = selections common to all curves (printed on the pad)
  std::vector<Variant> variants;
};
const std::vector<Comparison> &Comparisons()
{
  static const std::vector<Comparison> c = {
      {"mvtx", "|z_{true}|<8 cm, inclusive",
       {{"mvtx0_ztrue_incl", "MVTX #geq 0", "MVTX0", kBlack, 20},
        {"mvtx2_ztrue_incl", "MVTX #geq 2", "MVTX2", kAzure + 2, 21},
        {"mvtx3_ztrue_incl", "MVTX #geq 3", "MVTX3", kRed + 1, 22}}},
      {"trackable", "|z_{true}|<8 cm, MVTX#geq0",
       {{"mvtx0_ztrue_incl", "inclusive", "incl.", kBlack, 20},
        {"mvtx0_ztrue_track", "trackable (#geq18 lay.)", "track.", kGreen + 2, 25}}},
      {"vertex", "MVTX#geq0, inclusive",
       {{"mvtx0_ztrue_incl", "|z_{true}| < 8 cm", "z_{true}", kBlack, 20},
        {"mvtx0_zreco_incl", "|z_{reco}| < 8 cm (vertex)", "z_{reco}", kMagenta + 1, 33}}}};
  return c;
}

const char *kAllTags[12] = {"mvtx0_ztrue_incl", "mvtx2_ztrue_incl", "mvtx3_ztrue_incl", "mvtx0_ztrue_track",
                            "mvtx2_ztrue_track", "mvtx3_ztrue_track", "mvtx0_zreco_incl", "mvtx2_zreco_incl",
                            "mvtx3_zreco_incl", "mvtx0_zreco_track", "mvtx2_zreco_track", "mvtx3_zreco_track"};

// All efficiencies of all variants, detached from the (read-only) input files.
struct Inputs
{
  std::map<std::string, TEfficiency *> eff;  // key "<variant>/<eff name>"
  std::map<std::string, std::string> selection;  // variant -> TNamed "selection"
};

Inputs ReadVariants(const std::string &dir)
{
  Inputs in;
  for (const char *tag : kAllTags)
  {
    std::unique_ptr<TFile> f(TFile::Open(Form("%s/TrkEff_%s.root", dir.c_str(), tag), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error(std::string("missing variant ") + tag);
    auto *sel = dynamic_cast<TNamed *>(f->Get("selection"));
    in.selection[tag] = sel ? sel->GetTitle() : "?";
    for (const auto &s : Species())
      for (const char *kind : {"eff", "effQ"})
        for (const char *proj : {"ypt", "pt", "y"})
        {
          const std::string n = Form("%s_%s_%s", kind, s.tag, proj);
          auto *e = dynamic_cast<TEfficiency *>(f->Get(n.c_str()));
          if (!e) throw std::runtime_error("missing " + n + " in " + tag);
          e = static_cast<TEfficiency *>(e->Clone());
          e->SetDirectory(nullptr);
          in.eff[std::string(tag) + "/" + n] = e;
        }
  }
  return in;
}

// One panel: the curves of one comparison for one species versus pT.
std::vector<TGraphAsymmErrors *> DrawComparison(const Inputs &in, const Comparison &c, const SpeciesStyle &s,
                                                const std::string &sample)
{
  std::vector<Curve> cv;
  for (const auto &v : c.variants)
  {
    const std::string key = v.tag + "/eff_" + s.tag + "_pt";
    cv.push_back({in.eff.at(key), v.label, v.shortLabel, v.color, v.marker, 1,
                  "g_" + v.tag + "_eff_" + s.tag + "_pt"});
  }
  auto graphs = DrawEfficiencyCurves(cv, true, sample, true);
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.034);
  t.DrawLatex(gPad->GetLeftMargin() + 0.02, 0.86, Form("#bf{%s}:  %s", s.label, c.fixed));
  return graphs;
}

void WriteNumbers(const std::string &path, const Inputs &in)
{
  FILE *f = fopen(path.c_str(), "w");
  fprintf(f, "# variant  species  generated  reconstructed  reconstructed_correct_charge  eff  eff_lo  eff_up  effQ  effQ_lo  effQ_up   (integrated over the (y,pT) range)\n");
  for (const char *tag : kAllTags)
  {
    fprintf(f, "# %s : %s\n", tag, in.selection.at(tag).c_str());
    for (const auto &s : Species())
    {
      const Integrated a = Integrate(in.eff.at(std::string(tag) + "/eff_" + s.tag + "_ypt"));
      const Integrated q = Integrate(in.eff.at(std::string(tag) + "/effQ_" + s.tag + "_ypt"));
      fprintf(f, "%-18s %-5s %5.0f %5.0f %5.0f  %.4f %.4f %.4f  %.4f %.4f %.4f\n", tag, s.tag, a.total, a.passed, q.passed,
              a.eff, a.low, a.up, q.eff, q.low, q.up);
    }
  }
  fclose(f);
}

// Compact row label of a variant tag, e.g. "mvtx2_zreco_track" -> "|z_reco|<8 (vtx) | MVTX>=2 | trackable".
std::string RowLabel(const std::string &tag)
{
  const bool reco = tag.find("zreco") != std::string::npos, track = tag.find("track") != std::string::npos;
  return std::string(reco ? "|z_{reco}|<8 (vtx)" : "|z_{true}|<8") + "   MVTX#geq" + tag.substr(4, 1) + "   " +
         (track ? "trackable" : "inclusive");
}

// Page 1: integrated efficiency (N/D = eps) of every variant and species.
void TablePage(TCanvas *c, const Inputs &in, const std::string &sample)
{
  c->cd();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.030);
  t.DrawLatex(0.03, 0.95, ("#bf{Tracking-efficiency variants}   #bf{#it{sPHENIX}} Simulation   " + sample).c_str());
  t.SetTextSize(0.022);
  t.DrawLatex(0.03, 0.91, "Integrated over the (y, p_{T}) range: #varepsilon = reconstructed/generated (Clopper-Pearson 68.3%). TPC #geq 18, p_{T} > 0.1 GeV/c, "
                          "|#eta| < 1.1 for all. Trackable: #geq 18 truth TPC layers (conditional efficiency).");
  const double x0 = 0.03, xs[7] = {0.03, 0.30, 0.415, 0.53, 0.645, 0.76, 0.875};
  double y = 0.85;
  t.SetTextSize(0.021);
  t.SetTextFont(62);
  t.DrawLatex(xs[0], y, "event sel.   MVTX   generated set");
  for (size_t k = 0; k < Species().size(); ++k) t.DrawLatex(xs[k + 1], y, Species()[k].label);
  t.SetTextFont(42);
  for (const char *tag : kAllTags)
  {
    y -= 0.062;
    t.DrawLatex(x0, y, RowLabel(tag).c_str());
    for (size_t k = 0; k < Species().size(); ++k)
    {
      const Integrated a = Integrate(in.eff.at(std::string(tag) + "/eff_" + Species()[k].tag + "_ypt"));
      t.DrawLatex(xs[k + 1], y, a.defined() ? Form("%.0f/%.0f = %.3f^{+%.3f}_{-%.3f}", a.passed, a.total, a.eff, a.up - a.eff, a.eff - a.low)
                                            : "undefined");
    }
  }
}
}  // namespace cmpvar

void CompareTrkEffVariants(const char *variantDir, const char *outDir, const char *sampleLabel = "",
                           const char *pdfName = "TrkEffVariants.pdf")
{
  using namespace cmpvar;
  SetTrkEffStyle();
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  const Inputs in = ReadVariants(variantDir);
  const std::string dir = outDir, pdf = dir + "/" + pdfName, sample = sampleLabel;
  gSystem->mkdir((dir + "/figures").c_str(), true);
  WriteNumbers(dir + "/TrkEffVariants_numbers.txt", in);
  std::unique_ptr<TFile> objects(TFile::Open((dir + "/TrkEffVariants_plots.root").c_str(), "RECREATE"));
  std::set<std::string> written;

  // Pass 1: multipage PDF (+ objects); pass 2: single figures. ROOT has one active
  // PDF stream, so the two are never written at the same time.
  const std::vector<std::pair<const char *, std::vector<int>>> groups = {{"pions", {0, 1}}, {"kaons", {2, 3}}, {"protons", {4, 5}}};
  for (const bool pagesPass : {true, false})
  {
    {
      TCanvas c("table", "", 1600, 900);
      TablePage(&c, in, sample);
      if (pagesPass)
      {
        c.Print((pdf + "[").c_str());
        c.Print(pdf.c_str(), "Title:table");
      }
      else
      {
        c.SaveAs(Form("%s/figures/page_table.pdf", dir.c_str()));
        c.SaveAs(Form("%s/figures/page_table.png", dir.c_str()));
      }
    }
    for (const auto &g : groups)
    {
      TCanvas page(Form("compare_%s", g.first), "", 1600, 900);
      page.Divide(3, 2, 0.002, 0.002);
      for (int row = 0; row < 2; ++row)
        for (size_t k = 0; k < Comparisons().size(); ++k)
        {
          const SpeciesStyle &s = Species()[g.second[row]];
          page.cd(row * 3 + k + 1);
          for (TGraphAsymmErrors *gr : DrawComparison(in, Comparisons()[k], s, sample))
            if (pagesPass && written.insert(gr->GetName()).second)
            {
              objects->cd();
              gr->Write();
            }
          if (!pagesPass)
          {
            TCanvas one("one", "", 900, 750);
            DrawComparison(in, Comparisons()[k], s, sample);
            const std::string name = Form("%s/figures/eff_vs_pt_%s_%s", dir.c_str(), Comparisons()[k].name, s.tag);
            one.SaveAs((name + ".pdf").c_str());
            one.SaveAs((name + ".png").c_str());
          }
        }
      if (pagesPass) page.Print(pdf.c_str(), Form("Title:%s", g.first));
      else
      {
        page.SaveAs(Form("%s/figures/page_%s.pdf", dir.c_str(), g.first));
        page.SaveAs(Form("%s/figures/page_%s.png", dir.c_str(), g.first));
      }
    }
    if (pagesPass)
    {
      TCanvas c("c", "", 1600, 900);
      c.Print((pdf + "]").c_str());
    }
  }
  objects->Close();
  printf("CompareTrkEffVariants: 4 pages -> %s\n", pdf.c_str());
}
