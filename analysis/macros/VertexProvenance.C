// =============================================================================
// VertexProvenance.C - the event vertex used by TrkEffAnalysis (use_reco_zvtx = 1)
// is the final tracking vertex of the DST (SvtxVertexMap)
//
// For each DST (one per job), event by event (DST entry i = skim event i):
//   - number of SvtxVertexMap vertices      == skim EventTree nVtx (finite rows)
//   - a vertex exists                        == skim hasVtx
//   - vertex with the most tracks: z, ntracks == skim zReco, vtxNtracks (on a tie
//     for the most tracks any tied vertex is accepted: DST map order and evaluator
//     row order can differ; the number of tied events is printed)
//   - every SvtxVertexMap vertex has >= 2 tracks (the finder builds vertices from
//     track pairs)
// SiliconSvtxVertexMap (silicon-only finder, used only by the TPC crossing
// finder) is counted for information.
// Read-only: the DST is opened with TFile; nothing is written to it.
//
// Usage: root -l -b -q 'VertexProvenance.C+("dst.list","SKIMDIR","log.txt"[,maxFiles])'
// =============================================================================
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <regex>
#include <string>

void VertexProvenance(const char *dstList, const char *skimDir, const char *logPath, int maxFiles = 0)
{
  FILE *log = fopen(logPath, "w");
  auto out = [&](const std::string &s) {
    printf("%s\n", s.c_str());
    if (log) fprintf(log, "%s\n", s.c_str());
  };
  out("=== VertexProvenance ===");
  out(std::string("DST list: ") + dstList + "\nskims   : " + skimDir);
  long nEv = 0, badN = 0, badHas = 0, badZ = 0, badNtrk = 0, smallVtx = 0, nVtxTot = 0, nSilVtx = 0, nFiles = 0, badFiles = 0, tiedEv = 0;
  std::ifstream in(dstList);
  std::string dst;
  while (std::getline(in, dst) && (maxFiles <= 0 || nFiles < maxFiles))
  {
    if (dst.empty()) continue;
    std::smatch m;
    if (!std::regex_search(dst, m, std::regex("_(\\d{6})\\.root$")))
    {
      out("cannot read the job number from " + dst);
      ++badFiles;
      continue;
    }
    const std::string job = m[1];
    TFile fd(dst.c_str(), "READ"), fs((std::string(skimDir) + "/skim_" + job + ".root").c_str(), "READ");
    auto *T = fd.Get<TTree>("T");
    auto *E = fs.Get<TTree>("EventTree");
    if (!T || !E || T->GetEntries() != E->GetEntries())
    {
      out("missing tree or entry-count mismatch for job " + job);
      ++badFiles;
      continue;
    }
    ++nFiles;
    SvtxVertexMap *vmap = nullptr, *smap = nullptr;
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("DST#SVTX#SvtxVertexMap*", 1);
    T->SetBranchStatus("DST#SVTX#SiliconSvtxVertexMap*", 1);
    T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vmap);
    T->SetBranchAddress("DST#SVTX#SiliconSvtxVertexMap", &smap);
    short hasVtx, nVtx, vtxNtracks;
    float zReco;
    E->SetBranchAddress("hasVtx", &hasVtx);
    E->SetBranchAddress("nVtx", &nVtx);
    E->SetBranchAddress("vtxNtracks", &vtxNtracks);
    E->SetBranchAddress("zReco", &zReco);
    for (long i = 0; i < T->GetEntries(); ++i)
    {
      T->GetEntry(i);
      E->GetEntry(i);
      ++nEv;
      const int n = vmap ? static_cast<int>(vmap->size()) : 0;
      nVtxTot += n;
      nSilVtx += smap ? smap->size() : 0;
      const SvtxVertex *best = nullptr;
      for (const auto &[key, v] : *vmap)
      {
        if (v->size_tracks() < 2) ++smallVtx;
        if (!best || v->size_tracks() > best->size_tracks()) best = v;  // first on ties, as EvalSkim
      }
      badN += n != nVtx;
      badHas += (n > 0) != (hasVtx != 0);
      if (best)
      {
        // ties (several vertices with the maximal ntracks): the DST map order (key) and the evaluator
        // row order (used by EvalSkim) can differ, so any tied vertex is accepted
        bool zMatch = false;
        int nTied = 0;
        for (const auto &[key, v] : *vmap)
          if (v->size_tracks() == best->size_tracks())
          {
            ++nTied;
            zMatch |= std::isfinite(zReco) && std::fabs(v->get_z() - zReco) < 1e-4;
          }
        tiedEv += nTied > 1;
        badZ += !zMatch;
        badNtrk += static_cast<int>(best->size_tracks()) != vtxNtracks;
      }
    }
  }
  auto line = [&](bool ok, const std::string &s) { out(std::string(ok ? "PASS  " : "FAIL  ") + s); };
  out(Form("files %ld, events %ld; SvtxVertexMap vertices %ld; SiliconSvtxVertexMap vertices %ld (info)", nFiles, nEv,
           nVtxTot, nSilVtx));
  line(badFiles == 0, Form("all DSTs matched to a skim with the same number of events (%ld problems)", badFiles));
  line(badN == 0, Form("number of SvtxVertexMap vertices == skim nVtx in every event (%ld differ)", badN));
  line(badHas == 0, Form("vertex present in the DST == skim hasVtx (%ld differ)", badHas));
  line(badZ == 0 && badNtrk == 0,
       Form("vertex with the most tracks (any of them on a tie): z and ntracks == skim zReco, vtxNtracks (%ld z, %ld ntracks differ)",
            badZ, badNtrk));
  out(Form("events with a tie for the most tracks: %ld (info; event vertex = first tied row of the evaluator)", tiedEv));
  line(smallVtx == 0, Form("every SvtxVertexMap vertex has >= 2 tracks (%ld with fewer)", smallVtx));
  const bool ok = badFiles == 0 && badN == 0 && badHas == 0 && badZ == 0 && badNtrk == 0 && smallVtx == 0 && nFiles > 0;
  out(ok ? "VALIDATION PASSED" : "VALIDATION FAILED");
  if (log) fclose(log);
}
