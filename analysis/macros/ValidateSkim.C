// =============================================================================
// ValidateSkim.C - independent validation of the EvalSkim output
//
// Re-reads the ORIGINAL SvtxEvaluator files (own code, does not use
// EvalCommon.h / EvalSkim.cc) and compares them with the skim files, job by job.
//
// Usage:
//   root -l -b -q 'ValidateSkim.C+("eval.list","SKIMDIR","report.log",NEXP,NZ8,NZ8VTX,"VtxEff_TPC.root")'
//     eval.list   evaluator files used for the skim
//     SKIMDIR     directory with skim_NNNNNN.root
//     NEXP        expected number of events (-1 = do not check)
//     NZ8,NZ8VTX  expected events / vertexed events with |z_true| < 8 cm (-1 = skip)
//     VtxEff_TPC.root  optional: original VtxEffAnalysis output (TPC selection,
//                 evt tree) for an event-by-event cross-check ("" = skip)
//
// Checks (PASS/FAIL), see the report: A event count, B ntp_info events once,
// C vertexless events, D |z|<8 vertex efficiency, E key uniqueness, F PDG/charge,
// G rapidity, H reco rows vs ntp_track, I matching multiplicities, J NaN/Inf,
// K duplicate jobs/events, X cross-check with the original VtxEff evt tree,
// R truth retention recomputed. DCA comparisons are reported as information.
// =============================================================================
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace vs
{
FILE *gOut = nullptr;
std::map<std::string, std::pair<long, long>> gChk;  // name -> (checked, failed)
std::map<std::string, std::string> gFirstFail;

void Log(const char *s)
{
  printf("%s\n", s);
  if (gOut) fprintf(gOut, "%s\n", s);
}
// count one test of check "name"; ok=false records a failure (first message kept)
void Test(const std::string &name, bool ok, const std::string &msg = "")
{
  auto &c = gChk[name];
  c.first++;
  if (!ok)
  {
    c.second++;
    if (!gFirstFail.count(name)) gFirstFail[name] = msg;
  }
}
bool Fin(double x) { return std::isfinite(x); }
int Q(int pdg)
{
  if (std::abs(pdg) >= 1000000000) return 9999;  // nuclei: not tested here
  TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(pdg);
  return p ? (int)std::lround(p->Charge() / 3.) : 0;
}
double M(int pdg)
{
  if (std::abs(pdg) >= 1000000000) return NAN;
  TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(pdg);
  return p ? p->Mass() : NAN;
}
// independent rapidity formula: 0.5 ln((E+pz)/(E-pz))
double Yref(double px, double py, double pz, double m)
{
  if (!Fin(m)) return NAN;
  const double E = std::sqrt(px * px + py * py + pz * pz + m * m);
  return 0.5 * std::log((E + pz) / (E - pz));
}
int Job(const std::string &p)
{
  std::smatch m;
  static const std::regex re("_([0-9]+)\\.root$");
  return std::regex_search(p, m, re) ? std::stoi(m[1]) : -1;
}
// read selected float branches of a TNtuple
struct Nt
{
  TTree *t;
  std::vector<std::string> n;
  std::vector<float> v;
  Nt(TTree *tt, std::vector<std::string> names) : t(tt), n(names), v(names.size())
  {
    t->SetBranchStatus("*", 0);
    for (size_t i = 0; i < n.size(); ++i)
    {
      t->SetBranchStatus(n[i].c_str(), 1);
      t->SetBranchAddress(n[i].c_str(), &v[i]);
    }
  }
  float operator()(const char *s) const
  {
    for (size_t i = 0; i < n.size(); ++i)
      if (n[i] == s) return v[i];
    printf("no column %s\n", s);
    return NAN;
  }
};
}  // namespace vs

void ValidateSkim(const char *evalList, const char *skimDir, const char *report, long nExp = -1, long nZ8 = -1,
                  long nZ8Vtx = -1, const char *vtxEffRef = "")
{
  using namespace vs;
  gOut = fopen(report, "w");
  Log("=== ValidateSkim ===");
  Log(Form("evaluator list: %s", evalList));
  Log(Form("skim dir      : %s", skimDir));

  // ---- K: job IDs ------------------------------------------------------------------
  std::vector<std::string> files;
  {
    std::ifstream in(evalList);
    std::string l;
    while (std::getline(in, l))
      if (!l.empty()) files.push_back(l);
  }
  std::set<int> jobs;
  for (auto &f : files)
  {
    const int j = Job(f);
    Test("K duplicate job IDs in the input", j >= 0 && jobs.insert(j).second, f);
  }

  // global accumulators
  long nEvtSkim = 0, nEvtInfo = 0, nVtxSkim = 0, z8 = 0, z8vtx = 0, nNoVtx = 0, nRandomSeed = 0;
  std::set<std::pair<int, int>> evKeys;
  std::map<std::pair<int, int>, std::tuple<int, int, int, int, int, int, int, float, float, int>> evSummary;
  long dcaN[3] = {0, 0, 0}, dcaRefAttached = 0, dcaRefAttachedReco = 0;
  std::vector<double> dXYtrue, dZtrue, dXYreco, dZreco;
  long nTruthRows = 0, nRecoRows = 0, nYnan = 0, nYnanNucleus = 0, nYnanUnknown = 0;
  std::map<int, long> speciesQ;

  for (auto &f : files)
  {
    const int job = Job(f);
    TFile *E = TFile::Open(f.c_str());
    TFile *S = TFile::Open(Form("%s/skim_%06d.root", skimDir, job));
    Test("K one skim file per job", S && !S->IsZombie(), f);
    if (!E || !S || S->IsZombie()) continue;
    auto *tE = (TTree *)S->Get("EventTree");
    auto *tT = (TTree *)S->Get("TruthTree");
    auto *tR = (TTree *)S->Get("RecoTree");

    // --- evaluator: info, gpoint, vertex --------------------------------------------
    std::map<int, int> infoRows, infoSeed, infoNtrk;
    {
      Nt n((TTree *)E->Get("ntp_info"), {"event", "seed", "ntrk"});
      for (long i = 0; i < n.t->GetEntries(); ++i)
      {
        n.t->GetEntry(i);
        int e = std::lround(n("event"));
        infoRows[e]++;
        // random-seed mode (PHRandomSeed from std::random_device): ntp_info.seed is NaN and the
        // skim stores -1 ("no fixed seed")
        const float sd = n("seed");
        infoSeed[e] = std::isfinite(sd) ? static_cast<int>(std::lround(sd)) : -1;
        nRandomSeed += !std::isfinite(sd);
        infoNtrk[e] = std::lround(n("ntrk"));
      }
    }
    nEvtInfo += infoRows.size();
    std::map<int, std::vector<float>> gp;
    {
      Nt n((TTree *)E->Get("ntp_gpoint"), {"event", "gvx", "gvy", "gvz"});
      for (long i = 0; i < n.t->GetEntries(); ++i)
      {
        n.t->GetEntry(i);
        gp[std::lround(n("event"))] = {n("gvx"), n("gvy"), n("gvz")};
      }
    }
    std::map<int, int> nFiniteVtx;
    std::map<int, float> bestVz;
    std::map<int, int> bestNtr;
    {
      Nt n((TTree *)E->Get("ntp_vertex"), {"event", "vz", "vx", "vy", "ntracks"});
      for (long i = 0; i < n.t->GetEntries(); ++i)
      {
        n.t->GetEntry(i);
        int e = std::lround(n("event"));
        if (!(Fin(n("vx")) && Fin(n("vy")) && Fin(n("vz")))) continue;
        nFiniteVtx[e]++;
        int nt = Fin(n("ntracks")) ? std::lround(n("ntracks")) : 0;
        if (!bestNtr.count(e) || nt > bestNtr[e])
        {
          bestNtr[e] = nt;
          bestVz[e] = n("vz");
        }
      }
    }

    // --- evaluator gtrack (all) and track (all) ------------------------------------
    struct GR
    {
      int pdg, prim, parent, tid;
      float pt, eta, px, py, pz;
    };
    std::map<std::pair<int, int>, GR> gtr;  // (event, gtrackID)
    {
      Nt n((TTree *)E->Get("ntp_gtrack"),
           {"event", "gtrackID", "gflavor", "gprimary", "gparentid", "trackID", "gpt", "geta", "gpx", "gpy", "gpz"});
      for (long i = 0; i < n.t->GetEntries(); ++i)
      {
        n.t->GetEntry(i);
        GR g{(int)std::lround(n("gflavor")), (int)std::lround(n("gprimary")), (int)std::lround(n("gparentid")),
             Fin(n("trackID")) ? (int)std::lround(n("trackID")) : -1, n("gpt"), n("geta"), n("gpx"), n("gpy"), n("gpz")};
        gtr[{(int)std::lround(n("event")), (int)std::lround(n("gtrackID"))}] = g;
      }
    }
    struct TR
    {
      float pt, eta, q, ntpc, nmaps, dca3dxy, dca3dz, pcaz, sm;
      int gid;
    };
    std::map<std::tuple<int, int>, TR> trk;      // (event, trackID)
    std::map<std::pair<int, int>, int> nMatchEval;  // (event, gtrackID) -> #tracks
    std::map<int, int> nTrkEval;
    {
      Nt n((TTree *)E->Get("ntp_track"), {"event", "trackID", "pt", "eta", "charge", "ntpc", "nmaps", "dca3dxy",
                                          "dca3dz", "pcaz", "singlematch", "gtrackID"});
      for (long i = 0; i < n.t->GetEntries(); ++i)
      {
        n.t->GetEntry(i);
        int e = std::lround(n("event"));
        int gid = Fin(n("gtrackID")) ? (int)std::lround(n("gtrackID")) : 0;
        TR r{n("pt"), n("eta"), n("charge"), n("ntpc"), n("nmaps"), n("dca3dxy"), n("dca3dz"), n("pcaz"), n("singlematch"), gid};
        bool fresh = trk.emplace(std::make_tuple(e, (int)std::lround(n("trackID"))), r).second;
        Test("H evaluator (event,trackID) unique", fresh);
        nTrkEval[e]++;
        if (gid != 0) nMatchEval[{e, gid}]++;
      }
    }

    // --- expected truth retention (recomputed independently, default skim.cfg) -----
    std::set<std::pair<int, int>> expKeep, v0;
    for (auto &kv : gtr)
    {
      const GR &g = kv.second;
      if (g.prim != 1) continue;
      const double y = Yref(g.px, g.py, g.pz, M(g.pdg));
      const bool broad = g.pt >= 0. && (std::fabs(g.eta) < 1.5 || (Fin(y) && std::fabs(y) < 1.5));
      if (Q(g.pdg) != 0 && broad) expKeep.insert(kv.first);
      if (broad && (g.pdg == 310 || std::abs(g.pdg) == 3122))
      {
        expKeep.insert(kv.first);
        v0.insert(kv.first);
      }
    }
    for (auto &kv : gtr)
      if (kv.second.prim != 1 && v0.count({kv.first.first, kv.second.parent})) expKeep.insert(kv.first);
    for (auto &kv : nMatchEval) expKeep.insert(kv.first);

    // --- EventTree ------------------------------------------------------------------
    std::map<int, int> seenEvt;
    {
      int job_, event, seed, vtxId;
      short hasVtx, nVtx, nTrkAll, nTrkGood, m1, m2, m3, Nch, Npip, Npim, NKp, NKm, Np, Npbar;
      float xT, yT, zT, xR, yR, zR, chi2;
      tE->SetBranchAddress("job", &job_);
      tE->SetBranchAddress("event", &event);
      tE->SetBranchAddress("seed", &seed);
      tE->SetBranchAddress("vtxId", &vtxId);
      tE->SetBranchAddress("hasVtx", &hasVtx);
      tE->SetBranchAddress("nVtx", &nVtx);
      tE->SetBranchAddress("nTrkAll", &nTrkAll);
      tE->SetBranchAddress("nTrkGood", &nTrkGood);
      tE->SetBranchAddress("nGoodMVTX1", &m1);
      tE->SetBranchAddress("nGoodMVTX2", &m2);
      tE->SetBranchAddress("nGoodMVTX3", &m3);
      tE->SetBranchAddress("NchPrim", &Nch);
      tE->SetBranchAddress("Npip", &Npip);
      tE->SetBranchAddress("Npim", &Npim);
      tE->SetBranchAddress("NKp", &NKp);
      tE->SetBranchAddress("NKm", &NKm);
      tE->SetBranchAddress("Np", &Np);
      tE->SetBranchAddress("Npbar", &Npbar);
      tE->SetBranchAddress("xTrue", &xT);
      tE->SetBranchAddress("yTrue", &yT);
      tE->SetBranchAddress("zTrue", &zT);
      tE->SetBranchAddress("xReco", &xR);
      tE->SetBranchAddress("yReco", &yR);
      tE->SetBranchAddress("zReco", &zR);
      tE->SetBranchAddress("vtxChi2", &chi2);
      for (long i = 0; i < tE->GetEntries(); ++i)
      {
        tE->GetEntry(i);
        nEvtSkim++;
        Test("K EventTree job == file job", job_ == job);
        Test("K no duplicate (job,event)", evKeys.insert({job_, event}).second, Form("job %d event %d", job_, event));
        seenEvt[event]++;
        Test("B event is an ntp_info event", infoRows.count(event) == 1, Form("job %d event %d", job, event));
        Test("B seed equals ntp_info.seed (NaN -> -1)", infoSeed[event] == seed);
        const bool evVtx = nFiniteVtx[event] > 0;
        Test("C hasVtx equals evaluator (finite ntp_vertex row)", (hasVtx != 0) == evVtx, Form("job %d ev %d", job, event));
        Test("C nVtx equals finite ntp_vertex rows", nVtx == nFiniteVtx[event]);
        if (!hasVtx)
        {
          nNoVtx++;
          Test("C vertexless: reco vertex and chi2 are NaN, vtxId = -1",
               !Fin(xR) && !Fin(yR) && !Fin(zR) && !Fin(chi2) && vtxId == -1);
        }
        else
        {
          nVtxSkim++;
          Test("J vertexed: reco vertex finite", Fin(xR) && Fin(yR) && Fin(zR));
          Test("C zReco equals evaluator vertex with most tracks", zR == bestVz[event]);
        }
        Test("J true vertex finite", Fin(xT) && Fin(yT) && Fin(zT));
        Test("D zTrue equals ntp_gpoint.gvz", gp.count(event) && gp[event][2] == zT && gp[event][0] == xT && gp[event][1] == yT);
        Test("H nTrkAll == ntp_track rows == ntp_info.ntrk", nTrkAll == nTrkEval[event] && nTrkAll == infoNtrk[event],
             Form("job %d ev %d: %d %d %d", job, event, nTrkAll, nTrkEval[event], infoNtrk[event]));
        if (Fin(zT) && std::fabs(zT) < 8)
        {
          z8++;
          if (hasVtx) z8vtx++;
        }
        evSummary[{job_, event}] = std::make_tuple((int)nTrkGood, (int)Nch, (int)m1, (int)m2, (int)m3, (int)hasVtx,
                                                   (int)nVtx, zT, zR, Npip * 1000000 + Npim * 10000 + NKp * 100 + NKm);
        (void)Np;
        (void)Npbar;
      }
      for (auto &kv : infoRows) Test("B every ntp_info event appears exactly once", seenEvt[kv.first] == 1, Form("job %d ev %d", job, kv.first));
    }

    // --- TruthTree ------------------------------------------------------------------
    std::set<std::pair<int, int>> skimKeep;
    std::map<std::pair<int, int>, int> truthNReco;
    {
      int job_, event, id, pdg, bestTid;
      short q, isPrim, keep, nReco, pdgKnown;
      float vx, vy, vz, px, py, pz, pt, eta, phi, y, mass;
      tT->SetBranchAddress("job", &job_);
      tT->SetBranchAddress("event", &event);
      tT->SetBranchAddress("gtrackID", &id);
      tT->SetBranchAddress("pdg", &pdg);
      tT->SetBranchAddress("q", &q);
      tT->SetBranchAddress("pdgKnown", &pdgKnown);
      tT->SetBranchAddress("isPrimary", &isPrim);
      tT->SetBranchAddress("keep", &keep);
      tT->SetBranchAddress("bestTrackID", &bestTid);
      tT->SetBranchAddress("nRecoMatches", &nReco);
      tT->SetBranchAddress("vx", &vx);
      tT->SetBranchAddress("vy", &vy);
      tT->SetBranchAddress("vz", &vz);
      tT->SetBranchAddress("px", &px);
      tT->SetBranchAddress("py", &py);
      tT->SetBranchAddress("pz", &pz);
      tT->SetBranchAddress("pt", &pt);
      tT->SetBranchAddress("eta", &eta);
      tT->SetBranchAddress("phi", &phi);
      tT->SetBranchAddress("y", &y);
      tT->SetBranchAddress("mass", &mass);
      for (long i = 0; i < tT->GetEntries(); ++i)
      {
        tT->GetEntry(i);
        nTruthRows++;
        const std::pair<int, int> k{event, id};
        Test("E truth key (job,event,gtrackID) unique", skimKeep.insert(k).second, Form("job %d ev %d id %d", job, event, id));
        Test("E gtrackID != 0 and |gtrackID| < 2^24", id != 0 && std::abs(id) < (1 << 24));
        auto it = gtr.find(k);
        Test("E truth key exists in ntp_gtrack", it != gtr.end());
        if (it == gtr.end()) continue;
        const GR &g = it->second;
        Test("F pdg equals rounded ntp_gtrack.gflavor", pdg == g.pdg);
        Test("F isPrimary equals ntp_gtrack.gprimary", isPrim == g.prim);
        const int qr = Q(pdg);
        if (qr != 9999) Test("F charge from PDG (TDatabasePDG, independent)", q == qr, Form("pdg %d q %d ref %d", pdg, q, qr));
        speciesQ[pdg] = q;
        if (std::abs(pdg) < 1000000000 && TDatabasePDG::Instance()->GetParticle(pdg))
        {
          Test("G mass equals PDG mass", std::fabs(mass - M(pdg)) < 1e-5);
          const double yr = Yref(px, py, pz, M(pdg));
          Test("G rapidity equals 0.5 ln((E+pz)/(E-pz)) with PDG mass", Fin(y) && std::fabs(y - yr) < 2e-4 * std::max(1.0, std::fabs(yr)),
               Form("pdg %d y %g ref %g", pdg, y, yr));
          Test("G |y| <= |eta| (massive particles)", !(mass > 0) || std::fabs(y) <= std::fabs(eta) + 1e-4);
        }
        if (!Fin(y))
        {
          nYnan++;
          if (std::abs(pdg) >= 1000000000) nYnanNucleus++;
          else if (!pdgKnown) nYnanUnknown++;
        }
        Test("J truth kinematics and production point finite",
             Fin(vx) && Fin(vy) && Fin(vz) && Fin(px) && Fin(py) && Fin(pz) && Fin(pt) && Fin(eta) && Fin(phi));
        Test("I bestTrackID equals ntp_gtrack.trackID", bestTid == g.tid);
        const int ne = nMatchEval.count(k) ? nMatchEval[k] : 0;
        Test("I nRecoMatches equals number of ntp_track rows with this gtrackID", nReco == ne, Form("ev %d id %d: %d vs %d", event, id, nReco, ne));
        truthNReco[k] = nReco;
        // retention bits
        const bool prim = g.prim == 1;
        const double yy = Yref(g.px, g.py, g.pz, M(g.pdg));
        const bool broad = g.pt >= 0. && (std::fabs(g.eta) < 1.5 || (Fin(yy) && std::fabs(yy) < 1.5));
        Test("R keep bit 1 correct", ((keep & 1) != 0) == (prim && Q(g.pdg) != 0 && Q(g.pdg) != 9999 && broad) || std::abs(g.pdg) >= 1000000000);
        Test("R keep bit 8 correct", ((keep & 8) != 0) == (ne > 0));
      }
      Test("R retained particle set equals independent recomputation", skimKeep == expKeep,
           Form("job %d: skim %zu expected %zu", job, skimKeep.size(), expKeep.size()));
    }

    // --- RecoTree -------------------------------------------------------------------
    {
      int job_, event, tid, gid, vertexID;
      short charge, ntpc, nmaps, matched, sm, best, nSame, rank, chargeOK, gq, dcaRef, truthKept;
      float pt, eta, px, py, pz, phi, dxy, dz, pcax, pcay, pcaz, evx, evy, evz, xT, zT, xR, zR, linxyT, linzT, linxyR,
          linzR, gpt, geta, quality;
      tR->SetBranchAddress("job", &job_);
      tR->SetBranchAddress("event", &event);
      tR->SetBranchAddress("trackID", &tid);
      tR->SetBranchAddress("gtrackID", &gid);
      tR->SetBranchAddress("vertexID", &vertexID);
      tR->SetBranchAddress("charge", &charge);
      tR->SetBranchAddress("ntpc", &ntpc);
      tR->SetBranchAddress("nmaps", &nmaps);
      tR->SetBranchAddress("matched", &matched);
      tR->SetBranchAddress("singlematch", &sm);
      tR->SetBranchAddress("isBestOfTruth", &best);
      tR->SetBranchAddress("nSameTruth", &nSame);
      tR->SetBranchAddress("rankSameTruth", &rank);
      tR->SetBranchAddress("chargeOK", &chargeOK);
      tR->SetBranchAddress("gq", &gq);
      tR->SetBranchAddress("dcaRef", &dcaRef);
      tR->SetBranchAddress("truthKept", &truthKept);
      tR->SetBranchAddress("pt", &pt);
      tR->SetBranchAddress("eta", &eta);
      tR->SetBranchAddress("px", &px);
      tR->SetBranchAddress("py", &py);
      tR->SetBranchAddress("pz", &pz);
      tR->SetBranchAddress("phi", &phi);
      tR->SetBranchAddress("quality", &quality);
      tR->SetBranchAddress("dca3dxy", &dxy);
      tR->SetBranchAddress("dca3dz", &dz);
      tR->SetBranchAddress("pcax", &pcax);
      tR->SetBranchAddress("pcay", &pcay);
      tR->SetBranchAddress("pcaz", &pcaz);
      tR->SetBranchAddress("evalVtxX", &evx);
      tR->SetBranchAddress("evalVtxY", &evy);
      tR->SetBranchAddress("evalVtxZ", &evz);
      tR->SetBranchAddress("xTrue", &xT);
      tR->SetBranchAddress("zTrue", &zT);
      tR->SetBranchAddress("xReco", &xR);
      tR->SetBranchAddress("zReco", &zR);
      tR->SetBranchAddress("linXYTrue", &linxyT);
      tR->SetBranchAddress("linZTrue", &linzT);
      tR->SetBranchAddress("linXYReco", &linxyR);
      tR->SetBranchAddress("linZReco", &linzR);
      tR->SetBranchAddress("gpt", &gpt);
      tR->SetBranchAddress("geta", &geta);
      std::set<std::pair<int, int>> keys;
      std::map<std::pair<int, int>, std::vector<int>> ranks;
      std::map<std::pair<int, int>, int> nSameSeen;
      for (long i = 0; i < tR->GetEntries(); ++i)
      {
        tR->GetEntry(i);
        nRecoRows++;
        Test("E reco key (job,event,trackID) unique", keys.insert({event, tid}).second);
        auto it = trk.find(std::make_tuple(event, tid));
        Test("H every RecoTree row is an ntp_track row", it != trk.end());
        if (it == trk.end()) continue;
        const TR &r = it->second;
        Test("H reco values identical to ntp_track (pt,eta,charge,ntpc,nmaps,dca3dxy,dca3dz,pcaz,gtrackID)",
             pt == r.pt && eta == r.eta && charge == std::lround(r.q) && ntpc == std::lround(r.ntpc) &&
                 nmaps == std::lround(r.nmaps) && dxy == r.dca3dxy && dz == r.dca3dz && pcaz == r.pcaz && gid == r.gid);
        Test("I isBestOfTruth equals evaluator singlematch", best == std::lround(r.sm) && sm == std::lround(r.sm));
        Test("J reco kinematics, PCA, DCA, quality finite",
             Fin(pt) && Fin(eta) && Fin(px) && Fin(py) && Fin(pz) && Fin(phi) && Fin(pcax) && Fin(pcay) && Fin(pcaz) &&
                 Fin(dxy) && Fin(dz) && Fin(quality));
        Test("J linXYTrue/linZTrue finite", Fin(linxyT) && Fin(linzT));
        const bool evVtx = nFiniteVtx[event] > 0;
        Test("J reco-vertex columns NaN exactly when the event has no vertex",
             evVtx == (Fin(xR) && Fin(zR) && Fin(linxyR) && Fin(linzR)));
        if (matched)
        {
          const std::pair<int, int> k{event, gid};
          ranks[k].push_back(rank);
          nSameSeen[k] = nSame;
          Test("I nSameTruth equals evaluator count of tracks with this gtrackID", nSame == nMatchEval[k]);
          Test("I matched truth particle retained in TruthTree (keep_matched)", truthKept == 1 && skimKeep.count(k) == 1);
          Test("F gq equals charge from PDG of the matched particle", gq == Q(gtr[k].pdg) || std::abs(gtr[k].pdg) >= 1000000000);
          Test("F chargeOK == (charge == gq)", chargeOK == (charge == gq ? 1 : 0));
          Test("J matched true kinematics finite", Fin(gpt) && Fin(geta));
        }
        else
          Test("I unmatched track has gtrackID 0 and chargeOK -1", gid == 0 && chargeOK == -1);
        // DCA reference (information + consistency)
        dcaN[dcaRef]++;
        if (vertexID >= 0)
        {
          dcaRefAttached++;
          if (dcaRef == 1) dcaRefAttachedReco++;
        }
        if (dcaRef == 2)
        {
          dXYtrue.push_back(std::fabs(std::fabs(linxyT) - std::fabs(dxy)));
          dZtrue.push_back(std::fabs(std::fabs(linzT) - std::fabs(dz)));
        }
        if (dcaRef == 1)
        {
          dXYreco.push_back(std::fabs(std::fabs(linxyR) - std::fabs(dxy)));
          dZreco.push_back(std::fabs(std::fabs(linzR) - std::fabs(dz)));
        }
        Test("J dcaRef is 1 (reco vertex) or 2 (true vertex)", dcaRef == 1 || dcaRef == 2);
        Test("J evalVtx equals the vertex named by dcaRef",
             (dcaRef == 2 && evx == xT && evz == zT) || (dcaRef == 1 && Fin(evx) && Fin(evz)));
      }
      Test("H RecoTree rows == ntp_track rows (file)", (long)keys.size() == (long)trk.size(),
           Form("job %d: %zu vs %zu", job, keys.size(), trk.size()));
      for (auto &kv : ranks)
      {
        auto v = kv.second;
        std::sort(v.begin(), v.end());
        bool perm = (int)v.size() == nSameSeen[kv.first];
        for (size_t j = 0; j < v.size(); ++j) perm = perm && v[j] == (int)j + 1;
        Test("I rankSameTruth is 1..nSameTruth for every matched particle", perm);
        Test("I TruthTree nRecoMatches equals RecoTree rows of that particle", truthNReco[kv.first] == (int)kv.second.size());
      }
    }
    E->Close();
    S->Close();
  }

  // ---- global checks ------------------------------------------------------------------
  Test("A EventTree events == ntp_info events", nEvtSkim == nEvtInfo, Form("%ld vs %ld", nEvtSkim, nEvtInfo));
  if (nExp >= 0) Test(Form("A EventTree events == %ld", nExp), nEvtSkim == nExp, Form("%ld", nEvtSkim));
  Test("C vertexless events present", nNoVtx > 0, Form("%ld", nNoVtx));
  if (nZ8 >= 0) Test(Form("D |z_true| < 8 cm: events == %ld", nZ8), z8 == nZ8, Form("%ld", z8));
  if (nZ8Vtx >= 0) Test(Form("D |z_true| < 8 cm: vertexed == %ld", nZ8Vtx), z8vtx == nZ8Vtx, Form("%ld", z8vtx));
  // species charge table
  const int sp[] = {211, -211, 321, -321, 2212, -2212, 11, -11, 13, -13, 22, 310, 3122, -3122};
  const int qs[] = {1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0};
  for (int i = 0; i < 14; ++i)
    if (speciesQ.count(sp[i])) Test("F charge table (pi,K,p,e,mu,gamma,K0S,Lambda)", speciesQ[sp[i]] == qs[i], Form("pdg %d", sp[i]));

  // ---- X: cross-check with the original VtxEffAnalysis evt tree ------------------------
  if (vtxEffRef && *vtxEffRef)
  {
    TFile *V = TFile::Open(vtxEffRef);
    TTree *t = V ? (TTree *)V->Get("evt") : nullptr;
    Test("X original VtxEff evt tree readable", t != nullptr, vtxEffRef);
    if (t)
    {
      int job, event, n, N, hasVtx, nVtx, m1, m2, m3, Npip, Npim, NKp, NKm;
      float zT, zR;
      t->SetBranchAddress("job", &job);
      t->SetBranchAddress("event", &event);
      t->SetBranchAddress("n", &n);
      t->SetBranchAddress("N", &N);
      t->SetBranchAddress("hasVtx", &hasVtx);
      t->SetBranchAddress("nVtx", &nVtx);
      t->SetBranchAddress("nGoodMVTX", &m1);
      t->SetBranchAddress("nMVTX2", &m2);
      t->SetBranchAddress("nMVTX3", &m3);
      t->SetBranchAddress("Npip", &Npip);
      t->SetBranchAddress("Npim", &Npim);
      t->SetBranchAddress("NKp", &NKp);
      t->SetBranchAddress("NKm", &NKm);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);
      Test("X same number of events as original evt tree", t->GetEntries() == nEvtSkim, Form("%lld vs %ld", t->GetEntries(), nEvtSkim));
      for (long i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        auto it = evSummary.find({job, event});
        Test("X event present in skim", it != evSummary.end());
        if (it == evSummary.end()) continue;
        auto &s = it->second;
        Test("X n, N, nMVTX1/2/3, hasVtx, nVtx identical to original evt tree",
             std::get<0>(s) == n && std::get<1>(s) == N && std::get<2>(s) == m1 && std::get<3>(s) == m2 &&
                 std::get<4>(s) == m3 && std::get<5>(s) == hasVtx && std::get<6>(s) == nVtx,
             Form("job %d ev %d", job, event));
        Test("X zTrue, zReco identical (NaN-aware)", std::get<7>(s) == zT && (std::get<8>(s) == zR || (!Fin(zR) && !Fin(std::get<8>(s)))));
        Test("X species counts pi+-,K+- identical", std::get<9>(s) == Npip * 1000000 + Npim * 10000 + NKp * 100 + NKm);
      }
    }
  }

  // ---- report ---------------------------------------------------------------------
  Log("");
  Log("--- checks (tests / failures) ---");
  long nFailChecks = 0;
  for (auto &kv : gChk)
  {
    const bool ok = kv.second.second == 0;
    if (!ok) nFailChecks++;
    Log(Form("%s  %-95s %8ld tests %6ld failed%s", ok ? "PASS" : "FAIL", kv.first.c_str(), kv.second.first, kv.second.second,
             ok ? "" : Form("  (first: %s)", gFirstFail[kv.first].c_str())));
  }
  Log("");
  Log("--- information ---");
  Log(Form("events %ld (vertexed %ld, vertexless %ld); |z_true|<8: %ld events, %ld vertexed, eps_vtx = %.4f", nEvtSkim, nVtxSkim,
           nNoVtx, z8, z8vtx, z8 ? double(z8vtx) / z8 : 0.));
  Log(Form("events in random-seed mode (ntp_info.seed NaN): %ld of %ld", nRandomSeed, nEvtInfo));
  Log(Form("TruthTree rows %ld (%.1f / event), RecoTree rows %ld (%.1f / event)", nTruthRows, double(nTruthRows) / nEvtSkim,
           nRecoRows, double(nRecoRows) / nEvtSkim));
  Log(Form("TruthTree rows with y = NaN: %ld (nuclei %ld, unknown PDG %ld)", nYnan, nYnanNucleus, nYnanUnknown));
  Log(Form("evaluator DCA reference: reco vertex %ld tracks, TRUE vertex %ld tracks, other %ld; tracks with vertexID>=0: %ld, of which ref = reco vertex %ld",
           dcaN[1], dcaN[2], dcaN[0], dcaRefAttached, dcaRefAttachedReco));
  auto q = [](std::vector<double> v, double f) -> double {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
  };
  auto frac = [](const std::vector<double> &v, double c) {
    long k = 0;
    for (double x : v) k += x < c;
    return v.empty() ? NAN : double(k) / v.size();
  };
  Log(Form("straight-line approximation vs evaluator DCA, same vertex (|lin| - |eval|):"));
  Log(Form("  ref TRUE vertex (%zu tracks): xy median %.2e cm, 90%% %.2e cm, frac < 10 um %.3f ; z median %.2e cm, 90%% %.2e cm, frac < 10 um %.3f",
           dXYtrue.size(), q(dXYtrue, .5), q(dXYtrue, .9), frac(dXYtrue, 1e-3), q(dZtrue, .5), q(dZtrue, .9), frac(dZtrue, 1e-3)));
  Log(Form("  ref RECO vertex (%zu tracks): xy median %.2e cm, 90%% %.2e cm, frac < 10 um %.3f ; z median %.2e cm, 90%% %.2e cm, frac < 10 um %.3f",
           dXYreco.size(), q(dXYreco, .5), q(dXYreco, .9), frac(dXYreco, 1e-3), q(dZreco, .5), q(dZreco, .9), frac(dZreco, 1e-3)));
  Log("");
  Log(Form("SUMMARY: %zu checks, %ld failed", gChk.size(), nFailChecks));
  Log(nFailChecks == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED");
  fclose(gOut);
  gOut = nullptr;
}
