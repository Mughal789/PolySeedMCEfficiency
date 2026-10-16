// =============================================================================
// EvalSkim.cc - compact skim of SvtxEvaluator files (step 1 of the unified
// vertex + single-track efficiency framework)
//
// Input : SvtxEvaluator ntuples outputeval_NNNNNN.root (one file per job).
// Output: one file per input job, OUTDIR/skim_NNNNNN.root, with three trees.
//         No efficiency is computed here and no analysis cut is applied to the
//         trees; cuts are applied later on the skim.
//
//   EventTree  one row per simulated collision = every ntp_info row, also
//              collisions without a reconstructed vertex (reco vertex = NaN).
//   TruthTree  one row per RETAINED G4 particle (see "Truth retention").
//   RecoTree   one row per reconstructed track = every ntp_track row.
//
// Unique keys:  event  (job, event)
//               truth  (job, event, gtrackID)   - G4 ids restart every event
//               track  (job, event, trackID)
//   job = number in the file name; event = evaluator event counter (restarts
//   at 0 in every file); seed = ntp_info.seed (100000 + job for our
//   productions) is stored as an additional check.
//
// Truth retention (TruthTree.keep is a bit mask; a particle is written if any
// bit is set, all other G4 particles are discarded):
//   bit 1  charged G4 primary (gprimary == 1, charge from PDG != 0) with
//          gpt >= truth_pt_min AND (|geta| < truth_abseta_max OR
//          |y| < truth_absy_max). Defaults: no pT cut, |eta| < 1.5 OR |y| < 1.5.
//   bit 2  G4 primary K0S (310), Lambda (3122) or anti-Lambda (-3122) in the
//          same broad acceptance (keep_v0 = 1).
//   bit 4  direct daughter (gparentid) of a bit-2 particle, any kinematics,
//          any charge (keep_v0 = 1).
//   bit 8  any G4 particle (primary or secondary) that is the dominant truth
//          particle of at least one reconstructed track (keep_matched = 1),
//          so every matched track has its truth particle in the TruthTree.
//   Discarded: neutral primaries other than K0S/Lambda/anti-Lambda, charged
//   primaries outside the broad acceptance (unless matched), and all G4
//   secondaries that are neither V0 daughters nor matched to a track.
//   Note: "G4 primary" = particle handed to GEANT4 by PYTHIA. With
//   ParticleDecays:limitCylinder (1 nm), K0S/Lambda/Sigma+-/Xi/Omega are G4
//   primaries and their decay products are G4 secondaries. Lambda from Xi/Omega
//   decays are G4 secondaries and are therefore NOT bit 2.
//
// Truth matching (from SvtxTruthRecoTableEval, see TRKEFF_PLAN):
//   weight = number of clusters on the track to which the particle contributed;
//   no purity threshold; ties go to the lowest id.
//   reco -> truth : ntp_track.gtrackID = dominant particle of the track
//                   (RecoTree.gtrackID; 0 = no match)
//   truth -> reco : ntp_gtrack.trackID = track with most clusters of the
//                   particle (TruthTree.bestTrackID; -1 = none)
//   The two directions are not always symmetric. Both are kept, plus
//   nRecoMatches (tracks whose dominant particle is this one), and per track
//   nSameTruth, rankSameTruth, isBestOfTruth, singlematch, chargeOK.
//
// DCA: the evaluator values (TrackAnalysisUtils::get_dca) are computed from the
//   track state position pca = (track->get_x,y,z) and the vertex
//   evalVtx = GlobalVertexMap entry of track->get_vertex_id():
//     dca3dxy = component of (pca - evalVtx) transverse to the momentum (xy plane)
//     dca3dz  = pcaz - evalVtxZ   (NO extrapolation along the track)
//   evalVtx is the reconstructed vertex for tracks attached to it
//   (vertexID >= 0, dcaRef = 1) and the TRUE collision vertex for all other
//   tracks, also in vertexed events (dcaRef = 2). For attached tracks the state
//   is at the transverse closest approach to that vertex. For the others it is
//   at the transverse closest approach to the ORIGIN (0,0) - not to the beam
//   spot (-0.07, 0.141) cm and not to the true vertex - so their dca3dz is NOT
//   the z distance at the closest approach to the true vertex (K1: differs from
//   linZTrue by 0.6 mm on average). All original values are kept unchanged.
//   Additional APPROXIMATE variables lin{XY,Z}{True,Reco}: straight line through
//   pca along the reconstructed momentum, evaluated at the transverse closest
//   approach to the true / reconstructed (most tracks) vertex. linXY equals
//   dca3dxy exactly when the vertex is evalVtx; linZ is extrapolated. This is
//   NOT the helix DCA: the neglected curvature term is ~ d^2 / (2R) for a
//   transverse distance d between pca and the vertex, with
//   R[cm] ~ 333 pT[GeV/c] / B[T] (e.g. ~3 um for d = 1.6 mm, pT = 0.1 GeV/c);
//   for tracks far from the vertex (secondaries) it is only indicative.
//
// Usage: EvalSkim -c skim.cfg -o OUTDIR [-l list] [-n maxEventsPerFile] [files...]
// Build: analysis/build.sh
// =============================================================================

#include "EvalCommon.h"

#include <TFile.h>
#include <TNamed.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <memory>

using namespace evalcommon;

namespace
{
const char *kVersion = "EvalSkim v1 (2026-10-15)";

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------
struct SkimConfig
{
  std::string input_dir, input_list, output_dir;
  double truth_pt_min = 0., truth_abseta_max = 1.5, truth_absy_max = 1.5;
  bool keep_v0 = true, keep_matched = true;
  // event-level summary multiplicities (same definitions as VtxEffAnalysis)
  double good_pt_min = 0.1, good_abseta_max = 1.1;
  int good_ntpc_min = 18;
  double true_pt_min = 0.1, true_abseta_max = 1.1;
  long max_events_per_file = 0;
  int compression = 505;
  std::string dump;
};

SkimConfig ReadSkimConfig(const std::string &path)
{
  KeyValueConfig k;
  k.Declare("input_dir", "");
  k.Declare("input_list", "");
  k.Declare("output_dir", "");
  k.Declare("truth_pt_min", "0.0");
  k.Declare("truth_abseta_max", "1.5");
  k.Declare("truth_absy_max", "1.5");
  k.Declare("keep_v0", "1");
  k.Declare("keep_matched", "1");
  k.Declare("good_pt_min", "0.1");
  k.Declare("good_abseta_max", "1.1");
  k.Declare("good_ntpc_min", "18");
  k.Declare("true_pt_min", "0.1");
  k.Declare("true_abseta_max", "1.1");
  k.Declare("max_events_per_file", "0");
  k.Declare("compression", "505");
  k.Read(path);
  SkimConfig c;
  c.input_dir = k.Str("input_dir");
  c.input_list = k.Str("input_list");
  c.output_dir = k.Str("output_dir");
  c.truth_pt_min = k.Dbl("truth_pt_min");
  c.truth_abseta_max = k.Dbl("truth_abseta_max");
  c.truth_absy_max = k.Dbl("truth_absy_max");
  c.keep_v0 = k.Bool("keep_v0");
  c.keep_matched = k.Bool("keep_matched");
  c.good_pt_min = k.Dbl("good_pt_min");
  c.good_abseta_max = k.Dbl("good_abseta_max");
  c.good_ntpc_min = k.Int("good_ntpc_min");
  c.true_pt_min = k.Dbl("true_pt_min");
  c.true_abseta_max = k.Dbl("true_abseta_max");
  c.max_events_per_file = k.Int("max_events_per_file");
  c.compression = k.Int("compression");
  c.dump = k.Dump();
  return c;
}

// float -> integer with a sentinel for NaN
int ToInt(float x, int missing) { return std::isfinite(x) ? static_cast<int>(std::lround(x)) : missing; }
short ToShort(float x, short missing = -1) { return std::isfinite(x) ? static_cast<short>(std::lround(x)) : missing; }

// Straight-line "DCA" to vertex V: line through the PCA along the reconstructed
// momentum. dxy is signed (z component of (P - V) x p_hat_T); dz is the z
// distance at the transverse closest approach. NaN if any input is not finite.
void LinDca(float px, float py, float pz, float pcax, float pcay, float pcaz, float vx, float vy, float vz, float &dxy,
            float &dz)
{
  dxy = dz = NAN;
  const float in[] = {px, py, pz, pcax, pcay, pcaz, vx, vy, vz};
  for (float v : in)
    if (!std::isfinite(v)) return;
  const double pt = std::hypot(px, py);
  if (pt <= 0) return;
  const double ux = px / pt, uy = py / pt;
  const double dx = pcax - vx, dy = pcay - vy, dzz = pcaz - vz;
  dxy = static_cast<float>(dx * uy - dy * ux);
  const double s = -(dx * ux + dy * uy);  // transverse path length to the closest approach
  dz = static_cast<float>(dzz + s * pz / pt);
}

// -----------------------------------------------------------------------------
// In-memory rows of one file
// -----------------------------------------------------------------------------
struct VtxRow
{
  int id;
  float x, y, z, ntracks, chi2, ndof, crossing;
};

struct EventRec
{
  int nInfoRows = 0, nGpointRows = 0;
  int seed = -1, ntrkInfo = -1, gntrkallInfo = -1;
  float xT = NAN, yT = NAN, zT = NAN, tT = NAN;
  std::vector<VtxRow> vtx;  // finite ntp_vertex rows
  int nVtxRowsNaN = 0;
  int best = -1;            // index into vtx: most tracks (first on ties)
  std::vector<int> gRows;   // indices into the file's gtrack rows
  std::vector<int> tRows;   // indices into the file's track rows
};

// ntp_gtrack columns used
const std::vector<std::string> kGCols = {
    "event", "gtrackID", "gflavor", "gprimary", "gpx", "gpy", "gpz", "gpt", "geta", "gphi", "gvx", "gvy", "gvz", "gvt",
    "gnhits", "gnmaps", "gnintt", "gntpc", "gnmms", "gnlmaps", "gnlintt", "gnltpc", "gnlmms", "gparentid",
    "gparentflavor", "gprimaryid", "gprimaryflavor", "trackID"};
// ntp_track columns used
const std::vector<std::string> kTCols = {
    "event", "trackID", "crossing", "px", "py", "pz", "pt", "eta", "phi", "charge", "quality", "chisq", "ndf", "nhits",
    "nmaps", "nintt", "ntpc", "nmms", "nlmaps", "nlintt", "nltpc", "nlmms", "layers", "vertexID", "vx", "vy", "vz",
    "dca2d", "dca2dsigma", "dca3dxy", "dca3dxysigma", "dca3dz", "dca3dzsigma", "pcax", "pcay", "pcaz", "gtrackID",
    "singlematch", "gflavor", "gprimary", "gparentid", "gparentflavor", "gprimaryid", "gprimaryflavor", "gpx", "gpy",
    "gpz", "gpt", "geta", "gphi", "gvx", "gvy", "gvz", "nfromtruth", "nwrong", "ntrumaps", "nwrongmaps", "ntruintt",
    "nwrongintt", "ntrutpc", "nwrongtpc", "ntrumms", "nwrongmms", "layersfromtruth"};

// Read all rows of a TNtuple into a flat float array (row-major).
std::vector<float> ReadAll(TTree *t, const std::vector<std::string> &cols, long &nrows)
{
  FloatRow r(t, cols);
  nrows = t->GetEntries();
  std::vector<float> out(static_cast<size_t>(nrows) * cols.size());
  for (long i = 0; i < nrows; ++i)
  {
    t->GetEntry(i);
    for (size_t c = 0; c < cols.size(); ++c) out[i * cols.size() + c] = r[static_cast<int>(c)];
  }
  return out;
}

int Col(const std::vector<std::string> &cols, const char *n)
{
  auto it = std::find(cols.begin(), cols.end(), std::string(n));
  if (it == cols.end())
  {
    std::cerr << "FATAL: no column " << n << std::endl;
    std::exit(2);
  }
  return static_cast<int>(it - cols.begin());
}

// -----------------------------------------------------------------------------
// Output trees
// -----------------------------------------------------------------------------
struct EventOut
{
  Int_t job, event, seed;
  Float_t xTrue, yTrue, zTrue, tTrue;
  Short_t hasVtx, nVtx, nVtxRowsNaN;
  Int_t vtxId;
  Float_t xReco, yReco, zReco, vtxChi2, vtxNdf, vtxCrossing;
  Short_t vtxNtracks;
  Short_t ntrkInfo, gntrkallInfo;
  Short_t NprimAll, NchPrimAll, NchPrim, Npip, Npim, NKp, NKm, Np, Npbar;
  Short_t nTrkAll, nTrkGood, nGoodMVTX1, nGoodMVTX2, nGoodMVTX3;
  Short_t nTruthKept, nTrkMatched;

  void Branches(TTree *t)
  {
    t->Branch("job", &job, "job/I");
    t->Branch("event", &event, "event/I");
    t->Branch("seed", &seed, "seed/I");
    t->Branch("xTrue", &xTrue, "xTrue/F");
    t->Branch("yTrue", &yTrue, "yTrue/F");
    t->Branch("zTrue", &zTrue, "zTrue/F");
    t->Branch("tTrue", &tTrue, "tTrue/F");
    t->Branch("hasVtx", &hasVtx, "hasVtx/S");
    t->Branch("nVtx", &nVtx, "nVtx/S");
    t->Branch("nVtxRowsNaN", &nVtxRowsNaN, "nVtxRowsNaN/S");
    t->Branch("vtxId", &vtxId, "vtxId/I");
    t->Branch("xReco", &xReco, "xReco/F");
    t->Branch("yReco", &yReco, "yReco/F");
    t->Branch("zReco", &zReco, "zReco/F");
    t->Branch("vtxChi2", &vtxChi2, "vtxChi2/F");
    t->Branch("vtxNdf", &vtxNdf, "vtxNdf/F");
    t->Branch("vtxCrossing", &vtxCrossing, "vtxCrossing/F");
    t->Branch("vtxNtracks", &vtxNtracks, "vtxNtracks/S");
    t->Branch("ntrkInfo", &ntrkInfo, "ntrkInfo/S");
    t->Branch("gntrkallInfo", &gntrkallInfo, "gntrkallInfo/S");
    t->Branch("NprimAll", &NprimAll, "NprimAll/S");
    t->Branch("NchPrimAll", &NchPrimAll, "NchPrimAll/S");
    t->Branch("NchPrim", &NchPrim, "NchPrim/S");
    t->Branch("Npip", &Npip, "Npip/S");
    t->Branch("Npim", &Npim, "Npim/S");
    t->Branch("NKp", &NKp, "NKp/S");
    t->Branch("NKm", &NKm, "NKm/S");
    t->Branch("Np", &Np, "Np/S");
    t->Branch("Npbar", &Npbar, "Npbar/S");
    t->Branch("nTrkAll", &nTrkAll, "nTrkAll/S");
    t->Branch("nTrkGood", &nTrkGood, "nTrkGood/S");
    t->Branch("nGoodMVTX1", &nGoodMVTX1, "nGoodMVTX1/S");
    t->Branch("nGoodMVTX2", &nGoodMVTX2, "nGoodMVTX2/S");
    t->Branch("nGoodMVTX3", &nGoodMVTX3, "nGoodMVTX3/S");
    t->Branch("nTruthKept", &nTruthKept, "nTruthKept/S");
    t->Branch("nTrkMatched", &nTrkMatched, "nTrkMatched/S");
  }
};

struct TruthOut
{
  Int_t job, event, gtrackID, pdg;
  Short_t q, isPrimary, keep, pdgKnown;
  Int_t parentId, parentPdg, primaryId, primaryPdg;
  Float_t vx, vy, vz, vt, px, py, pz, pt, eta, phi, y, mass;
  Short_t gnhits, gnmaps, gnintt, gntpc, gnmms, gnlmaps, gnlintt, gnltpc, gnlmms;
  Int_t bestTrackID;
  Short_t nRecoMatches;

  void Branches(TTree *t)
  {
    t->Branch("job", &job, "job/I");
    t->Branch("event", &event, "event/I");
    t->Branch("gtrackID", &gtrackID, "gtrackID/I");
    t->Branch("pdg", &pdg, "pdg/I");
    t->Branch("q", &q, "q/S");
    t->Branch("pdgKnown", &pdgKnown, "pdgKnown/S");
    t->Branch("isPrimary", &isPrimary, "isPrimary/S");
    t->Branch("keep", &keep, "keep/S");
    t->Branch("parentId", &parentId, "parentId/I");
    t->Branch("parentPdg", &parentPdg, "parentPdg/I");
    t->Branch("primaryId", &primaryId, "primaryId/I");
    t->Branch("primaryPdg", &primaryPdg, "primaryPdg/I");
    t->Branch("vx", &vx, "vx/F");
    t->Branch("vy", &vy, "vy/F");
    t->Branch("vz", &vz, "vz/F");
    t->Branch("vt", &vt, "vt/F");
    t->Branch("px", &px, "px/F");
    t->Branch("py", &py, "py/F");
    t->Branch("pz", &pz, "pz/F");
    t->Branch("pt", &pt, "pt/F");
    t->Branch("eta", &eta, "eta/F");
    t->Branch("phi", &phi, "phi/F");
    t->Branch("y", &y, "y/F");
    t->Branch("mass", &mass, "mass/F");
    t->Branch("gnhits", &gnhits, "gnhits/S");
    t->Branch("gnmaps", &gnmaps, "gnmaps/S");
    t->Branch("gnintt", &gnintt, "gnintt/S");
    t->Branch("gntpc", &gntpc, "gntpc/S");
    t->Branch("gnmms", &gnmms, "gnmms/S");
    t->Branch("gnlmaps", &gnlmaps, "gnlmaps/S");
    t->Branch("gnlintt", &gnlintt, "gnlintt/S");
    t->Branch("gnltpc", &gnltpc, "gnltpc/S");
    t->Branch("gnlmms", &gnlmms, "gnlmms/S");
    t->Branch("bestTrackID", &bestTrackID, "bestTrackID/I");
    t->Branch("nRecoMatches", &nRecoMatches, "nRecoMatches/S");
  }
};

struct RecoOut
{
  Int_t job, event, trackID;
  Float_t px, py, pz, pt, eta, phi;
  Short_t charge;
  Float_t quality, chisq, ndf, crossing;
  Short_t nhits, nmaps, nintt, ntpc, nmms, nlmaps, nlintt, nltpc, nlmms, layers;
  Int_t vertexID;
  Float_t evalVtxX, evalVtxY, evalVtxZ;
  Short_t dcaRef;
  Float_t dca2d, dca2dsigma, dca3dxy, dca3dxysigma, dca3dz, dca3dzsigma, pcax, pcay, pcaz;
  Float_t xTrue, yTrue, zTrue, xReco, yReco, zReco;
  Float_t linXYTrue, linZTrue, linXYReco, linZReco;
  // truth match
  Int_t gtrackID;
  Short_t matched, singlematch, isBestOfTruth, nSameTruth, rankSameTruth, chargeOK, truthKept;
  Int_t gpdg;
  Short_t gq, gprimary;
  Int_t gparentId, gparentPdg, gprimaryId, gprimaryPdg;
  Float_t gpx, gpy, gpz, gpt, geta, gphi, gy, gvx, gvy, gvz;
  Short_t nfromtruth, nwrong, ntrumaps, nwrongmaps, ntruintt, nwrongintt, ntrutpc, nwrongtpc, ntrumms, nwrongmms,
      layersfromtruth;
  Float_t purity;

  void Branches(TTree *t)
  {
    t->Branch("job", &job, "job/I");
    t->Branch("event", &event, "event/I");
    t->Branch("trackID", &trackID, "trackID/I");
    t->Branch("px", &px, "px/F");
    t->Branch("py", &py, "py/F");
    t->Branch("pz", &pz, "pz/F");
    t->Branch("pt", &pt, "pt/F");
    t->Branch("eta", &eta, "eta/F");
    t->Branch("phi", &phi, "phi/F");
    t->Branch("charge", &charge, "charge/S");
    t->Branch("quality", &quality, "quality/F");
    t->Branch("chisq", &chisq, "chisq/F");
    t->Branch("ndf", &ndf, "ndf/F");
    t->Branch("crossing", &crossing, "crossing/F");
    t->Branch("nhits", &nhits, "nhits/S");
    t->Branch("nmaps", &nmaps, "nmaps/S");
    t->Branch("nintt", &nintt, "nintt/S");
    t->Branch("ntpc", &ntpc, "ntpc/S");
    t->Branch("nmms", &nmms, "nmms/S");
    t->Branch("nlmaps", &nlmaps, "nlmaps/S");
    t->Branch("nlintt", &nlintt, "nlintt/S");
    t->Branch("nltpc", &nltpc, "nltpc/S");
    t->Branch("nlmms", &nlmms, "nlmms/S");
    t->Branch("layers", &layers, "layers/S");
    t->Branch("vertexID", &vertexID, "vertexID/I");
    t->Branch("evalVtxX", &evalVtxX, "evalVtxX/F");
    t->Branch("evalVtxY", &evalVtxY, "evalVtxY/F");
    t->Branch("evalVtxZ", &evalVtxZ, "evalVtxZ/F");
    t->Branch("dcaRef", &dcaRef, "dcaRef/S");
    t->Branch("dca2d", &dca2d, "dca2d/F");
    t->Branch("dca2dsigma", &dca2dsigma, "dca2dsigma/F");
    t->Branch("dca3dxy", &dca3dxy, "dca3dxy/F");
    t->Branch("dca3dxysigma", &dca3dxysigma, "dca3dxysigma/F");
    t->Branch("dca3dz", &dca3dz, "dca3dz/F");
    t->Branch("dca3dzsigma", &dca3dzsigma, "dca3dzsigma/F");
    t->Branch("pcax", &pcax, "pcax/F");
    t->Branch("pcay", &pcay, "pcay/F");
    t->Branch("pcaz", &pcaz, "pcaz/F");
    t->Branch("xTrue", &xTrue, "xTrue/F");
    t->Branch("yTrue", &yTrue, "yTrue/F");
    t->Branch("zTrue", &zTrue, "zTrue/F");
    t->Branch("xReco", &xReco, "xReco/F");
    t->Branch("yReco", &yReco, "yReco/F");
    t->Branch("zReco", &zReco, "zReco/F");
    t->Branch("linXYTrue", &linXYTrue, "linXYTrue/F");
    t->Branch("linZTrue", &linZTrue, "linZTrue/F");
    t->Branch("linXYReco", &linXYReco, "linXYReco/F");
    t->Branch("linZReco", &linZReco, "linZReco/F");
    t->Branch("gtrackID", &gtrackID, "gtrackID/I");
    t->Branch("matched", &matched, "matched/S");
    t->Branch("singlematch", &singlematch, "singlematch/S");
    t->Branch("isBestOfTruth", &isBestOfTruth, "isBestOfTruth/S");
    t->Branch("nSameTruth", &nSameTruth, "nSameTruth/S");
    t->Branch("rankSameTruth", &rankSameTruth, "rankSameTruth/S");
    t->Branch("chargeOK", &chargeOK, "chargeOK/S");
    t->Branch("truthKept", &truthKept, "truthKept/S");
    t->Branch("gpdg", &gpdg, "gpdg/I");
    t->Branch("gq", &gq, "gq/S");
    t->Branch("gprimary", &gprimary, "gprimary/S");
    t->Branch("gparentId", &gparentId, "gparentId/I");
    t->Branch("gparentPdg", &gparentPdg, "gparentPdg/I");
    t->Branch("gprimaryId", &gprimaryId, "gprimaryId/I");
    t->Branch("gprimaryPdg", &gprimaryPdg, "gprimaryPdg/I");
    t->Branch("gpx", &gpx, "gpx/F");
    t->Branch("gpy", &gpy, "gpy/F");
    t->Branch("gpz", &gpz, "gpz/F");
    t->Branch("gpt", &gpt, "gpt/F");
    t->Branch("geta", &geta, "geta/F");
    t->Branch("gphi", &gphi, "gphi/F");
    t->Branch("gy", &gy, "gy/F");
    t->Branch("gvx", &gvx, "gvx/F");
    t->Branch("gvy", &gvy, "gvy/F");
    t->Branch("gvz", &gvz, "gvz/F");
    t->Branch("nfromtruth", &nfromtruth, "nfromtruth/S");
    t->Branch("nwrong", &nwrong, "nwrong/S");
    t->Branch("ntrumaps", &ntrumaps, "ntrumaps/S");
    t->Branch("nwrongmaps", &nwrongmaps, "nwrongmaps/S");
    t->Branch("ntruintt", &ntruintt, "ntruintt/S");
    t->Branch("nwrongintt", &nwrongintt, "nwrongintt/S");
    t->Branch("ntrutpc", &ntrutpc, "ntrutpc/S");
    t->Branch("nwrongtpc", &nwrongtpc, "nwrongtpc/S");
    t->Branch("ntrumms", &ntrumms, "ntrumms/S");
    t->Branch("nwrongmms", &nwrongmms, "nwrongmms/S");
    t->Branch("layersfromtruth", &layersfromtruth, "layersfromtruth/S");
    t->Branch("purity", &purity, "purity/F");
  }
};

struct FileResult
{
  int job = -1;
  long nEvents = 0, nVertexed = 0, nTruthRows = 0, nTruthKept = 0, nTracks = 0;
  long keepBits[4] = {0, 0, 0, 0};
  int errors = 0;
  double seconds = 0;
  std::map<std::string, long> warn;
};

// -----------------------------------------------------------------------------
// Process one evaluator file -> one skim file
// -----------------------------------------------------------------------------
FileResult ProcessFile(const std::string &path, const std::string &outPath, const SkimConfig &cfg, PdgLookup &pdgdb)
{
  FileResult R;
  R.job = JobFromName(path);
  TStopwatch sw;
  sw.Start();

  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    R.errors++;
    return R;
  }
  // TFile::Get returns the highest key cycle (the complete tree).
  auto *tInfo = dynamic_cast<TTree *>(f->Get("ntp_info"));
  auto *tGpoint = dynamic_cast<TTree *>(f->Get("ntp_gpoint"));
  auto *tVertex = dynamic_cast<TTree *>(f->Get("ntp_vertex"));
  auto *tTrack = dynamic_cast<TTree *>(f->Get("ntp_track"));
  auto *tGtrack = dynamic_cast<TTree *>(f->Get("ntp_gtrack"));
  if (!tInfo || !tGpoint || !tVertex || !tTrack || !tGtrack)
  {
    std::cerr << "ERROR: missing evaluator tree in " << path << std::endl;
    R.errors++;
    return R;
  }

  // ---- 1. ntp_info: the list of simulated collisions -------------------------
  const long nEvtFile = tInfo->GetEntries();
  const long nEvt = (cfg.max_events_per_file > 0 && cfg.max_events_per_file < nEvtFile) ? cfg.max_events_per_file : nEvtFile;
  std::vector<EventRec> ev(nEvt);
  {
    float event, seed, ntrk, gntrkall;
    Bind(tInfo, {{"event", &event}, {"seed", &seed}, {"ntrk", &ntrk}, {"gntrkall", &gntrkall}});
    for (long i = 0; i < nEvtFile; ++i)
    {
      tInfo->GetEntry(i);
      const long e = std::lround(event);
      if (e < 0 || e >= nEvtFile)
      {
        std::cerr << "ERROR: ntp_info event " << event << " outside [0," << nEvtFile << ")" << std::endl;
        R.errors++;
        continue;
      }
      if (e >= nEvt) continue;
      ev[e].nInfoRows++;
      ev[e].seed = ToInt(seed, -1);
      ev[e].ntrkInfo = ToInt(ntrk, -1);
      ev[e].gntrkallInfo = ToInt(gntrkall, -1);
    }
  }
  auto lookup = [&](float event, const char *tree) -> EventRec * {
    const long e = std::lround(event);
    if (e >= nEvt && e < nEvtFile) return nullptr;  // beyond the test limit
    if (e < 0 || e >= nEvt || ev[e].nInfoRows == 0)
    {
      R.warn[std::string("event not in ntp_info (") + tree + ")"]++;
      R.errors++;
      return nullptr;
    }
    return &ev[e];
  };

  // ---- 2. ntp_gpoint: true collision vertex ------------------------------------
  {
    float event, gvx, gvy, gvz, gvt;
    Bind(tGpoint, {{"event", &event}, {"gvx", &gvx}, {"gvy", &gvy}, {"gvz", &gvz}, {"gvt", &gvt}});
    for (long i = 0; i < tGpoint->GetEntries(); ++i)
    {
      tGpoint->GetEntry(i);
      EventRec *r = lookup(event, "ntp_gpoint");
      if (!r) continue;
      r->nGpointRows++;
      r->xT = gvx;
      r->yT = gvy;
      r->zT = gvz;
      r->tT = gvt;
    }
  }

  // ---- 3. ntp_vertex: reconstructed vertices (all-NaN row = no vertex) -------
  {
    float event, vid, vx, vy, vz, ntracks, chi2, ndof, crossing;
    Bind(tVertex, {{"event", &event}, {"vertexID", &vid}, {"vx", &vx}, {"vy", &vy}, {"vz", &vz}, {"ntracks", &ntracks},
                   {"chi2", &chi2}, {"ndof", &ndof}, {"crossing", &crossing}});
    for (long i = 0; i < tVertex->GetEntries(); ++i)
    {
      tVertex->GetEntry(i);
      EventRec *r = lookup(event, "ntp_vertex");
      if (!r) continue;
      if (!(Finite(vx) && Finite(vy) && Finite(vz)))
      {
        r->nVtxRowsNaN++;
        continue;
      }
      r->vtx.push_back({ToInt(vid, -1), vx, vy, vz, ntracks, chi2, ndof, crossing});
      const int k = static_cast<int>(r->vtx.size()) - 1;
      if (r->best < 0 || ToInt(ntracks, 0) > ToInt(r->vtx[r->best].ntracks, 0)) r->best = k;
    }
  }

  // ---- 4. ntp_gtrack (all G4 particles) and ntp_track (all tracks) -----------
  long nG = 0, nT = 0;
  const std::vector<float> G = ReadAll(tGtrack, kGCols, nG);
  const std::vector<float> T = ReadAll(tTrack, kTCols, nT);
  const size_t NG = kGCols.size(), NT = kTCols.size();
  auto g = [&](long row, int col) { return G[row * NG + col]; };
  auto t = [&](long row, int col) { return T[row * NT + col]; };
  R.nTruthRows = nG;
  R.nTracks = nT;

  // column indices
  const int gEv = Col(kGCols, "event"), gId = Col(kGCols, "gtrackID"), gFl = Col(kGCols, "gflavor"),
            gPr = Col(kGCols, "gprimary"), gPx = Col(kGCols, "gpx"), gPy = Col(kGCols, "gpy"), gPz = Col(kGCols, "gpz"),
            gPt = Col(kGCols, "gpt"), gEta = Col(kGCols, "geta"), gPhi = Col(kGCols, "gphi"), gVx = Col(kGCols, "gvx"),
            gVy = Col(kGCols, "gvy"), gVz = Col(kGCols, "gvz"), gVt = Col(kGCols, "gvt"),
            gNhits = Col(kGCols, "gnhits"), gNmaps = Col(kGCols, "gnmaps"), gNintt = Col(kGCols, "gnintt"),
            gNtpc = Col(kGCols, "gntpc"), gNmms = Col(kGCols, "gnmms"), gNlmaps = Col(kGCols, "gnlmaps"),
            gNlintt = Col(kGCols, "gnlintt"), gNltpc = Col(kGCols, "gnltpc"), gNlmms = Col(kGCols, "gnlmms"),
            gPar = Col(kGCols, "gparentid"), gParFl = Col(kGCols, "gparentflavor"), gPri = Col(kGCols, "gprimaryid"),
            gPriFl = Col(kGCols, "gprimaryflavor"), gTid = Col(kGCols, "trackID");
  auto tc = [&](const char *n) { return Col(kTCols, n); };

  for (long i = 0; i < nG; ++i)
  {
    EventRec *r = lookup(g(i, gEv), "ntp_gtrack");
    if (r) r->gRows.push_back(static_cast<int>(i));
  }
  for (long i = 0; i < nT; ++i)
  {
    EventRec *r = lookup(t(i, tc("event")), "ntp_track");
    if (r) r->tRows.push_back(static_cast<int>(i));
  }

  // ---- 5. output file and trees ---------------------------------------------
  TFile out(outPath.c_str(), "RECREATE", "", cfg.compression);
  if (out.IsZombie())
  {
    std::cerr << "ERROR: cannot create " << outPath << std::endl;
    R.errors++;
    return R;
  }
  auto *tE = new TTree("EventTree", "One row per simulated collision (job, event); reco vertex NaN if none");
  auto *tTr = new TTree("TruthTree", "Retained G4 particles, key (job, event, gtrackID); see 'retention'");
  auto *tR = new TTree("RecoTree", "All reconstructed tracks (ntp_track), key (job, event, trackID); no cuts");
  EventOut E{};
  TruthOut P{};
  RecoOut K{};
  E.Branches(tE);
  P.Branches(tTr);
  K.Branches(tR);

  // track column indices (resolved once)
  const int tTid = tc("trackID"), tCross = tc("crossing"), tPx = tc("px"), tPy = tc("py"), tPz = tc("pz"),
            tPt = tc("pt"), tEta = tc("eta"), tPhi = tc("phi"), tQ = tc("charge"), tQual = tc("quality"),
            tChi = tc("chisq"), tNdf = tc("ndf"), tNhits = tc("nhits"), tNmaps = tc("nmaps"), tNintt = tc("nintt"),
            tNtpc = tc("ntpc"), tNmms = tc("nmms"), tNlmaps = tc("nlmaps"), tNlintt = tc("nlintt"),
            tNltpc = tc("nltpc"), tNlmms = tc("nlmms"), tLayers = tc("layers"), tVid = tc("vertexID"),
            tVx = tc("vx"), tVy = tc("vy"), tVz = tc("vz"), tD2 = tc("dca2d"), tD2s = tc("dca2dsigma"),
            tDxy = tc("dca3dxy"), tDxys = tc("dca3dxysigma"), tDz = tc("dca3dz"), tDzs = tc("dca3dzsigma"),
            tPcax = tc("pcax"), tPcay = tc("pcay"), tPcaz = tc("pcaz"), tGid = tc("gtrackID"),
            tSm = tc("singlematch"), tGfl = tc("gflavor"), tGpr = tc("gprimary"), tGpar = tc("gparentid"),
            tGparFl = tc("gparentflavor"), tGpri = tc("gprimaryid"), tGpriFl = tc("gprimaryflavor"),
            tGpx = tc("gpx"), tGpy = tc("gpy"), tGpz = tc("gpz"), tGpt = tc("gpt"), tGeta = tc("geta"),
            tGphi = tc("gphi"), tGvx = tc("gvx"), tGvy = tc("gvy"), tGvz = tc("gvz"), tNft = tc("nfromtruth"),
            tNw = tc("nwrong"), tNtm = tc("ntrumaps"), tNwm = tc("nwrongmaps"), tNti = tc("ntruintt"),
            tNwi = tc("nwrongintt"), tNtt = tc("ntrutpc"), tNwt = tc("nwrongtpc"), tNtmm = tc("ntrumms"),
            tNwmm = tc("nwrongmms"), tLft = tc("layersfromtruth");

  for (long e = 0; e < nEvt; ++e)
  {
    EventRec &r = ev[e];
    if (r.nInfoRows != 1)
    {
      R.warn["ntp_info rows != 1 for an event"]++;
      R.errors++;
    }
    if (r.nGpointRows != 1)
    {
      R.warn["ntp_gpoint rows != 1 for an event"]++;
      R.errors++;
    }
    if (r.ntrkInfo != static_cast<int>(r.tRows.size()))
    {
      R.warn["ntp_info.ntrk != number of ntp_track rows"]++;
      R.errors++;
    }
    const bool hasVtx = !r.vtx.empty();
    const VtxRow *bv = hasVtx ? &r.vtx[r.best] : nullptr;

    // -- truth particles of this event: id -> gtrack row ---------------------
    std::unordered_map<int, int> rowOfId;
    rowOfId.reserve(r.gRows.size() * 2);
    for (int i : r.gRows)
    {
      const int id = ToInt(g(i, gId), 0);
      if (!rowOfId.emplace(id, i).second)
      {
        R.warn["duplicated (event, gtrackID) row in ntp_gtrack"]++;
        R.errors++;
      }
    }

    // -- reco -> truth multiplicities -----------------------------------------
    // tracks per dominant particle, ranked by nfromtruth (desc), then trackID
    std::unordered_map<int, std::vector<int>> tracksOfId;
    for (int i : r.tRows)
    {
      const int id = ToInt(t(i, tGid), 0);
      if (id != 0) tracksOfId[id].push_back(i);
    }
    std::unordered_map<int, int> rankOfRow;
    for (auto &kv : tracksOfId)
    {
      auto &v = kv.second;
      std::stable_sort(v.begin(), v.end(), [&](int a, int b) {
        const float na = Finite(t(a, tNft)) ? t(a, tNft) : -1.f, nb = Finite(t(b, tNft)) ? t(b, tNft) : -1.f;
        if (na != nb) return na > nb;
        return t(a, tTid) < t(b, tTid);
      });
      for (size_t k = 0; k < v.size(); ++k) rankOfRow[v[k]] = static_cast<int>(k) + 1;
    }

    // -- truth retention -------------------------------------------------------
    std::unordered_map<int, short> keepOf;  // gtrackID -> mask
    std::set<int> v0Ids;
    E.NprimAll = E.NchPrimAll = E.NchPrim = E.Npip = E.Npim = E.NKp = E.NKm = E.Np = E.Npbar = 0;
    for (int i : r.gRows)
    {
      const int id = ToInt(g(i, gId), 0);
      const long long pdg = std::llround(g(i, gFl));
      const PdgInfo &pi = pdgdb.Get(pdg);
      if (!pi.known) R.warn["unknown PDG code " + std::to_string(pdg) + " (charge 0, mass NaN)"]++;
      const bool prim = ToInt(g(i, gPr), 0) == 1;
      const double y = Rapidity(g(i, gPx), g(i, gPy), g(i, gPz), pi.mass);
      const bool broad = Finite(g(i, gPt)) && g(i, gPt) >= cfg.truth_pt_min &&
                         ((Finite(g(i, gEta)) && std::fabs(g(i, gEta)) < cfg.truth_abseta_max) ||
                          (std::isfinite(y) && std::fabs(y) < cfg.truth_absy_max));
      short m = 0;
      if (prim)
      {
        E.NprimAll++;
        if (pi.charge != 0)
        {
          E.NchPrimAll++;
          if (broad) m |= 1;
          if (Finite(g(i, gPt)) && Finite(g(i, gEta)) && g(i, gPt) > cfg.true_pt_min &&
              std::fabs(g(i, gEta)) < cfg.true_abseta_max)
          {
            E.NchPrim++;
            switch (pdg)
            {
            case 211: E.Npip++; break;
            case -211: E.Npim++; break;
            case 321: E.NKp++; break;
            case -321: E.NKm++; break;
            case 2212: E.Np++; break;
            case -2212: E.Npbar++; break;
            default: break;
            }
          }
        }
        if (cfg.keep_v0 && broad && (pdg == 310 || pdg == 3122 || pdg == -3122))
        {
          m |= 2;
          v0Ids.insert(id);
        }
      }
      if (cfg.keep_matched && tracksOfId.count(id)) m |= 8;
      if (m) keepOf[id] = m;
    }
    if (cfg.keep_v0)
      for (int i : r.gRows)
        if (ToInt(g(i, gPr), 0) != 1 && v0Ids.count(ToInt(g(i, gPar), 0))) keepOf[ToInt(g(i, gId), 0)] |= 4;

    // -- write TruthTree rows (gtrack order) -----------------------------------
    for (int i : r.gRows)
    {
      const int id = ToInt(g(i, gId), 0);
      auto kit = keepOf.find(id);
      if (kit == keepOf.end()) continue;
      const long long pdg = std::llround(g(i, gFl));
      const PdgInfo &pi = pdgdb.Get(pdg);
      P.job = R.job;
      P.event = static_cast<int>(e);
      P.gtrackID = id;
      P.pdg = static_cast<int>(pdg);
      P.q = static_cast<short>(pi.charge);
      P.pdgKnown = pi.known;
      P.isPrimary = ToShort(g(i, gPr), -1);
      P.keep = kit->second;
      P.parentId = ToInt(g(i, gPar), 0);
      P.parentPdg = ToInt(g(i, gParFl), 0);
      P.primaryId = ToInt(g(i, gPri), 0);
      P.primaryPdg = ToInt(g(i, gPriFl), 0);
      P.vx = g(i, gVx);
      P.vy = g(i, gVy);
      P.vz = g(i, gVz);
      P.vt = g(i, gVt);
      P.px = g(i, gPx);
      P.py = g(i, gPy);
      P.pz = g(i, gPz);
      P.pt = g(i, gPt);
      P.eta = g(i, gEta);
      P.phi = g(i, gPhi);
      P.mass = static_cast<float>(pi.mass);
      P.y = static_cast<float>(Rapidity(P.px, P.py, P.pz, pi.mass));
      P.gnhits = ToShort(g(i, gNhits));
      P.gnmaps = ToShort(g(i, gNmaps));
      P.gnintt = ToShort(g(i, gNintt));
      P.gntpc = ToShort(g(i, gNtpc));
      P.gnmms = ToShort(g(i, gNmms));
      P.gnlmaps = ToShort(g(i, gNlmaps));
      P.gnlintt = ToShort(g(i, gNlintt));
      P.gnltpc = ToShort(g(i, gNltpc));
      P.gnlmms = ToShort(g(i, gNlmms));
      P.bestTrackID = ToInt(g(i, gTid), -1);
      auto tit = tracksOfId.find(id);
      P.nRecoMatches = tit == tracksOfId.end() ? 0 : static_cast<short>(tit->second.size());
      tTr->Fill();
      R.nTruthKept++;
      for (int b = 0; b < 4; ++b)
        if (P.keep & (1 << b)) R.keepBits[b]++;
    }

    // -- write RecoTree rows (ntp_track order) ---------------------------------
    E.nTrkGood = E.nGoodMVTX1 = E.nGoodMVTX2 = E.nGoodMVTX3 = E.nTrkMatched = 0;
    for (int i : r.tRows)
    {
      K.job = R.job;
      K.event = static_cast<int>(e);
      K.trackID = ToInt(t(i, tTid), -1);
      K.px = t(i, tPx);
      K.py = t(i, tPy);
      K.pz = t(i, tPz);
      K.pt = t(i, tPt);
      K.eta = t(i, tEta);
      K.phi = t(i, tPhi);
      K.charge = ToShort(t(i, tQ), 0);
      K.quality = t(i, tQual);
      K.chisq = t(i, tChi);
      K.ndf = t(i, tNdf);
      K.crossing = t(i, tCross);
      K.nhits = ToShort(t(i, tNhits), 0);
      K.nmaps = ToShort(t(i, tNmaps), 0);
      K.nintt = ToShort(t(i, tNintt), 0);
      K.ntpc = ToShort(t(i, tNtpc), 0);
      K.nmms = ToShort(t(i, tNmms), 0);
      K.nlmaps = ToShort(t(i, tNlmaps), 0);
      K.nlintt = ToShort(t(i, tNlintt), 0);
      K.nltpc = ToShort(t(i, tNltpc), 0);
      K.nlmms = ToShort(t(i, tNlmms), 0);
      K.layers = ToShort(t(i, tLayers), 0);
      K.vertexID = ToInt(t(i, tVid), -1);
      K.evalVtxX = t(i, tVx);
      K.evalVtxY = t(i, tVy);
      K.evalVtxZ = t(i, tVz);
      K.dcaRef = 0;
      for (const auto &v : r.vtx)
        if (std::fabs(K.evalVtxX - v.x) < 1e-4f && std::fabs(K.evalVtxY - v.y) < 1e-4f &&
            std::fabs(K.evalVtxZ - v.z) < 1e-4f)
          K.dcaRef = 1;
      if (K.dcaRef == 0 && std::fabs(K.evalVtxX - r.xT) < 1e-4f && std::fabs(K.evalVtxY - r.yT) < 1e-4f &&
          std::fabs(K.evalVtxZ - r.zT) < 1e-4f)
        K.dcaRef = 2;
      if (K.dcaRef == 0) R.warn["track DCA reference vertex is neither a reco nor the true vertex"]++;
      K.dca2d = t(i, tD2);
      K.dca2dsigma = t(i, tD2s);
      K.dca3dxy = t(i, tDxy);
      K.dca3dxysigma = t(i, tDxys);
      K.dca3dz = t(i, tDz);
      K.dca3dzsigma = t(i, tDzs);
      K.pcax = t(i, tPcax);
      K.pcay = t(i, tPcay);
      K.pcaz = t(i, tPcaz);
      K.xTrue = r.xT;
      K.yTrue = r.yT;
      K.zTrue = r.zT;
      K.xReco = bv ? bv->x : NAN;
      K.yReco = bv ? bv->y : NAN;
      K.zReco = bv ? bv->z : NAN;
      LinDca(K.px, K.py, K.pz, K.pcax, K.pcay, K.pcaz, K.xTrue, K.yTrue, K.zTrue, K.linXYTrue, K.linZTrue);
      LinDca(K.px, K.py, K.pz, K.pcax, K.pcay, K.pcaz, K.xReco, K.yReco, K.zReco, K.linXYReco, K.linZReco);

      // truth match
      const int id = ToInt(t(i, tGid), 0);
      K.gtrackID = id;
      K.matched = id != 0;
      K.singlematch = ToShort(t(i, tSm), -1);
      K.isBestOfTruth = 0;
      K.nSameTruth = 0;
      K.rankSameTruth = 0;
      K.chargeOK = -1;
      K.truthKept = 0;
      K.gpdg = 0;
      K.gq = 0;
      K.gprimary = -1;
      K.gy = NAN;
      if (K.matched)
      {
        const long long pdg = std::llround(t(i, tGfl));
        const PdgInfo &pi = pdgdb.Get(pdg);
        K.gpdg = static_cast<int>(pdg);
        K.gq = static_cast<short>(pi.charge);
        K.chargeOK = (pi.known && pi.charge == K.charge) ? 1 : 0;
        K.gprimary = ToShort(t(i, tGpr), -1);
        K.gy = static_cast<float>(Rapidity(t(i, tGpx), t(i, tGpy), t(i, tGpz), pi.mass));
        const auto &same = tracksOfId[id];
        K.nSameTruth = static_cast<short>(same.size());
        K.rankSameTruth = static_cast<short>(rankOfRow[i]);
        auto rit = rowOfId.find(id);
        if (rit == rowOfId.end())
        {
          R.warn["matched G4 particle not found in ntp_gtrack"]++;
          R.errors++;
        }
        else
        {
          K.isBestOfTruth = ToInt(g(rit->second, gTid), -1) == K.trackID;
        }
        K.truthKept = keepOf.count(id) ? 1 : 0;
        E.nTrkMatched++;
      }
      K.gparentId = ToInt(t(i, tGpar), 0);
      K.gparentPdg = ToInt(t(i, tGparFl), 0);
      K.gprimaryId = ToInt(t(i, tGpri), 0);
      K.gprimaryPdg = ToInt(t(i, tGpriFl), 0);
      K.gpx = t(i, tGpx);
      K.gpy = t(i, tGpy);
      K.gpz = t(i, tGpz);
      K.gpt = t(i, tGpt);
      K.geta = t(i, tGeta);
      K.gphi = t(i, tGphi);
      K.gvx = t(i, tGvx);
      K.gvy = t(i, tGvy);
      K.gvz = t(i, tGvz);
      K.nfromtruth = ToShort(t(i, tNft));
      K.nwrong = ToShort(t(i, tNw));
      K.ntrumaps = ToShort(t(i, tNtm));
      K.nwrongmaps = ToShort(t(i, tNwm));
      K.ntruintt = ToShort(t(i, tNti));
      K.nwrongintt = ToShort(t(i, tNwi));
      K.ntrutpc = ToShort(t(i, tNtt));
      K.nwrongtpc = ToShort(t(i, tNwt));
      K.ntrumms = ToShort(t(i, tNtmm));
      K.nwrongmms = ToShort(t(i, tNwmm));
      K.layersfromtruth = ToShort(t(i, tLft));
      K.purity = (K.nhits > 0 && K.nfromtruth >= 0) ? static_cast<float>(K.nfromtruth) / K.nhits : NAN;
      tR->Fill();

      // event-level summary multiplicities (same definition as VtxEffAnalysis n)
      if (Finite(K.pt) && Finite(K.eta) && K.pt > cfg.good_pt_min && std::fabs(K.eta) < cfg.good_abseta_max &&
          K.ntpc >= cfg.good_ntpc_min)
      {
        E.nTrkGood++;
        if (K.nmaps >= 1) E.nGoodMVTX1++;
        if (K.nmaps >= 2) E.nGoodMVTX2++;
        if (K.nmaps >= 3) E.nGoodMVTX3++;
      }
    }

    // -- EventTree row -----------------------------------------------------------
    E.job = R.job;
    E.event = static_cast<int>(e);
    E.seed = r.seed;
    E.xTrue = r.xT;
    E.yTrue = r.yT;
    E.zTrue = r.zT;
    E.tTrue = r.tT;
    E.hasVtx = hasVtx;
    E.nVtx = static_cast<short>(r.vtx.size());
    E.nVtxRowsNaN = static_cast<short>(r.nVtxRowsNaN);
    E.vtxId = bv ? bv->id : -1;
    E.xReco = bv ? bv->x : NAN;
    E.yReco = bv ? bv->y : NAN;
    E.zReco = bv ? bv->z : NAN;
    E.vtxChi2 = bv ? bv->chi2 : NAN;
    E.vtxNdf = bv ? bv->ndof : NAN;
    E.vtxCrossing = bv ? bv->crossing : NAN;
    E.vtxNtracks = bv ? ToShort(bv->ntracks, -1) : -1;
    E.ntrkInfo = static_cast<short>(r.ntrkInfo);
    E.gntrkallInfo = static_cast<short>(r.gntrkallInfo);
    E.nTrkAll = static_cast<short>(r.tRows.size());
    E.nTruthKept = static_cast<short>(keepOf.size());
    tE->Fill();
    R.nEvents++;
    R.nVertexed += hasVtx;
  }

  // provenance
  out.cd();
  TNamed("input_file", path.c_str()).Write();
  TNamed("program", kVersion).Write();
  TNamed("config", cfg.dump.c_str()).Write();
  TNamed("retention",
         Form("keep bits: 1 = charged G4 primary, gpt >= %g and (|eta| < %g or |y| < %g); "
              "2 = G4 primary K0S/Lambda/anti-Lambda in the same acceptance (keep_v0=%d); "
              "4 = direct daughter of a bit-2 particle (keep_v0=%d); "
              "8 = dominant truth particle of >= 1 reco track (keep_matched=%d). All other G4 particles discarded.",
              cfg.truth_pt_min, cfg.truth_abseta_max, cfg.truth_absy_max, cfg.keep_v0, cfg.keep_v0, cfg.keep_matched))
      .Write();
  TNamed("dca_note",
         "dca2d/dca3dxy/dca3dz/pca: original evaluator values w.r.t. evalVtx (dcaRef 1 = reco vertex, 2 = TRUE vertex). "
         "dca3dz = pcaz - evalVtxZ without extrapolation; for dcaRef 2 the track state (pca) is the transverse closest "
         "approach to the ORIGIN (0,0), so dca3dz is not the z distance at the closest approach to the true vertex. "
         "lin*: straight line through pca along the reco momentum, evaluated at the transverse closest approach to the "
         "true/reco vertex - an approximation, NOT the helix DCA.")
      .Write();
  tE->Write();
  tTr->Write();
  tR->Write();
  out.Close();

  sw.Stop();
  R.seconds = sw.RealTime();
  return R;
}

void Usage()
{
  std::cout << "Usage: EvalSkim -c skim.cfg [-o OUTDIR] [-l file.list] [-n maxEventsPerFile] [file1.root ...]\n"
               "  -c  configuration file (required)\n"
               "  -o  output directory (overrides output_dir); writes OUTDIR/skim_NNNNNN.root per input job\n"
               "  -l  overrides input_list from the config\n"
               "  -n  process only the first N events of each file (tests; 0 = all)\n"
               "  extra arguments are input files (override the list)\n";
}
}  // namespace

int main(int argc, char **argv)
{
  std::string cfgPath, listOverride, outOverride;
  long maxEv = -1;
  std::vector<std::string> files;
  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "-c" && i + 1 < argc) cfgPath = argv[++i];
    else if (a == "-l" && i + 1 < argc) listOverride = argv[++i];
    else if (a == "-o" && i + 1 < argc) outOverride = argv[++i];
    else if (a == "-n" && i + 1 < argc) maxEv = std::stol(argv[++i]);
    else if (a == "-h" || a == "--help")
    {
      Usage();
      return 0;
    }
    else files.push_back(a);
  }
  if (cfgPath.empty())
  {
    Usage();
    return 1;
  }
  SkimConfig cfg = ReadSkimConfig(cfgPath);
  if (!listOverride.empty()) cfg.input_list = listOverride;
  if (!outOverride.empty()) cfg.output_dir = outOverride;
  if (maxEv >= 0) cfg.max_events_per_file = maxEv;
  if (files.empty()) files = ReadFileList(cfg.input_list, cfg.input_dir);
  if (files.empty() || cfg.output_dir.empty())
  {
    std::cerr << "FATAL: no input files or no output directory" << std::endl;
    return 1;
  }
  if (!cfg.input_dir.empty() && cfg.output_dir.rfind(cfg.input_dir, 0) == 0)
  {
    std::cerr << "FATAL: output directory is inside input_dir - refusing" << std::endl;
    return 1;
  }
  gSystem->mkdir(cfg.output_dir.c_str(), true);

  std::cout << "=== " << kVersion << " ===\n"
            << "config  : " << cfgPath << "\n"
            << "output  : " << cfg.output_dir << "\n"
            << "files   : " << files.size() << "\n"
            << "truth   : charged G4 primaries gpt >= " << cfg.truth_pt_min << ", |eta| < " << cfg.truth_abseta_max
            << " or |y| < " << cfg.truth_absy_max << "; keep_v0 = " << cfg.keep_v0
            << "; keep_matched = " << cfg.keep_matched << "\n"
            << "events/file : " << (cfg.max_events_per_file > 0 ? std::to_string(cfg.max_events_per_file) : std::string("all"))
            << "\n";

  PdgLookup pdgdb;
  std::set<int> jobsSeen;
  int totalErrors = 0;
  long totEv = 0, totVtx = 0, totTruth = 0, totKept = 0, totTrk = 0;
  double totSec = 0;
  std::map<std::string, long> warn;
  for (const auto &fpath : files)
  {
    const int job = JobFromName(fpath);
    if (job < 0)
    {
      std::cerr << "ERROR: cannot parse job number from " << fpath << " - skipped" << std::endl;
      totalErrors++;
      continue;
    }
    if (!jobsSeen.insert(job).second)
    {
      std::cerr << "ERROR: job " << job << " appears twice in the input - skipped" << std::endl;
      totalErrors++;
      continue;
    }
    const std::string outPath = cfg.output_dir + Form("/skim_%06d.root", job);
    FileResult R = ProcessFile(fpath, outPath, cfg, pdgdb);
    totalErrors += R.errors;
    totEv += R.nEvents;
    totVtx += R.nVertexed;
    totTruth += R.nTruthRows;
    totKept += R.nTruthKept;
    totTrk += R.nTracks;
    totSec += R.seconds;
    for (const auto &w : R.warn) warn[w.first] += w.second;
    std::cout << Form("[file] job=%06d events=%ld vertexed=%ld gtrackRows=%ld truthKept=%ld (bit1=%ld bit2=%ld bit4=%ld "
                      "bit8=%ld) tracks=%ld errors=%d time=%.2fs -> %s",
                      R.job, R.nEvents, R.nVertexed, R.nTruthRows, R.nTruthKept, R.keepBits[0], R.keepBits[1],
                      R.keepBits[2], R.keepBits[3], R.nTracks, R.errors, R.seconds, outPath.c_str())
              << std::endl;
  }
  std::cout << Form("[total] files=%zu events=%ld vertexed=%ld gtrackRows=%ld truthKept=%ld tracks=%ld time=%.1fs",
                    jobsSeen.size(), totEv, totVtx, totTruth, totKept, totTrk, totSec)
            << std::endl;
  std::cout << "--- warnings / diagnostics ---" << std::endl;
  if (warn.empty()) std::cout << "(none)" << std::endl;
  for (const auto &w : warn) std::cout << "  " << w.second << " x " << w.first << std::endl;
  std::cout << "structural errors: " << totalErrors << std::endl;
  return totalErrors == 0 ? 0 : 3;
}
