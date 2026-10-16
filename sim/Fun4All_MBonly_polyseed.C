// =============================================================================
// Fun4All_MBonly_polyseed.C
//
// Pileup-free p+p 200 GeV PYTHIA8 minimum-bias (Detroit tune), full GEANT4
// simulation of the tracking detectors, polyseed TPC tracking, ONE tracking
// vertex finder, SvtxEvaluator ntuples and a DST.
//
// Based on M. Peters' Fun4All_HFG_MBonly_polyseed.C (version of 2026-10-12)
// and on the official polyseed chain in macros/TrackingProduction/
// Fun4All_PRDF_SAseed_run3pp.C (master after PR #1393). Every difference from
// Peters' macro is listed in CHANGES.md and marked "CHANGE" below.
//
// Target release: first release containing coresoftware PR #4473 (merge
// 3e4855d97) and macros PR #1393 - the rebuilt new.1 of 2026-10-14
// (coresoftware bb9d98e). Detected by API, not by release name (names cycle). The MC timing fix (TPC ADC clock from the simulated geometry,
// t0 = 1.5 time bins, 387.5 V/cm) is applied by the shared PHGarfield with
// SetIsNonDistortedMC(true) - nothing is set by hand on the polyseed modules.
//
// Usage:
//   root -l -b -q 'Fun4All_MBonly_polyseed.C(nEvents, jobId, "outDir"[, manualClockFix, randomSeed])'
//   randomSeed = false (default): fixed seed 100000 + jobId (reproducible; K1, 10k samples)
//   randomSeed = true : no RANDOMSEED flag, as in Peters' production - PHRandomSeed draws
//                       every seed from std::random_device and prints it to the log
//                       ("PHRandomSeed::GetSeed() seed: N"); ntp_info.seed is then 0
// Output (in outDir):
//   DST/MBonly_polyseed_DST_<jobId>.root
//   evaluator/outputeval_<jobId>.root
// =============================================================================
#ifndef MACRO_FUN4ALL_MBONLY_POLYSEED_C
#define MACRO_FUN4ALL_MBONLY_POLYSEED_C

#include <G4_Input.C>  // CHANGE: release copy (Peters' local copy is identical to new.20)
#include <G4_Global.C>
#include <G4Setup_sPHENIX.C>
#include <SaveGitTags.C>

#include <Trkr_Clustering.C>
#include <Trkr_Eval.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <phgarfield/PHGarfield.h>
#include <tpctrackreco/TpcCrossingFinder.h>
#include <tpctrackreco/TpcPolyClusterTrkrClusterConverter.h>
#include <tpctrackreco/TpcPolyTrackSeedConverter.h>
#include <tpctrackreco/Tpc_AssembledTrackReco.h>
#include <tpctrackreco/Tpc_ModuleTrackReco.h>
#include <tpctrackreco/Tpc_PolyClusterizer.h>
#include <tpctrackreco/Tpc_PolyTrackReco.h>
#include <tpctrackreco/Tpc_PolyTrackVertexer.h>

#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libPHGarfield.so)
R__LOAD_LIBRARY(libtpctrackreco.so)

namespace MBPOLY
{
  // --- generator -------------------------------------------------------------
  // Detroit tune, copied unchanged from the card Peters used:
  // /sphenix/user/cdean/.../steeringCards/pythia8_MB_Detroit_Tony.cfg
  std::string pythia_card = "pythia8_MB_Detroit_Tony.cfg";  // resolved next to this macro

  // --- conditions --------------------------------------------------------------
  // CHANGE (2026-10-14): official MC tag "MDC2" (Peters: data tag "newcdbtag").
  // Experts: do not use the data cdbtag for simulations. newcdbtag's TRACKINGALIGNMENT
  // is the 2025 DATA alignment and displaced the reconstructed silicon in MC by
  // 0.3-0.6 mm (fix_comparison/audit_mvtx/); MDC2's is all zero, matching the
  // simulated geometry. Test mdc2_test/ (K1, same seeds): eps_vtx(|z|<8) 0.654 ->
  // 0.735, identical to newcdbtag + zero alignment. MDC2 has no polyseed payloads:
  // the field maps are given explicitly below; the gas falls back to
  // /cvmfs/sphenix.sdcc.bnl.gov/files/gasfiles (same v_drift 7.551 cm/us at 387.5 V/cm).
  const std::string cdb_tag = "MDC2";
  int run_number = 29;  // MC run number (sets the simulated TPC ADC clock 53.326 ns)

  // The four polyseed field maps, pinned to the files M. Peters' production
  // loaded (they are what newcdbtag returns; they do not exist under MDC2).
  // In non-distorted MC their contributions are scaled to zero (space charge 0,
  // frame charge 0), so only their presence is needed.
  const std::string cdbdir = "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/";
  const std::string efield3d_side0 = cdbdir + "Tpc_PolySeeding_EField3D_Side0/ba/76/ba7618a6247314a270a9ee4b95530f39_phgarfield_ibf_3d_map_gainOn_alpha_1p8_flatphi_side0.root";
  const std::string efield3d_side1 = cdbdir + "Tpc_PolySeeding_EField3D_Side1/6b/fd/6bfd3df3e1ef1dc50931e28171240431_phgarfield_ibf_3d_map_gainOn_alpha_1p8_flatphi_side1.root";
  const std::string frames3d_side0 = cdbdir + "Tpc_PolySeeding_ModuleFrames3D_Side0/05/0c/050c53906b9f27aa2cb69733a3b863df_phgarfield_tpc_moduleframes_3d_side0.root";
  const std::string frames3d_side1 = cdbdir + "Tpc_PolySeeding_ModuleFrames3D_Side1/8d/c2/8dc290303f3c6caabf4a48fea0c52925_phgarfield_tpc_moduleframes_3d_side1.root";

  // --- detector options ----------------------------------------------------------
  bool mvtx_misalignment = true;  // as Peters and the release default (open decision)

  // --- release feature detection ---------------------------------------------------
  // Release names cycle (new.1 ... new.20), so the polyseed fix is detected by
  // its API, not by name: PR #4473 added PHGarfield::GetNonDistortedMCT0Tbins().
  // true  -> shared-PHGarfield design (fix built in)
  // false -> legacy release (e.g. new.20 of 2026-10-13, coresoftware 2cb2898)
  bool releaseHasPolyseedFix()
  {
    const char *offline = gSystem->Getenv("OFFLINE_MAIN");
    if (!offline) return false;
    std::ifstream header(std::string(offline) + "/include/phgarfield/PHGarfield.h");
    std::string line;
    while (std::getline(header, line))
    {
      if (line.find("GetNonDistortedMCT0Tbins") != std::string::npos) return true;
    }
    return false;
  }
}  // namespace MBPOLY

// Silicon seeding and the final tracking vertex use the RELEASE functions
// Tracking_Reco_SiliconSeed_run2pp() and Tracking_Reco_Vertex_run2pp()
// (Trkr_Reco.C), unchanged. With the MC tag MDC2, TrackingInit() sets
// CDB::is_data_reco = false, so these give the official MC reconstruction:
// silicon-seeding beam spot (0,0) (Trkr_Reco.C:214,231) and the
// PHSimpleVertexFinder default beam-spot box |x|,|y| < 0.2 cm around (0,0).
// The simulated collisions are at the official p+p beam spot (-0.07, 0.141) cm
// (Input::ApplysPHENIXBeamParameter, G4_Input.C). No beam numbers are set in
// this macro. (Until 2026-10-15 the macro set the simulated beam spot in the
// seeding and centred the vertex box on it; test fix_comparison/beamspot_release_test/.)
// CHANGE: Peters' production registered a SECOND, loose PHSimpleVertexFinder
// (dca 1, no MVTX requirement) after the release one; it is not used here.

int Fun4All_MBonly_polyseed(const int nEvents = 10,
                            const int jobId = 0,
                            const std::string &outDir = "./",
                            // legacy releases only: apply the TPC timing fix by hand (ADC clock
                            // 53.326 ns of the simulation, t0 = 1.5 time bins, drift
                            // start +-102.605 cm). Fixed releases apply it themselves (ignored there).
                            const bool manualClockFix = false,
                            // true: random seeds as in Peters' production (see Usage)
                            const bool randomSeed = false)
{
  // ---------------------------------------------------------------------------
  // 0. Bookkeeping
  // ---------------------------------------------------------------------------
  std::ostringstream jobStr;
  jobStr << std::setw(6) << std::setfill('0') << jobId;
  const std::string tag = jobStr.str();
  gSystem->mkdir((outDir + "/DST").c_str(), true);
  gSystem->mkdir((outDir + "/evaluator").c_str(), true);
  const std::string macroDir = gSystem->DirName(__FILE__);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();

  // Seeds. Default: fixed, reproducible seed per job (CHANGE w.r.t. Peters). PHRandomSeed
  // (and therefore PYTHIA8, GEANT4, digitisation) uses this flag.
  // randomSeed = true: as Peters - no flag, every seed from std::random_device, printed.
  const int seed = randomSeed ? 0 : 100000 + jobId;  // 0 = random (no RANDOMSEED flag)
  if (!randomSeed)
  {
    rc->set_IntFlag("RANDOMSEED", seed);
    std::cout << "Fun4All_MBonly_polyseed: fixed seed " << seed << " (100000 + jobId)" << std::endl;
  }
  else
    std::cout << "Fun4All_MBonly_polyseed: RANDOM seeds (no RANDOMSEED flag, as Peters); every seed is printed by PHRandomSeed" << std::endl;
  PHRandomSeed::Verbosity(1);

  Enable::VERBOSITY = 0;
  Input::VERBOSITY = 0;

  // ---------------------------------------------------------------------------
  // 1. Generator: PYTHIA8 minimum bias, Detroit tune, p+p 200 GeV
  // ---------------------------------------------------------------------------
  Input::READHITS = false;
  Input::SIMPLE = false;  // no embedded particles (MB only)
  Input::PYTHIA8 = true;
  PYTHIA8::config_file[0] = macroDir + "/" + MBPOLY::pythia_card;
  Input::BEAM_CONFIGURATION = Input::pp_COLLISION;

  InputInit();
  // sPHENIX p+p beam: vertex mean (-0.07, 0.141, 0) cm, sigma_z 16 cm,
  // crossing angle 1.5 mrad (G4_Input.C:146-158)
  Input::ApplysPHENIXBeamParameter(INPUTGENERATOR::Pythia8[0]);
  InputRegister();

  // ---------------------------------------------------------------------------
  // 2. Conditions database and run flags
  // ---------------------------------------------------------------------------
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", MBPOLY::cdb_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);
  rc->set_IntFlag("RUNNUMBER", MBPOLY::run_number);

  Enable::MVTX_APPLYMISALIGNMENT = MBPOLY::mvtx_misalignment;  // same value in G4 and Acts
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  // CHANGE: no SAMPA t0 bias. It is applied only in reconstruction
  // (ActsGeometry, poly-cluster converter); with -65 ns the converter leaves
  // 3 x bias = -195 ns (~1.5 cm) on TPC side 0 only. 0 is consistent on both sides.
  G4TPC::sampa_tzero_bias = 0.;

  auto *flag = new FlagHandler();
  se->registerSubsystem(flag);

  // ---------------------------------------------------------------------------
  // 3. Detector simulation
  // ---------------------------------------------------------------------------
  Enable::MBDFAKE = true;
  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;
  TRACKING::streaming_mode = true;  // pp mode of the vertex finder (as Peters)

  InputManagers();
  G4Init();
  MagnetInit();
  MagnetFieldInit();
  G4Setup();

  Mbd_Reco();
  Mvtx_Cells();
  Intt_Cells();
  TPC_Cells();
  Micromegas_Cells();

  // ---------------------------------------------------------------------------
  // 4. Tracking: geometry, silicon/TPOT clustering
  // ---------------------------------------------------------------------------
  TrackingInit();

  Mvtx_Clustering();
  Intt_Clustering();
  // CHANGE: no TPC_Clustering(). TpcPolyClusterTrkrClusterConverter removes all
  // standard TPC clusters before writing the poly clusters
  // (clearOutputTpcClusters), as in the official polyseed chain.
  Micromegas_Clustering();

  // ---------------------------------------------------------------------------
  // 5. Silicon seeds and silicon-only vertex (input to the crossing finder)
  // ---------------------------------------------------------------------------
  Tracking_Reco_SiliconSeed_run2pp();  // release function (beam spot: release default)

  auto *converter = new TrackSeedTrackMapConverter("SiliconSeedToSvtxTrackMap");
  converter->setTrackSeedName("SiliconTrackSeedContainer");
  converter->setTrackMapName("SiliconSvtxTrackMap");
  converter->setClusterMapName("TRKR_CLUSTER");
  se->registerSubsystem(converter);

  auto *finder_svx = new PHSimpleVertexFinder("SiliconVertexFinder");
  finder_svx->Verbosity(0);
  finder_svx->setDcaCut(0.1);
  finder_svx->setTrackPtCut(0.2);
  finder_svx->setBeamLineCut(1);
  finder_svx->setTrackQualityCut(500);
  finder_svx->setNmvtxRequired(3);
  finder_svx->setOutlierPairCut(0.1);
  finder_svx->setTrackMapName("SiliconSvtxTrackMap");
  finder_svx->setVertexMapName("SiliconSvtxVertexMap");
  se->registerSubsystem(finder_svx);

  // ---------------------------------------------------------------------------
  // 6. Polyseed TPC tracking
  // ---------------------------------------------------------------------------
  // CHANGE (PR #4473 + PR #1393): one shared PHGarfield, registered before the
  // crossing finder and clusterizer. SetIsNonDistortedMC(true) gives: zero
  // space charge and field-cage distortions, 387.5 V/cm, t0 = 1.5 time bins and
  // the TPC ADC clock from the simulated geometry (53.326 ns). The modules read
  // these in InitRun and print "non-distorted MC: t0=... ADC clock=...".
  // CHANGE: no setT0(0), setCMVoltageDefault, setIsNonDistortedMC or
  // setUseSurveyGeometry on the modules (Peters set them; with the shared
  // PHGarfield they are ignored, and setT0(0) would override t0 = 1.5).
  auto *garfield = new PHGarfield();
  garfield->SetUseSurveyGeometry(true);  // simulation places the TPC at the survey position (G4TPC::tpc_survey_position = true)
  garfield->SetIsNonDistortedMC(true);
  garfield->SetElectricFieldMap3D(MBPOLY::efield3d_side0, MBPOLY::efield3d_side1);  // see MBPOLY (not in MDC2)
  garfield->SetFrameElectricFieldMap3D(MBPOLY::frames3d_side0, MBPOLY::frames3d_side1);
  se->registerSubsystem(garfield);

  se->registerSubsystem(new Tpc_ModuleTrackReco());     // makes TPC_MODULETRACKS
  se->registerSubsystem(new Tpc_AssembledTrackReco());  // makes TPC_ASSEMBLEDTRACKS

  // Releases before PR #4473 (e.g. new.20, coresoftware 2cb2898) do not have the shared-PHGarfield
  // design: the clusterizer builds a PRIVATE PHGarfield from its own setters,
  // with distortions ON by default. For those releases only, the polyseed
  // modules are configured exactly as in Peters' macro (setT0(0), 387.5 V/cm,
  // non-distorted MC, no survey geometry) plus the field maps. The ADC clock
  // stays at the old default 56.881 ns (PRE-FIX behaviour) unless
  // manualClockFix is set, which applies the values of the PR #4473 fix by hand.
  const bool legacyPolyseed = !MBPOLY::releaseHasPolyseedFix();
  const bool applyClockFix = legacyPolyseed && manualClockFix;
  const double mcAdcClock = 53.326184;  // ns, simulated TPC ADC clock (PHG4TpcDetector "TBinWidth 53.3262")
  const double mcT0 = 1.5;              // time bins, PR #4473 non-distorted MC value
  std::cout << "Fun4All_MBonly_polyseed: polyseed configuration = "
            << (!legacyPolyseed ? "shared PHGarfield (PR #4473 release; clock fix built in)"
                                : (applyClockFix ? "LEGACY (pre-PR #4473) + MANUAL CLOCK FIX (53.326 ns, t0 1.5)"
                                                 : "LEGACY (pre-PR #4473, as Peters; pre-fix clock 56.881 ns, t0 0)"))
            << std::endl;
  if (!legacyPolyseed && manualClockFix)
    std::cout << "Fun4All_MBonly_polyseed: manualClockFix ignored (not needed for this release)" << std::endl;

  auto *crossingFinder = new TpcCrossingFinder();
  crossingFinder->Verbosity(0);
  crossingFinder->setTriggeredMode(true);  // single collision per event (as Peters)
  crossingFinder->setInputNodeName("TPC_ASSEMBLEDTRACKS");
  crossingFinder->setOutputNodeName("TPC_CROSSING_DECISIONS");
  crossingFinder->setVertexMapNodeName("SiliconSvtxVertexMap");
  if (legacyPolyseed)
  {
    crossingFinder->setCMVoltageDefault(387.5);
    crossingFinder->setT0(applyClockFix ? mcT0 : 0);
    crossingFinder->setIsNonDistortedMC(true);
    if (applyClockFix)
    {
      crossingFinder->setTpcAdcClock(mcAdcClock);  // default 56.881262 = data clock
      crossingFinder->setStartZ(-102.605, 102.605);  // geometry drift start (default +-102.325)
    }
  }
  se->registerSubsystem(crossingFinder);

  auto *cluster = new Tpc_PolyClusterizer();  // makes TPC_POLYCLUSTERS (fixed releases: uses the shared PHGarfield)
  if (legacyPolyseed)
  {
    cluster->setIsNonDistortedMC(true);
    cluster->setCMVoltageDefault(387.5);
    cluster->setT0(applyClockFix ? mcT0 : 0);
    if (applyClockFix) cluster->setTpcAdcClock(mcAdcClock);  // drift start already +-102.605 by default
    cluster->setUseSurveyGeometry(false);
    cluster->setElectricFieldMap3DSide0(MBPOLY::efield3d_side0);
    cluster->setElectricFieldMap3DSide1(MBPOLY::efield3d_side1);
    cluster->setFrameElectricFieldMap3DSide0(MBPOLY::frames3d_side0);
    cluster->setFrameElectricFieldMap3DSide1(MBPOLY::frames3d_side1);
  }
  se->registerSubsystem(cluster);
  se->registerSubsystem(new Tpc_PolyTrackReco());      // makes TPC_POLYTRACKS
  se->registerSubsystem(new Tpc_PolyTrackVertexer());  // makes TPC_POLYTRACKVERTICES

  se->registerSubsystem(new TpcPolyTrackSeedConverter());           // TPC_POLYTRACKS -> TpcTrackSeedContainer
  se->registerSubsystem(new TpcPolyClusterTrkrClusterConverter());  // TPC_POLYCLUSTERS -> TRKR_CLUSTER

  // ---------------------------------------------------------------------------
  // 7. Silicon-TPC matching, track fit, ONE tracking vertex finder
  // ---------------------------------------------------------------------------
  Tracking_Reco_TrackMatching_run2pp();
  // CHANGE: no DSTClusterPruning (not in the official chain; Peters added it).
  Tracking_Reco_TrackFit_run2pp();
  Tracking_Reco_Vertex_run2pp();  // release function (beam-spot box: release default)

  // ---------------------------------------------------------------------------
  // 8. Global vertex, evaluator (unchanged SvtxEvaluator settings = baseline)
  // ---------------------------------------------------------------------------
  Global_Reco();
  // CHANGE: Tracking_Eval() already calls build_truthreco_tables()
  // (Trkr_Eval.C:25); Peters called it a second time.
  Tracking_Eval(outDir + "/evaluator/outputeval_" + tag + ".root");

  // ---------------------------------------------------------------------------
  // 9. DST
  // ---------------------------------------------------------------------------
  // CHANGE: G4 hits, hit sets and truth associations are KEPT (Peters stripped
  // them) so that the evaluator can be re-run on this DST for the fix studies.
  // CHANGE: no KFParticle V0 reconstruction (HF_selections.C) in this macro.
  SaveGitTags();  // release/commit tags into the run node
  const std::string dstFile = outDir + "/DST/MBonly_polyseed_DST_" + tag + ".root";
  auto *out = new Fun4AllDstOutputManager("DSTOUT", dstFile);
  out->StripNode("G4HIT_PIPE");
  out->StripNode("G4HIT_SVTXSUPPORT");
  se->registerOutputManager(out);

  std::cout << "Fun4All_MBonly_polyseed: job " << tag << " seed " << (randomSeed ? std::string("random") : std::to_string(seed)) << " events " << nEvents
            << " CDB " << MBPOLY::cdb_tag << " run " << MBPOLY::run_number
            << " sampa_tzero_bias " << G4TPC::sampa_tzero_bias << " ns"
            << " manualClockFix=" << manualClockFix
            << " OFFLINE_MAIN " << gSystem->Getenv("OFFLINE_MAIN") << std::endl;

  se->run(nEvents);
  se->End();
  std::cout << "Fun4All_MBonly_polyseed: done, DST " << dstFile << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}

#endif
