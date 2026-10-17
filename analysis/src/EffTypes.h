// =============================================================================
// EffTypes.h - the species and pair types of the efficiency module, in one place
//
// Included by the C++ programs (src/) and by the ROOT macros (macros/).
//
// Species: one efficiency map each.
//   kTrack  identified charged hadron (pi, K, p), binned in TRUE rapidity y
//   kHadron every charged primary, PID ignored (Corral h+, h-), binned in TRUE
//           pseudorapidity eta, pion pT cut and binning
//   kV0     K0S, Lambda, Lbar: their maps come from the separate V0 module;
//           this module skips them but keeps them so the pair indices match
//
// Pair types: the SAME index, order and y bin width as Corral src/PairTypes.h
// (49 pair types, 2026-10-11). Pair type i here is pair type i in Corral.
// Adding a pair type: add the same line here and in Corral's PairTypes.h.
//
// Windows follow Corral src/fluct_common.h: Species_yu (pi 1.0, K 0.7, p 0.6,
// V0 0.7), Hadron_etau 1.1, Species_ptmin (p, pbar 0.4), Species_pmax_pid
// (pi 5.0, K 0.6, p 1.0 GeV/c; the dE/dx PID ptot caps, for reference only:
// the PID efficiency is NOT part of these maps).
// =============================================================================
#ifndef EFFTYPES_H
#define EFFTYPES_H

#include <Rtypes.h>  // ROOT colours and markers

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace efftypes
{
enum Kind
{
  kTrack,
  kHadron,
  kV0
};

struct Species
{
  const char *tag;    // histogram names (pip, pim, ..., hp, hm, k0s, lam, lbar)
  const char *label;  // ROOT LaTeX label
  int pdg;            // signed PDG code; 0 for kHadron (any charged primary)
  int charge;         // electric charge in units of e
  Kind kind;
  bool useEta;          // true: binned in eta (kHadron), false: in y
  double absMax;        // window |y| (or |eta|) < absMax
  double ptMin;         // window pT >= ptMin (GeV/c)
  double pidPmax;       // dE/dx PID ptot cap in Corral (GeV/c); 0 = no PID
  double ptAnaMax;      // integrated numbers use pT < ptAnaMax (kaons: QA bin above 0.6)
  std::vector<double> y;   // map edges in y (or eta); empty for kV0
  std::vector<double> pt;  // map edges in pT; empty for kV0
  int color, marker;
  const char *text;  // window description for plots and logs
};

// species index (the order of Species() below)
enum SpeciesIndex
{
  iPip,
  iPim,
  iKp,
  iKm,
  iP,
  iPbar,
  iHp,
  iHm,
  iK0s,
  iLam,
  iLbar,
  NSpecies
};

inline const std::vector<Species> &AllSpecies()
{
  static const std::vector<double> yPi = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  static const std::vector<double> ptPi = {0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60,
                                           0.80, 1.00, 1.25, 1.50, 2.00, 3.00, 5.00};
  static const std::vector<double> yK = {-0.7, -0.4, 0.0, 0.4, 0.7};
  static const std::vector<double> ptK = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0};
  static const std::vector<double> yP = {-0.6, -0.3, 0.0, 0.3, 0.6};
  static const std::vector<double> ptP = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
  static const std::vector<double> etaH = {-1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1};
  static const std::vector<double> none;
  static const std::vector<Species> s = {
      {"pip", "#pi^{+}", 211, +1, kTrack, false, 1.0, 0.1, 5.0, 99, yPi, ptPi, kAzure + 2, 20, "|y| < 1.0, pT >= 0.1 GeV/c"},
      {"pim", "#pi^{-}", -211, -1, kTrack, false, 1.0, 0.1, 5.0, 99, yPi, ptPi, kAzure + 2, 24, "|y| < 1.0, pT >= 0.1 GeV/c"},
      {"kp", "K^{+}", 321, +1, kTrack, false, 0.7, 0.1, 0.6, 0.6, yK, ptK, kRed + 1, 21,
       "|y| < 0.7, 0.1 <= pT < 0.6 GeV/c (0.6-1.0: QA bin)"},
      {"km", "K^{-}", -321, -1, kTrack, false, 0.7, 0.1, 0.6, 0.6, yK, ptK, kRed + 1, 25,
       "|y| < 0.7, 0.1 <= pT < 0.6 GeV/c (0.6-1.0: QA bin)"},
      {"p", "p", 2212, +1, kTrack, false, 0.6, 0.4, 1.0, 99, yP, ptP, kGreen + 2, 22, "|y| < 0.6, 0.4 <= pT < 1.0 GeV/c"},
      {"pbar", "#bar{p}", -2212, -1, kTrack, false, 0.6, 0.4, 1.0, 99, yP, ptP, kGreen + 2, 26,
       "|y| < 0.6, 0.4 <= pT < 1.0 GeV/c"},
      {"hp", "h^{+}", 0, +1, kHadron, true, 1.1, 0.1, 0, 99, etaH, ptPi, kGray + 2, 33, "|#eta| < 1.1, pT >= 0.1 GeV/c"},
      {"hm", "h^{-}", 0, -1, kHadron, true, 1.1, 0.1, 0, 99, etaH, ptPi, kGray + 2, 27, "|#eta| < 1.1, pT >= 0.1 GeV/c"},
      {"k0s", "K^{0}_{S}", 310, 0, kV0, false, 0.7, 0.1, 0, 99, none, none, kMagenta + 1, 29, "V0 module"},
      {"lam", "#Lambda", 3122, 0, kV0, false, 0.7, 0.1, 0, 99, none, none, kOrange + 7, 34, "V0 module"},
      {"lbar", "#bar{#Lambda}", -3122, 0, kV0, false, 0.7, 0.1, 0, 99, none, none, kOrange + 7, 28, "V0 module"}};
  return s;
}

// the species this module measures (kTrack, and kHadron once enabled by the caller)
inline std::vector<const Species *> MeasuredSpecies(bool withHadrons)
{
  std::vector<const Species *> v;
  for (const auto &s : AllSpecies())
    if (s.kind == kTrack || (withHadrons && s.kind == kHadron)) v.push_back(&s);
  return v;
}

inline const Species *FindByTag(const char *tag)
{
  for (const auto &s : AllSpecies())
    if (std::string(s.tag) == tag) return &s;
  return nullptr;
}

// -----------------------------------------------------------------------------
// Pair types: {species 1, species 2, y bin width in units of 0.01} - Corral
// PairTypes_Info[i][0], [1], [3] (field [2], N_min, is Corral's own).
// -----------------------------------------------------------------------------
struct PairType
{
  int s1, s2;  // SpeciesIndex
  int bw;      // y (eta) bin width in units of 0.01, for y1, y2 and dy
};

static const int NPairTypes = 49;
inline const PairType *PairTypes()
{
  static const PairType p[NPairTypes] = {
      {iPip, iPim, 4},    // 0  pi+ pi-
      {iPim, iPip, 4},    // 1  pi- pi+
      {iPip, iPip, 4},    // 2  pi+ pi+
      {iPim, iPim, 4},    // 3  pi- pi-
      {iKp, iKm, 10},     // 4  K+ K- (phi)
      {iKp, iKp, 10},     // 5  K+ K+
      {iKm, iKm, 10},     // 6  K- K-
      {iP, iPbar, 8},     // 7  p pbar
      {iP, iP, 6},        // 8  p p
      {iPbar, iPbar, 8},  // 9  pbar pbar
      {iPip, iKm, 4},     // 10 pi+ K-
      {iPim, iKp, 4},     // 11 pi- K+
      {iPip, iKp, 4},     // 12 pi+ K+
      {iPim, iKm, 4},     // 13 pi- K-
      {iPip, iPbar, 4},   // 14 pi+ pbar
      {iPim, iP, 4},      // 15 pi- p (Lambda, Delta0)
      {iPip, iP, 4},      // 16 pi+ p (Delta++)
      {iPim, iPbar, 4},   // 17 pi- pbar
      {iKp, iPbar, 10},   // 18 K+ pbar
      {iKm, iP, 10},      // 19 K- p
      {iKp, iP, 10},      // 20 K+ p
      {iKm, iPbar, 10},   // 21 K- pbar
      {iKp, iLam, 20},    // 22 K+ Lambda
      {iKp, iLbar, 20},   // 23 K+ Lbar
      {iKm, iLam, 20},    // 24 K- Lambda
      {iKm, iLbar, 20},   // 25 K- Lbar
      {iP, iLam, 20},     // 26 p Lambda
      {iP, iLbar, 20},    // 27 p Lbar
      {iPbar, iLam, 20},  // 28 pbar Lambda
      {iPbar, iLbar, 20}, // 29 pbar Lbar
      {iPip, iK0s, 5},    // 30 pi+ K0s (K*+)
      {iPim, iK0s, 5},    // 31 pi- K0s (K*-)
      {iPip, iLam, 5},    // 32 pi+ Lambda (Sigma*+)
      {iPim, iLam, 5},    // 33 pi- Lambda (Sigma*-, Xi-)
      {iPip, iLbar, 5},   // 34 pi+ Lbar
      {iPim, iLbar, 5},   // 35 pi- Lbar
      {iKp, iK0s, 20},    // 36 K+ K0s
      {iKm, iK0s, 20},    // 37 K- K0s
      {iP, iK0s, 20},     // 38 p K0s
      {iPbar, iK0s, 20},  // 39 pbar K0s
      {iK0s, iK0s, 35},   // 40 K0s K0s
      {iLam, iLam, 35},   // 41 Lambda Lambda
      {iLbar, iLbar, 35}, // 42 Lbar Lbar
      {iK0s, iLam, 35},   // 43 K0s Lambda
      {iK0s, iLbar, 35},  // 44 K0s Lbar
      {iLam, iLbar, 35},  // 45 Lambda Lbar
      {iHp, iHm, 4},      // 46 h+ h-   (|eta| < 1.1)
      {iHp, iHp, 4},      // 47 h+ h+   (|eta| < 1.1)
      {iHm, iHm, 4}};     // 48 h- h-   (|eta| < 1.1)
  return p;
}

// true when both legs have maps in this module (no V0 leg)
inline bool PairMeasurable(int ip)
{
  const auto &s = AllSpecies();
  const PairType &p = PairTypes()[ip];
  return s[p.s1].kind != kV0 && s[p.s2].kind != kV0;
}

// Corral's startup rule: each species window (full width 2 absMax) holds a whole
// number (>= 2) of bins of the pair's width. Returns the number of bad pair types.
inline int CheckPairTypes(bool verbose = true)
{
  int bad = 0;
  const auto &s = AllSpecies();
  for (int ip = 0; ip < NPairTypes; ++ip)
  {
    const PairType &p = PairTypes()[ip];
    for (int leg : {p.s1, p.s2})
    {
      const double nb = 2. * s[leg].absMax / (0.01 * p.bw);
      if (std::fabs(nb - std::lround(nb)) > 1e-6 || std::lround(nb) < 2)
      {
        ++bad;
        if (verbose)
          std::printf("EffTypes: pair type %d (%s %s): bin width %.2f does not divide the %s window %.2f\n", ip,
                      s[p.s1].tag, s[p.s2].tag, 0.01 * p.bw, s[leg].tag, 2. * s[leg].absMax);
      }
    }
  }
  return bad;
}
}  // namespace efftypes

#endif
