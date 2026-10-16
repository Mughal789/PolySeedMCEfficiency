// =============================================================================
// TrkEffConfigUtil.h - interpretation of the TrkEffAnalysis "config" record
// shared by the macros (one place for rules that changed between versions)
//
// EffectiveTruthEtaMax: the truth |eta| threshold actually applied, <= 0 = no
// truth-eta cut. Same rule as TrkEffAnalysis v3 (ReadTrkEffConfig):
//   use_truth_eta_cut absent or "auto" -> legacy: truth_abseta_max (> 0 = cut on)
//   use_truth_eta_cut = 1              -> truth_eta_cut_max
//   use_truth_eta_cut = 0              -> no cut
// =============================================================================
#ifndef TRKEFF_CONFIG_UTIL_H
#define TRKEFF_CONFIG_UTIL_H

#include <map>
#include <string>

namespace trkeffcfg
{
inline double EffectiveTruthEtaMax(const std::map<std::string, std::string> &cfg)
{
  auto get = [&](const char *k, const char *def) {
    auto it = cfg.find(k);
    return it == cfg.end() ? std::string(def) : it->second;
  };
  const std::string mode = get("use_truth_eta_cut", "auto");
  if (mode == "auto") return std::stod(get("truth_abseta_max", "1.1"));
  return mode == "1" ? std::stod(get("truth_eta_cut_max", "1.1")) : 0.0;
}
}  // namespace trkeffcfg

#endif
