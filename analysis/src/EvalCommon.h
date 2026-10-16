// =============================================================================
// EvalCommon.h - utilities shared by the programs that read SvtxEvaluator files
//
// The helpers below are the ones used by VtxEffAnalysis.cc, made reusable:
//   Trim, key = value configuration reader   (VtxEffAnalysis.cc ReadConfig)
//   ChargeLookup -> PdgLookup (charge + mass)(VtxEffAnalysis.cc ChargeLookup)
//   Bind                                      (VtxEffAnalysis.cc Bind)
//   JobFromName                               (VtxEffAnalysis.cc JobFromName)
//   Finite, ReadFileList                      (VtxEffAnalysis.cc main)
// The logic is unchanged. VtxEffAnalysis.cc itself is deliberately left
// byte-identical (it does not include this header yet), so the validated
// vertex-efficiency results cannot change. Switching it to this header is a
// possible later step that would need a regression check.
//
// Needs ROOT only (-lEG for TDatabasePDG).
// =============================================================================
#ifndef EVALCOMMON_H
#define EVALCOMMON_H

#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace evalcommon
{

// Trim leading/trailing white space.
inline std::string Trim(const std::string &s)
{
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

// -----------------------------------------------------------------------------
// key = value configuration, '#' starts a comment. Every key must be declared
// with a default before Read(); unknown keys in the file are fatal (typos can
// never silently fall back to a default).
// -----------------------------------------------------------------------------
class KeyValueConfig
{
 public:
  void Declare(const std::string &key, const std::string &def) { v_[key] = def; }

  void Read(const std::string &path)
  {
    std::ifstream in(path);
    if (!in)
    {
      std::cerr << "FATAL: cannot open config file " << path << std::endl;
      std::exit(1);
    }
    std::string line;
    int lineno = 0;
    while (std::getline(in, line))
    {
      ++lineno;
      const auto hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      line = Trim(line);
      if (line.empty()) continue;
      const auto eq = line.find('=');
      if (eq == std::string::npos)
      {
        std::cerr << "FATAL: " << path << ":" << lineno << " has no '='" << std::endl;
        std::exit(1);
      }
      const std::string key = Trim(line.substr(0, eq));
      if (!v_.count(key))
      {
        std::cerr << "FATAL: unknown config key '" << key << "' (" << path << ":" << lineno << ")" << std::endl;
        std::exit(1);
      }
      v_[key] = Trim(line.substr(eq + 1));
    }
  }

  void Set(const std::string &key, const std::string &val)
  {
    if (!v_.count(key))
    {
      std::cerr << "FATAL: unknown config key '" << key << "'" << std::endl;
      std::exit(1);
    }
    v_[key] = val;
  }
  std::string Str(const std::string &k) const { return v_.at(k); }
  double Dbl(const std::string &k) const { return std::stod(v_.at(k)); }
  long Int(const std::string &k) const { return std::stol(v_.at(k)); }
  bool Bool(const std::string &k) const { return std::stol(v_.at(k)) != 0; }
  // all keys and values, one "key = value" per line (stored in the output)
  std::string Dump() const
  {
    std::string s;
    for (const auto &kv : v_) s += kv.first + " = " + kv.second + "\n";
    return s;
  }

 private:
  std::map<std::string, std::string> v_;
};

// -----------------------------------------------------------------------------
// PDG code -> electric charge (units of e) and mass (GeV/c^2)
// -----------------------------------------------------------------------------
// The evaluator stores the PDG code as a Float_t (gflavor). A float holds
// integers exactly only up to 2^24 = 16,777,216. Ordinary hadron/lepton codes
// are far below that, but nuclear codes 10LZZZAAAI (~1e9) are NOT exact:
// e.g. 1000010020 (deuteron) can come back off by up to +-64. Therefore:
//   - round the float to the nearest integer;
//   - for nuclei (|pdg| >= 1e9), take Z from digits ZZZ after rounding to the
//     nearest multiple of 10^4, which is immune to the +-64 float error:
//       Z = ((|pdg| + 5000) / 10000) % 1000
//     The mass number A (digits AAA) is NOT recoverable, so the mass of a
//     nucleus is reported as unknown (NaN);
//   - otherwise ask ROOT's TDatabasePDG (charge is stored in units of |e|/3).
// Unknown codes: charge 0, mass NaN, known = false.
struct PdgInfo
{
  int charge = 0;
  double mass = NAN;
  bool known = false;   // charge known
  bool nucleus = false;
};

class PdgLookup
{
 public:
  const PdgInfo &Get(long long pdg)
  {
    auto it = cache_.find(pdg);
    if (it != cache_.end()) return it->second;
    PdgInfo p;
    const long long apdg = std::llabs(pdg);
    if (apdg >= 1000000000LL)
    {
      const long long Z = ((apdg + 5000LL) / 10000LL) % 1000LL;
      p.charge = static_cast<int>(pdg > 0 ? Z : -Z);
      p.known = true;
      p.nucleus = true;
    }
    else
    {
      TParticlePDG *d = TDatabasePDG::Instance()->GetParticle(static_cast<int>(pdg));
      if (d)
      {
        p.charge = static_cast<int>(std::lround(d->Charge() / 3.0));
        p.mass = d->Mass();
        p.known = true;
      }
    }
    return cache_.emplace(pdg, p).first->second;
  }
  int Charge(long long pdg, bool &known)
  {
    const PdgInfo &p = Get(pdg);
    known = p.known;
    return p.charge;
  }

 private:
  std::unordered_map<long long, PdgInfo> cache_;
};

// Rapidity from momentum and mass: y = asinh(pz / mT), mT = sqrt(m^2 + pT^2).
// (Equivalent to 0.5 ln((E+pz)/(E-pz)) but numerically stable.) NaN if the
// mass is unknown.
inline double Rapidity(double px, double py, double pz, double mass)
{
  if (!std::isfinite(mass)) return NAN;
  const double mt = std::sqrt(mass * mass + px * px + py * py);
  if (mt <= 0) return NAN;
  return std::asinh(pz / mt);
}

// -----------------------------------------------------------------------------
// TNtuple helpers
// -----------------------------------------------------------------------------
inline bool Finite(float x) { return std::isfinite(x); }

// Enable only the listed branches and bind them to float variables. Exit with
// an error if a branch is missing (better than silently reading zeros).
inline void Bind(TTree *t, const std::vector<std::pair<const char *, float *>> &br)
{
  t->SetBranchStatus("*", 0);
  for (const auto &b : br)
  {
    if (!t->GetBranch(b.first))
    {
      std::cerr << "FATAL: branch " << b.first << " missing in tree " << t->GetName() << std::endl;
      std::exit(2);
    }
    t->SetBranchStatus(b.first, 1);
    t->SetBranchAddress(b.first, b.second);
  }
}

// Bind many branches by name into one float buffer; returns name -> index.
// Use the index (resolved once) inside event loops.
class FloatRow
{
 public:
  FloatRow(TTree *t, const std::vector<std::string> &names) : buf_(names.size(), NAN), names_(names)
  {
    std::vector<std::pair<const char *, float *>> br;
    for (size_t i = 0; i < names_.size(); ++i)
    {
      idx_[names_[i]] = static_cast<int>(i);
      br.emplace_back(names_[i].c_str(), &buf_[i]);
    }
    Bind(t, br);
  }
  int Index(const std::string &n) const
  {
    auto it = idx_.find(n);
    if (it == idx_.end())
    {
      std::cerr << "FATAL: FloatRow has no column " << n << std::endl;
      std::exit(2);
    }
    return it->second;
  }
  float operator[](int i) const { return buf_[i]; }

 private:
  std::vector<float> buf_;
  std::vector<std::string> names_;
  std::map<std::string, int> idx_;
};

// Extract the job number from ".../outputeval_NNNNNN.root". -1 if not found.
inline int JobFromName(const std::string &path)
{
  static const std::regex re("outputeval_([0-9]+)\\.root$");
  std::smatch m;
  if (std::regex_search(path, m, re)) return std::stoi(m[1]);
  return -1;
}

// Comma- or space-separated list of numbers, e.g. "0.1, 0.2, 0.5". Exits on
// a malformed entry (a silent default would hide a typo in a bin edge).
inline std::vector<double> ParseNumberList(const std::string &text, const std::string &what)
{
  std::vector<double> out;
  std::string token;
  auto flush = [&]() {
    const std::string t = Trim(token);
    token.clear();
    if (t.empty()) return;
    try
    {
      size_t used = 0;
      out.push_back(std::stod(t, &used));
      if (used != t.size()) throw std::invalid_argument(t);
    }
    catch (const std::exception &)
    {
      std::cerr << "FATAL: '" << t << "' in " << what << " is not a number" << std::endl;
      std::exit(1);
    }
  };
  for (char c : text)
  {
    if (c == ',' || c == ' ' || c == '\t') flush();
    else token += c;
  }
  flush();
  return out;
}

// Parse a "key = value" text block (one pair per line), e.g. the configuration
// that EvalSkim stores in every skim file as TNamed "config".
inline std::map<std::string, std::string> ParseKeyValueText(const std::string &text)
{
  std::map<std::string, std::string> kv;
  size_t start = 0;
  while (start <= text.size())
  {
    size_t end = text.find('\n', start);
    if (end == std::string::npos) end = text.size();
    const std::string line = text.substr(start, end - start);
    const auto eq = line.find('=');
    if (eq != std::string::npos) kv[Trim(line.substr(0, eq))] = Trim(line.substr(eq + 1));
    start = end + 1;
  }
  return kv;
}

// Plain text list: one file per line ('#' comments); relative names are
// prefixed with dir (if not empty).
inline std::vector<std::string> ReadFileList(const std::string &list, const std::string &dir)
{
  std::vector<std::string> files;
  std::ifstream in(list);
  if (!in)
  {
    std::cerr << "FATAL: cannot open input list " << list << std::endl;
    std::exit(1);
  }
  std::string line;
  while (std::getline(in, line))
  {
    const auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = Trim(line);
    if (line.empty()) continue;
    files.push_back(line[0] == '/' || dir.empty() ? line : dir + "/" + line);
  }
  return files;
}

}  // namespace evalcommon

#endif
