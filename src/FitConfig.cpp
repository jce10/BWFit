#include "FitConfig.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>

namespace bwfit {

namespace {
std::string label;

template <class T>
void ReadLabeled(std::istream& in, T& value) {
  in >> label >> value;
  if (!in) throw std::runtime_error("Failed while reading value after label: " + label);
}

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}
}

std::string ToString(BackgroundType type) {
  switch (type) {
    case BackgroundType::Linear: return "linear";
    case BackgroundType::Quadratic: return "quadratic";
  }
  return "unknown";
}

BackgroundType ParseBackgroundType(const std::string& text) {
  const std::string s = Lower(text);
  if (s == "linear" || s == "lin" || s == "1") return BackgroundType::Linear;
  if (s == "quadratic" || s == "quad" || s == "parabolic" || s == "2") return BackgroundType::Quadratic;
  throw std::runtime_error("Unknown BACKGROUND_TYPE: " + text + " (expected linear or quadratic)");
}

FitConfig ReadConfig(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open config file: " + path);

  FitConfig cfg;
  int num_states = 0;
  int rel = 0;
  int interf = 1;
  int fit_phase = 1;

  // read config file line by line
  ReadLabeled(input, cfg.angle_deg);
  input >> label >> cfg.rootfile;
  input >> label >> cfg.output_rootfile;
  input >> label >> cfg.tree_name;
  input >> label >> cfg.branch_name;
  input >> label >> cfg.title;
  ReadLabeled(input, cfg.calib_slope);
  ReadLabeled(input, cfg.calib_intercept);

  ReadLabeled(input, num_states); 
  ReadLabeled(input, rel);
  ReadLabeled(input, interf);
  ReadLabeled(input, fit_phase);
  ReadLabeled(input, cfg.width_perc);
  ReadLabeled(input, cfg.trapped_index);
  ReadLabeled(input, cfg.superrad_index);

  ReadLabeled(input, cfg.fullmin);
  ReadLabeled(input, cfg.fullmax);
  ReadLabeled(input, cfg.fitmin);
  ReadLabeled(input, cfg.fitmax);
  ReadLabeled(input, cfg.fit_bins);
  ReadLabeled(input, cfg.graph_points);

  cfg.use_relativistic = (rel != 0); // relativistic 0 = off, 1 = on
  cfg.use_interference = (interf != 0); // interference 0 = off, 1 = on
  cfg.fit_phase = (fit_phase != 0); // phase 0 = fixed, 1 = floating

  ReadLabeled(input, cfg.phase.value);
  ReadLabeled(input, cfg.phase.min);
  ReadLabeled(input, cfg.phase.max);
  cfg.phase.value *= M_PI / 2.0;
  cfg.phase.min *= M_PI / 2.0;
  cfg.phase.max *= M_PI / 2.0;

  std::string bg_type;
  input >> label >> bg_type;
  if (!input) throw std::runtime_error("Failed while reading BACKGROUND_TYPE");
  cfg.background_type = ParseBackgroundType(bg_type);

  ReadLabeled(input, cfg.bg_e.value);
  ReadLabeled(input, cfg.bg_e.min);
  ReadLabeled(input, cfg.bg_e.max);

  ReadLabeled(input, cfg.bg_a0.value);
  ReadLabeled(input, cfg.bg_a0.min);
  ReadLabeled(input, cfg.bg_a0.max);

  ReadLabeled(input, cfg.bg_a1.value);
  ReadLabeled(input, cfg.bg_a1.min);
  ReadLabeled(input, cfg.bg_a1.max);

  if (cfg.background_type == BackgroundType::Quadratic) {
    ReadLabeled(input, cfg.bg_a2.value);
    ReadLabeled(input, cfg.bg_a2.min);
    ReadLabeled(input, cfg.bg_a2.max);
  } else {
    cfg.bg_a2.value = 0.0;
    cfg.bg_a2.min = 0.0;
    cfg.bg_a2.max = 0.0;
  }

  cfg.states.resize(num_states);
  std::string state_label;
  input >> state_label;
  for (int i = 0; i < num_states; ++i) {
    const std::string expected = "STATE" + std::to_string(i);
    // Preserve positional state indexing used by old cards, including cards
    // whose printed STATE labels are duplicated or out of sequence.
    if (state_label.compare(0, 5, "STATE") != 0)
      throw std::runtime_error("Expected a STATE block, found " + state_label);

    auto& st = cfg.states[i];
    st.M_range = -1.0; // Omitted: retain the legacy full-fit-window limits.
    std::set<std::string> seen;
    state_label.clear();
    std::string key;
    while (input >> key) {
      if (key[0] == '#') {
        std::string comment;
        std::getline(input, comment);
        continue;
      }
      if (key.compare(0, 5, "STATE") == 0) {
        state_label = key;
        break;
      }
      if (!key.empty() && key.back() == ':') key.pop_back();
      if (key != "FAC" && key != "M" && key != "M_RANGE" && key != "G" &&
          key != "IS_BW" && key != "FIT_M" && key != "FIT_G")
        throw std::runtime_error(expected + ": unknown field " + key);
      if (!seen.insert(key).second)
        throw std::runtime_error(expected + ": duplicate field " + key);

      std::string token;
      if (!(input >> token))
        throw std::runtime_error(expected + ": missing value for " + key);
      std::size_t used = 0;
      double value = 0.0;
      try { value = std::stod(token, &used); }
      catch (const std::exception&) {
        throw std::runtime_error(expected + ": invalid value for " + key + ": " + token);
      }
      if (used != token.size() || !std::isfinite(value))
        throw std::runtime_error(expected + ": invalid value for " + key + ": " + token);
      if (key == "FAC") st.fac = value;
      else if (key == "M") st.M = value;
      else if (key == "G") st.G = value;
      else if (key == "M_RANGE") {
        if (value <= 0.0)
          throw std::runtime_error(expected + ": M_RANGE must be positive (MeV); omit it for legacy limits");
        st.M_range = value;
      }
      else {
        if (std::trunc(value) != value)
          throw std::runtime_error(expected + ": " + key + " must be an integer");
        if (key == "IS_BW") st.is_bw = (value != 0.0);
        else if (key == "FIT_M") st.fit_M = (value != 0.0);
        else st.fit_G = (value != 0.0);
      }
    }
    for (const char* required : {"FAC", "M", "G", "IS_BW", "FIT_M", "FIT_G"})
      if (!seen.count(required))
        throw std::runtime_error(expected + ": missing field " + required);
    if (st.fit_M && st.M_range > 0.0) {
      const double lo = std::max(cfg.fitmin, st.M - st.M_range);
      const double hi = std::min(cfg.fitmax, st.M + st.M_range);
      if (!(lo < hi) || st.M < lo || st.M > hi)
        throw std::runtime_error(expected + ": centroid range/start is incompatible with FITMIN/FITMAX");
    }
  }
  if (!state_label.empty())
    throw std::runtime_error("Unexpected extra state: " + state_label);

  return cfg;
}

void PrintConfig(const FitConfig& cfg) {
  std::cout << "Input read:\n"
            << "  ANGLE = " << cfg.angle_deg << " deg\n"
            << "  ROOTFILE = " << cfg.rootfile << "\n"
            << "  TREE/BRANCH = " << cfg.tree_name << "/" << cfg.branch_name << "\n"
            << "  NUM_STATES = " << cfg.states.size() << " (numparams = " << cfg.NumParams() << ")\n"
            << "  USE_RELATIVISTIC = " << cfg.use_relativistic << "\n"
            << "  USE_INTERFERENCE = " << cfg.use_interference << "\n"
            << "  FIT_PHASE = " << cfg.fit_phase << "\n"
            << "  TI/SI = " << cfg.trapped_index << "/" << cfg.superrad_index << "\n"
            << "  WIDTH_PERC = " << 100.0 * cfg.width_perc << "%\n"
            << "  FULL = [" << cfg.fullmin << ", " << cfg.fullmax << "] MeV\n"
            << "  FIT  = [" << cfg.fitmin << ", " << cfg.fitmax << "] MeV\n"
            << "  PHASE = " << cfg.phase.value * 2.0 / M_PI
            << " * pi/2 [" << (cfg.fit_phase ? "FLOAT" : "FIXED") << "]\n"
            << "  BACKGROUND_TYPE = " << ToString(cfg.background_type) << "\n";

  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    const auto& st = cfg.states[i];
    std::cout << "  STATE" << i << " centroid: ";
    if (!st.fit_M) {
      std::cout << "FIXED at " << st.M << " MeV\n";
    } else {
      const double lo = st.M_range > 0.0 ? std::max(cfg.fitmin, st.M - st.M_range) : cfg.fitmin;
      const double hi = st.M_range > 0.0 ? std::min(cfg.fitmax, st.M + st.M_range) : cfg.fitmax;
      std::cout << "FLOAT in [" << lo << ", " << hi << "] MeV";
      if (st.M_range > 0.0) std::cout << " (M_RANGE = +/-" << st.M_range << " MeV)";
      std::cout << "\n";
    }
  }

  if (cfg.background_type == BackgroundType::Linear) {
    std::cout << "  BG(E)=A0+A1(E-E0)\n"
              << "    E0=" << cfg.bg_e.value << " A0=" << cfg.bg_a0.value
              << " A1=" << cfg.bg_a1.value << std::endl << std::endl;
  } else {
    std::cout << "  BG(E)=A0+A1(E-E0)+A2(E-E0)^2\n"
              << "    E0=" << cfg.bg_e.value << " A0=" << cfg.bg_a0.value
              << " A1=" << cfg.bg_a1.value << " A2=" << cfg.bg_a2.value << std::endl << std::endl;
  }
}

} // namespace bwfit
