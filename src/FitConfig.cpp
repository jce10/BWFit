#include "FitConfig.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
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
  for (int i = 0; i < num_states; ++i) {
    int is_bw = 0, fit_M = 0, fit_G = 0;
    input >> label; // STATEi
    ReadLabeled(input, cfg.states[i].fac);
    ReadLabeled(input, cfg.states[i].M);
    ReadLabeled(input, cfg.states[i].G);
    ReadLabeled(input, is_bw);
    ReadLabeled(input, fit_M);
    ReadLabeled(input, fit_G);
    cfg.states[i].is_bw = (is_bw != 0);
    cfg.states[i].fit_M = (fit_M != 0);
    cfg.states[i].fit_G = (fit_G != 0);
  }

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
            << "  TI/SI = " << cfg.trapped_index << "/" << cfg.superrad_index << "\n"
            << "  WIDTH_PERC = " << 100.0 * cfg.width_perc << "%\n"
            << "  FULL = [" << cfg.fullmin << ", " << cfg.fullmax << "] MeV\n"
            << "  FIT  = [" << cfg.fitmin << ", " << cfg.fitmax << "] MeV\n"
            << "  PHASE = " << cfg.phase.value * 2.0 / M_PI << " * pi/2\n"
            << "  BACKGROUND_TYPE = " << ToString(cfg.background_type) << "\n";

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
