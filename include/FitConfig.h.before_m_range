#pragma once

#include <string>
#include <vector>

namespace bwfit {

struct StateConfig {
  double fac = 0.0;
  double M = 0.0;
  double G = 0.0;
  bool is_bw = false;
  bool fit_M = false;
  bool fit_G = false;
};

struct BoundedValue {
  double value = 0.0;
  double min = 0.0;
  double max = 0.0;
};

enum class BackgroundType {
  Linear,
  Quadratic
};

std::string ToString(BackgroundType type);
BackgroundType ParseBackgroundType(const std::string& text);

struct FitConfig {
  double angle_deg = 0.0;
  std::string rootfile;
  std::string output_rootfile;
  std::string tree_name = "singles_tree";
  std::string branch_name = "singles_fppos";

  std::string title;
  double calib_slope = 1.0;      // MeV per branch unit
  double calib_intercept = 0.0;  // MeV

  int fit_bins = 250;
  int graph_points = 5000;
  double fullmin = 0.0;
  double fullmax = 20.0;
  double fitmin = 0.0;
  double fitmax = 20.0;

  bool use_relativistic = false;
  bool use_interference = true;
  bool fit_phase = true;
  double width_perc = 0.0;
  int trapped_index = 4;
  int superrad_index = 5;

  BoundedValue phase; // input as multiples of pi/2 in config; converted to radians by reader

  BackgroundType background_type = BackgroundType::Quadratic;
  BoundedValue bg_e;   // Center energy E0
  BoundedValue bg_a0;  // Constant term
  BoundedValue bg_a1;  // Linear term
  BoundedValue bg_a2;  // Quadratic term; used only for BackgroundType::Quadratic

  std::vector<StateConfig> states;

  int NumBackgroundShapeParams() const {
    return background_type == BackgroundType::Quadratic ? 4 : 3; // E0,A0,A1,(A2)
  }

  int NumParams() const {
    return 3 * static_cast<int>(states.size()) + 1 + NumBackgroundShapeParams();
  }

  int PhaseIndex() const { return 3 * static_cast<int>(states.size()); }
  int BgEIndex() const { return PhaseIndex() + 1; }
  int BgA0Index() const { return PhaseIndex() + 2; }
  int BgA1Index() const { return PhaseIndex() + 3; }
  int BgA2Index() const {
    return background_type == BackgroundType::Quadratic ? PhaseIndex() + 4 : -1;
  }
};

FitConfig ReadConfig(const std::string& path);
void PrintConfig(const FitConfig& cfg);

} // namespace bwfit
