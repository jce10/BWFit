#include "FitResults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


namespace bwfit {

namespace {

std::string JsonEscape(const std::string& text) {
  std::string out;
  out.reserve(text.size());

  for (char c : text) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\b': out += "\\b"; break;
      case '\f': out += "\\f"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out += c; break;
    }
  }

  return out;
}

const char* BoolText(bool value) {
  return value ? "true" : "false";
}

void WriteValueBoundsJSON(std::ostream& out,
                          const std::string& name,
                          double value,
                          double min,
                          double max,
                          int indent = 4) {
  const std::string pad(indent, ' ');
  out << pad << "\"" << name << "\": {\n"
      << pad << "  \"value\": " << value << ",\n"
      << pad << "  \"min\": " << min << ",\n"
      << pad << "  \"max\": " << max << "\n"
      << pad << "}";
}

void WriteFitParamJSON(std::ostream& out,
                       const std::string& name,
                       double set,
                       double fit,
                       double err,
                       const std::string& status,
                       int indent = 6) {
  const std::string pad(indent, ' ');
  out << pad << "\"" << name << "\": {\n"
      << pad << "  \"set\": " << set << ",\n"
      << pad << "  \"fit\": " << fit << ",\n"
      << pad << "  \"err\": " << err << ",\n"
      << pad << "  \"status\": \"" << status << "\"\n"
      << pad << "}";
}

} // namespace

void PrintFitParameters(const TF1& fit, const FitConfig& cfg) {

  std::cout << "\n=== Fit Parameters ===\n";
  std::cout << std::fixed << std::showpoint << std::setprecision(4);

  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    const int i_fac = 3 * i + 0;
    const int i_M   = 3 * i + 1;
    const int i_G   = 3 * i + 2;

    // state config info
    const auto& st = cfg.states[i];
    const char* M_status = st.fit_M ? "FLOAT" : "FIXED";
    const char* G_status = st.fit_G ? "FLOAT" : "FIXED";


    std::cout << "State " << i
              << (st.is_bw ? " (Breit-Wigner)" : " (Gaussian)") << ":\n";

    std::cout << "  Yield: "
              << "set = " << st.fac
              << " | fit = " << fit.GetParameter(i_fac)
              << " +/- " << fit.GetParError(i_fac) << "\n";

    std::cout << "  M:     "
              << "set = " << st.M << " MeV"
              << " | fit = " << fit.GetParameter(i_M)
              << " +/- " << fit.GetParError(i_M) << " MeV"
              << " [" << M_status << "]" << "\n";

    std::cout << "  G:     "
              << "set = " << st.G << " MeV"
              << " | fit = " << fit.GetParameter(i_G)
              << " +/- " << fit.GetParError(i_G) << " MeV"
              << " [" << G_status << "]" << "\n\n";
  }

  // std::cout << "\nInterference phase, δ:\n";
  // std::cout << " | fit = " << fit.GetParameter(cfg.PhaseIndex())
  //           << " +/- " << fit.GetParError(cfg.PhaseIndex()) << "\n";

  // Interference phase, δ, parameters
  const double phase_set_rad = cfg.phase.value;
  const double phase_set_pi2 = phase_set_rad * 2.0 / M_PI;

  const double phase_fit_rad = fit.GetParameter(cfg.PhaseIndex());
  const double phase_err_rad = fit.GetParError(cfg.PhaseIndex());

  const double phase_fit_pi2 = phase_fit_rad * 2.0 / M_PI;
  const double phase_err_pi2 = phase_err_rad * 2.0 / M_PI;
  const char* phase_status = cfg.fit_phase ? "FLOAT" : "FIXED";

  std::cout << "\nInterference phase, delta:\n";
  std::cout << "  set = " << phase_set_rad
            << " rad = " << phase_set_pi2 << " * pi/2"
            << " [" << phase_status << "]\n";

  std::cout << "  fit = " << phase_fit_rad
            << " +/- " << phase_err_rad << " rad\n";

  std::cout << "      = " << phase_fit_pi2
            << " +/- " << phase_err_pi2 << " * pi/2\n";


  // Background parameters
  std::cout << "\nBackground (" << ToString(cfg.background_type) << "):\n";

  std::cout << "  E0: set = " << cfg.bg_e.value
            << " | fit = " << fit.GetParameter(cfg.BgEIndex())
            << " +/- " << fit.GetParError(cfg.BgEIndex()) << "\n";

  std::cout << "  A0: set = " << cfg.bg_a0.value
            << " | fit = " << fit.GetParameter(cfg.BgA0Index())
            << " +/- " << fit.GetParError(cfg.BgA0Index()) << "\n";

  std::cout << "  A1: set = " << cfg.bg_a1.value
            << " | fit = " << fit.GetParameter(cfg.BgA1Index())
            << " +/- " << fit.GetParError(cfg.BgA1Index()) << "\n";

  if (cfg.background_type == BackgroundType::Quadratic) {
    std::cout << "  A2: set = " << cfg.bg_a2.value
              << " | fit = " << fit.GetParameter(cfg.BgA2Index())
              << " +/- " << fit.GetParError(cfg.BgA2Index()) << "\n";
  }

  const double chi2 = fit.GetChisquare();
  const int ndf = fit.GetNDF();

  std::cout << "\nFit quality:\n";
  std::cout << "  chi2/NDF = " << chi2 << " / " << ndf;
  if (ndf > 0) std::cout << " = " << chi2 / ndf;
  std::cout << "\n\n";
}

void WriteFitSummaryJSON(const TF1& fit,
                         const FitConfig& cfg,
                         const std::string& output_json,
                         int fit_status,
                         const std::string& config_path) {
  namespace fs = std::filesystem;

  fs::path outpath(output_json);
  if (outpath.extension().empty()) outpath += ".json";
  if (outpath.has_parent_path()) fs::create_directories(outpath.parent_path());

  std::ofstream out(outpath);
  if (!out) {
    std::cerr << "[error] Could not open JSON file: " << outpath << "\n";
    return;
  }

  out << std::fixed << std::showpoint << std::setprecision(10);

  const double chi2 = fit.GetChisquare();
  const int ndf = fit.GetNDF();
  const double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;

  const double phase_set_rad = cfg.phase.value;
  const double phase_set_pi2 = phase_set_rad * 2.0 / M_PI;
  const double phase_fit_rad = fit.GetParameter(cfg.PhaseIndex());
  const double phase_err_rad = fit.GetParError(cfg.PhaseIndex());
  const double phase_fit_pi2 = phase_fit_rad * 2.0 / M_PI;
  const double phase_err_pi2 = phase_err_rad * 2.0 / M_PI;

  out << "{\n";

  out << "  \"fit_config\": {\n"
      << "    \"source_config_path\": \"" << JsonEscape(config_path) << "\",\n"
      << "    \"angle_deg\": " << cfg.angle_deg << ",\n"
      << "    \"rootfile\": \"" << JsonEscape(cfg.rootfile) << "\",\n"
      << "    \"output_rootfile\": \"" << JsonEscape(cfg.output_rootfile) << "\",\n"
      << "    \"tree_name\": \"" << JsonEscape(cfg.tree_name) << "\",\n"
      << "    \"branch_name\": \"" << JsonEscape(cfg.branch_name) << "\",\n"
      << "    \"title\": \"" << JsonEscape(cfg.title) << "\",\n"
      << "    \"calibration\": {\n"
      << "      \"slope\": " << cfg.calib_slope << ",\n"
      << "      \"intercept\": " << cfg.calib_intercept << "\n"
      << "    },\n"
      << "    \"model_settings\": {\n"
      << "      \"num_states\": " << cfg.states.size() << ",\n"
      << "      \"num_params\": " << cfg.NumParams() << ",\n"
      << "      \"use_relativistic\": " << BoolText(cfg.use_relativistic) << ",\n"
      << "      \"use_interference\": " << BoolText(cfg.use_interference) << ",\n"
      << "      \"fit_phase\": " << BoolText(cfg.fit_phase) << ",\n"
      << "      \"width_perc\": " << cfg.width_perc << ",\n"
      << "      \"trapped_index\": " << cfg.trapped_index << ",\n"
      << "      \"superrad_index\": " << cfg.superrad_index << "\n"
      << "    },\n"
      << "    \"ranges\": {\n"
      << "      \"full_min_MeV\": " << cfg.fullmin << ",\n"
      << "      \"full_max_MeV\": " << cfg.fullmax << ",\n"
      << "      \"fit_min_MeV\": " << cfg.fitmin << ",\n"
      << "      \"fit_max_MeV\": " << cfg.fitmax << ",\n"
      << "      \"fit_bins\": " << cfg.fit_bins << ",\n"
      << "      \"graph_points\": " << cfg.graph_points << "\n"
      << "    },\n"
      << "    \"phase\": {\n"
      << "      \"value_rad\": " << cfg.phase.value << ",\n"
      << "      \"min_rad\": " << cfg.phase.min << ",\n"
      << "      \"max_rad\": " << cfg.phase.max << ",\n"
      << "      \"value_pi_over_2_units\": " << cfg.phase.value * 2.0 / M_PI << ",\n"
      << "      \"min_pi_over_2_units\": " << cfg.phase.min * 2.0 / M_PI << ",\n"
      << "      \"max_pi_over_2_units\": " << cfg.phase.max * 2.0 / M_PI << "\n"
      << "    },\n"
      << "    \"background\": {\n"
      << "      \"type\": \"" << ToString(cfg.background_type) << "\",\n";

  WriteValueBoundsJSON(out, "E0_MeV", cfg.bg_e.value, cfg.bg_e.min, cfg.bg_e.max, 6);
  out << ",\n";
  WriteValueBoundsJSON(out, "A0", cfg.bg_a0.value, cfg.bg_a0.min, cfg.bg_a0.max, 6);
  out << ",\n";
  WriteValueBoundsJSON(out, "A1", cfg.bg_a1.value, cfg.bg_a1.min, cfg.bg_a1.max, 6);
  if (cfg.background_type == BackgroundType::Quadratic) {
    out << ",\n";
    WriteValueBoundsJSON(out, "A2", cfg.bg_a2.value, cfg.bg_a2.min, cfg.bg_a2.max, 6);
  }
  out << "\n    },\n";

  out << "    \"states\": [\n";
  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    const auto& st = cfg.states[i];

    out << "      {\n"
        << "        \"index\": " << i << ",\n"
        << "        \"shape\": \"" << (st.is_bw ? "Breit-Wigner" : "Gaussian") << "\",\n"
        << "        \"is_bw\": " << BoolText(st.is_bw) << ",\n"
        << "        \"yield\": { \"value\": " << st.fac << " },\n"
        << "        \"centroid_MeV\": { \"value\": " << st.M
        << ", \"vary\": " << BoolText(st.fit_M) << " },\n"
        << "        \"width_MeV\": { \"value\": " << st.G
        << ", \"vary\": " << BoolText(st.fit_G) << " }\n"
        << "      }" << (i + 1 < static_cast<int>(cfg.states.size()) ? "," : "") << "\n";
  }
  out << "    ]\n"
      << "  },\n";

  out << "  \"fit_result\": {\n"
      << "    \"status\": " << fit_status << ",\n"
      << "    \"quality\": {\n"
      << "      \"chi2\": " << chi2 << ",\n"
      << "      \"ndf\": " << ndf << ",\n"
      << "      \"chi2_ndf\": " << chi2_ndf << "\n"
      << "    },\n"
      << "    \"states\": [\n";

  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    const int i_fac = 3 * i + 0;
    const int i_M   = 3 * i + 1;
    const int i_G   = 3 * i + 2;

    const auto& st = cfg.states[i];
    const char* M_status = st.fit_M ? "FLOAT" : "FIXED";
    const char* G_status = st.fit_G ? "FLOAT" : "FIXED";

    out << "      {\n"
        << "        \"index\": " << i << ",\n"
        << "        \"shape\": \"" << (st.is_bw ? "Breit-Wigner" : "Gaussian") << "\",\n";

    WriteFitParamJSON(out, "yield", st.fac, fit.GetParameter(i_fac), fit.GetParError(i_fac), "FLOAT", 8);
    out << ",\n";
    WriteFitParamJSON(out, "centroid_MeV", st.M, fit.GetParameter(i_M), fit.GetParError(i_M), M_status, 8);
    out << ",\n";
    WriteFitParamJSON(out, "width_MeV", st.G, fit.GetParameter(i_G), fit.GetParError(i_G), G_status, 8);
    out << "\n      }" << (i + 1 < static_cast<int>(cfg.states.size()) ? "," : "") << "\n";
  }

  out << "    ],\n"
      << "    \"interference\": {\n"
      << "      \"enabled\": " << BoolText(cfg.use_interference) << ",\n"
      << "      \"status\": \"" << (cfg.fit_phase ? "FLOAT" : "FIXED") << "\",\n"
      << "      \"phase_set_rad\": " << phase_set_rad << ",\n"
      << "      \"phase_set_pi_over_2_units\": " << phase_set_pi2 << ",\n"
      << "      \"phase_fit_rad\": " << phase_fit_rad << ",\n"
      << "      \"phase_err_rad\": " << phase_err_rad << ",\n"
      << "      \"phase_fit_pi_over_2_units\": " << phase_fit_pi2 << ",\n"
      << "      \"phase_err_pi_over_2_units\": " << phase_err_pi2 << "\n"
      << "    },\n"
      << "    \"background\": {\n"
      << "      \"type\": \"" << ToString(cfg.background_type) << "\",\n";

  WriteFitParamJSON(out, "E0_MeV", cfg.bg_e.value, fit.GetParameter(cfg.BgEIndex()), fit.GetParError(cfg.BgEIndex()), "FLOAT", 6);
  out << ",\n";
  WriteFitParamJSON(out, "A0", cfg.bg_a0.value, fit.GetParameter(cfg.BgA0Index()), fit.GetParError(cfg.BgA0Index()), "FLOAT", 6);
  out << ",\n";
  WriteFitParamJSON(out, "A1", cfg.bg_a1.value, fit.GetParameter(cfg.BgA1Index()), fit.GetParError(cfg.BgA1Index()), "FLOAT", 6);
  if (cfg.background_type == BackgroundType::Quadratic) {
    out << ",\n";
    WriteFitParamJSON(out, "A2", cfg.bg_a2.value, fit.GetParameter(cfg.BgA2Index()), fit.GetParError(cfg.BgA2Index()), "FLOAT", 6);
  }

  out << "\n    }\n"
      << "  }\n"
      << "}\n";

  std::cout << "\nWrote JSON fit summary: " << outpath << "\n";
}

} // namespace bwfit
