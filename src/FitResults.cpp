#include "FitResults.h"
#include <iomanip>
#include <iostream>


namespace bwfit {
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

  std::cout << "\nInterference phase, delta:\n";
  std::cout << "  set = " << phase_set_rad
            << " rad = " << phase_set_pi2 << " * pi/2\n";

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
} // namespace bwfit