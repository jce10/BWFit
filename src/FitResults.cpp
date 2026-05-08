#include "FitResults.h"


namespace bwfit {
void PrintFitParameters(const TF1& fit, const FitConfig& cfg) {

  std::cout << "\n=== Fit Parameters ===\n";

  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    const int i_fac = 3 * i + 0;
    const int i_M   = 3 * i + 1;
    const int i_G   = 3 * i + 2;

    std::cout << "State " << i << ":\n";
    std::cout << "  Yield = " << fit.GetParameter(i_fac)
              << " +/- " << fit.GetParError(i_fac) << "\n";
    std::cout << "  M     = " << fit.GetParameter(i_M)
              << " +/- " << fit.GetParError(i_M) << " MeV\n";
    std::cout << "  G     = " << fit.GetParameter(i_G)
              << " +/- " << fit.GetParError(i_G) << " MeV\n";
  }

  std::cout << "\nInterference:\n";
  std::cout << "  phase = " << fit.GetParameter(cfg.PhaseIndex())
            << " +/- " << fit.GetParError(cfg.PhaseIndex()) << "\n";

  std::cout << "\nBackground:\n";
  std::cout << "  E0 = " << fit.GetParameter(cfg.BgEIndex())
            << " +/- " << fit.GetParError(cfg.BgEIndex()) << "\n";
  std::cout << "  A0 = " << fit.GetParameter(cfg.BgA0Index())
            << " +/- " << fit.GetParError(cfg.BgA0Index()) << "\n";
  std::cout << "  A1 = " << fit.GetParameter(cfg.BgA1Index())
            << " +/- " << fit.GetParError(cfg.BgA1Index()) << "\n";

  if (cfg.background_type == BackgroundType::Quadratic) {
    std::cout << "  A2 = " << fit.GetParameter(cfg.BgA2Index())
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