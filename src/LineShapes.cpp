#include "LineShapes.h"

#include <cmath>

namespace bwfit {

double Gaussian(double E, double M, double FWHM, double scale) {
  return scale * std::sqrt(4.0 * std::log(2.0) / M_PI) / FWHM *
         std::exp(-4.0 * std::log(2.0) * (E - M) * (E - M) / (FWHM * FWHM));
}

double BreitWigner(double E, double M, double Gamma, double scale) {
  return scale * (Gamma / (2.0 * M_PI)) /
         ((E - M) * (E - M) + Gamma * Gamma / 4.0);
}

double BreitWignerRel(double E, double M, double Gamma, double scale) {
  return scale * KFactor(M, Gamma) /
         ((E * E - M * M) * (E * E - M * M) + M * M * Gamma * Gamma);
}

double KFactor(double M, double Gamma) {
  double k = 2.0 * std::sqrt(2.0) * M * M * Gamma * std::sqrt(M * M + Gamma * Gamma);
  k /= M_PI * std::sqrt(M * M + M * std::sqrt(M * M + Gamma * Gamma));
  return k;
}

double GaussianRoot(double* x, double* p) {
  return Gaussian(x[0], p[1], p[2], p[0]);
}

double BreitWignerRoot(double* x, double* p) {
  return BreitWigner(x[0], p[1], p[2], p[0]);
}

double BreitWignerRelRoot(double* x, double* p) {
  return BreitWignerRel(x[0], p[1], p[2], p[0]);
}

double LinearBackground(double E, double E0, double A0, double A1) {
  const double x = E - E0;
  return A0 + A1 * x;
}

double QuadraticBackground(double E, double E0, double A0, double A1, double A2) {
  const double x = E - E0;
  return A0 + A1 * x + A2 * x * x;
}

} // namespace bwfit
