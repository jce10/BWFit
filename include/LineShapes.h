#pragma once

namespace bwfit {

double Gaussian(double E, double M, double FWHM, double scale);
double BreitWigner(double E, double M, double Gamma, double scale);
double BreitWignerRel(double E, double M, double Gamma, double scale);
double KFactor(double M, double Gamma);

double GaussianRoot(double* x, double* p);
double BreitWignerRoot(double* x, double* p);
double BreitWignerRelRoot(double* x, double* p);

double LinearBackground(double E, double E0, double A0, double A1);
double QuadraticBackground(double E, double E0, double A0, double A1, double A2);

} // namespace bwfit
