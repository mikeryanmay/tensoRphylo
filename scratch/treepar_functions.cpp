
#include <Rcpp.h>
using namespace Rcpp;

// Stadler 2010 3.1
// [[Rcpp::export]]
double c1C(double l, double m, double psi) {
  return sqrt(pow(l - m - psi, 2) + 4 * l * psi);
}
// [[Rcpp::export]]
double c2C(double l, double m, double psi, double rho = 0) {
  return -(l - m - psi - 2 * l * rho) / c1C(l, m, psi);
}

// [[Rcpp::export]]
double phat1(double t, double l, double m, double psi, double rho = 1, bool logprob = false) {
  double co1 = c1C(l, m, psi);
  double co2 = c2C(l, m, psi, rho);
  double res = 4 * rho / pow(exp(co1 * t) * (1 + co2) + (1 - co2), 2);
  if (logprob) res = log(res);
  return res;
}

// probability that an individual alive at time t before today has no sampled extinct or extant descendants.
// [[Rcpp::export]]
double p0sersampC(double t, double l, double m, double psi, double rho = 0, bool logprob = false) {
  double co1 = c1C(l, m, psi);
  double co2 = c2C(l, m, psi, rho);
  double tmp = exp(-co1 * t) * (1 - co2);
  double res = (l + m + psi + co1 * (tmp - (1 + co2)) / (tmp + (1 + co2))) / (2 * l);
  if (logprob) res = log(res);
  return res;
}

// rho / p1(t), where p1(t) is the probability that an individual alive at time t before today has precisely one sampled extant descendant and no sampled extinct descendant.
// [[Rcpp::export]]
double qfuncC(double t, double l, double m, double psi, double rho = 0, bool logprob = false) {
  double co1 = c1C(l, m, psi);
  double co2 = c2C(l, m, psi, rho);
  double tmp = ((1 - co2) * exp(-t * co1 / 2) + (1 + co2) * exp(t * co1 / 2));
  double res = logprob ? (2 * log(tmp) - log(4)) : (pow(tmp, 2) / 4);
  return res;
}
