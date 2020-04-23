#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double chisq_fix_reject(S4 design, double n1, double nuisance, String type) {
  double reject_prob = 0;
  double r = design.slot("r");
  double alpha = design.slot("alpha");
  double delta = design.slot("delta");
  String alternative = design.slot("alternative");
  double krit = R::qnorm(1 - alpha, 0, 1, TRUE, FALSE);
  double n_c = ceil(n1 / (r + 1));
  double n_e = ceil(n1 * r / (r + 1));
  double r_act = n_e / n_c;
  bool ind = false;
  double p1_e = 0;
  double p1_c = 0;

  if (type == "power") {
    p1_e = nuisance + delta / (1 + r_act);
    p1_c = p1_e - delta;
  }

  for(int i = 0; i <= n_c; ++i) {
    for(int j = 0; j <= n_e; ++j) {
      if (((i + j) == 0) | ((i + j) == (n_c + n_e))) {
        continue;
      }
      double p_c = i / n_c;
      double p_e = j / n_e;
      double p_hat = (i + j) / (n_c + n_e);
      double ts = sqrt(n_c * n_e / (n_c + n_e)) * (p_e - p_c) /
        sqrt(p_hat * (1 - p_hat));

      if (alternative == "greater") {
        ind = ts > krit;
      } else {
        ind = ts < -krit;
      }

      if (ind & (type == "size")) {
        double sum_prob = Rf_choose(n_c, i) * Rf_choose(n_e, j) *
          pow(nuisance, i + j) * pow(1 - nuisance, n_c + n_e - i - j);
        reject_prob += sum_prob;
      } else if (ind & (type == "power")) {
        double sum_prob = Rf_choose(n_c, i) * Rf_choose(n_e, j) * pow(p1_c, i) *
          pow(1 - p1_c, n_c - i) * pow(p1_e, j) * pow(1 - p1_e, n_e - j);
        reject_prob += sum_prob;
      }
    }
  }
  return reject_prob;
}

// [[Rcpp::export]]
double chisq_recalc_reject(S4 design, double n1, double nuisance,
                           String type, NumericMatrix nmat) {
  bool ind;
  double p1_e, p1_c, n2, p_c, p_e, p_hat, n_new, ts, x;
  double sum_prob, n_c2, n_e2, n_cdot, n_edot, n_tot, x_c, x_e;

  double reject_prob = 0;
  double r = design.slot("r");
  double alpha = design.slot("alpha");
  double delta = design.slot("delta");
  String alternative = design.slot("alternative");
  double krit = R::qnorm(1 - alpha, 0, 1, TRUE, FALSE);
  double n_c1 = ceil(n1 / (r + 1));
  double n_e1 = ceil(n1 * r / (r + 1));
  double r_act = n_e1 / n_c1;

  if (type == "power") {
    p1_e = nuisance + delta / (1 + r_act);
    p1_c = p1_e - delta;
  }

  for (int i = 0; i <= n_c1; ++i) {
    for (int j = 0; j <= n_e1; ++j) {
      Rcpp::checkUserInterrupt();
      p_c = i / n_c1;
      p_e = j / n_e1;
      p_hat = (i + j) / (n_c1 + n_e1);
      n_new = nmat((i + j * (n_c1 + 1)), 3);

      if (((n_new > 0) & (n_new <= n1)) | (n_new == -99)) {
        if ((p_hat == 0) | (p_hat == 1)) {
          continue;
        }
        ts = sqrt(n_c1 * n_e1 / (n_c1 + n_e1)) * (p_e - p_c) /
          sqrt(p_hat * (1 - p_hat));

        if (alternative == "greater") {
          ind = ts > krit;
        } else {
          ind = ts < -krit;
        }

        if (ind & (type == "size")) {
          sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
            pow(nuisance, i + j) * pow(1 - nuisance, n_e1 + n_c1 - i - j);
          reject_prob += sum_prob;
        } else if (ind & (type == "power")) {
          sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
            pow(p1_c, i) * pow(1 - p1_c, n_c1 - i) * pow(p1_e, j) *
            pow(1 - p1_e, n_e1 - j);
          reject_prob += sum_prob;
        }
      } else {
        n2 = n_new - (n_c1 + n_e1);
        n_c2 = ceil(n2 / (r + 1));
        n_e2 = ceil(n2 * r / (r + 1));

        for (int k = 0; k <= n_c2; ++k) {
          for (int l = 0; l <= n_e2; ++l) {
            n_cdot = n_c1 + n_c2;
            n_edot = n_e1 + n_e2;
            p_c = (i + k) / n_cdot;
            p_e = (j + l) / n_edot;
            p_hat = (i + j + k + l) / (n_cdot + n_edot);
            ts = sqrt(n_cdot * n_edot / (n_cdot + n_edot)) * (p_e - p_c) /
              sqrt(p_hat * (1 - p_hat));

            if (alternative == "greater") {
              ind = ts > krit;
            } else {
              ind = ts < -krit;
            }

            if (ind & (type == "size")) {
              x = i + j + k + l;
              n_tot = n_c1 + n_e1 + n_c2 + n_e2;
              sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
                Rf_choose(n_c2, k) * Rf_choose(n_e2, l) * pow(nuisance, x) *
                pow(1 - nuisance, n_tot - x);
              reject_prob += sum_prob;
            } else if (ind & (type == "power")) {
              x_c = i + k;
              x_e = j + l;
              sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
                Rf_choose(n_c2, k) * Rf_choose(n_e2, l) * pow(p1_c, x_c) *
                pow(1 - p1_c, n_cdot - x_c) * pow(p1_e, x_e) *
                pow(1 - p1_e, n_edot - x_e);
              reject_prob += sum_prob;
            }
          }
        }
      }
    }
  }
  return reject_prob;
}
