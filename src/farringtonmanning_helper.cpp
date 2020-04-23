#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector p_rml(double p_c, double p_e, double r, double margin) {
  double x;
  double pi = std::atan(1)*4;
  NumericVector pt(2);

  double a = 1 + (1 / r);
  double b = -(1 + (1 / r) + p_e + (1 / r) * p_c - margin * ((1 / r) + 2));
  double c = pow(margin, 2) - margin * (2 * p_e + (1 / r) + 1) + p_e + (1 / r) * p_c;
  double d = p_e * margin * (1 - margin);

  double v = pow(b, 3) / pow(3 * a, 3) - (b * c) / (6 * pow(a, 2)) + d / (2 * a);
  double u = ((v > 0) - (v < 0)) * sqrt(pow(b, 2) / pow(3 * a, 2) - c / (3 * a));

  if ((v == 0) & (u == 0)) {
    x = 0;
  } else if ((v / pow(u, 3)) > 1) {
    x = 1;
  } else {
    x = v / pow(u, 3);
  }

  double w = (pi + acos(x)) / 3;

  pt[1] = fmax(0, 2 * u * cos(w) - b / (3 * a));
  pt[0] = fmin(1, pt[1] + margin);

  return pt;
}

// [[Rcpp::export]]
double fm_fix_reject(S4 design, double n, double nuisance, String type) {
  double delta_NI = design.slot("delta_NI");
  double delta = design.slot("delta");
  double r = design.slot("r");
  double alpha = design.slot("alpha");
  double p1_e, p1_c, p_c, p_e, pt_c, pt_e, se, ts, sum_prob;
  NumericVector pt;
  bool ind;

  double reject_prob = 0;
  double krit = R::qnorm(1 - alpha, 0, 1, TRUE, FALSE);
  double n_c = ceil(n / (r + 1));
  double n_e = ceil(n * r / (r + 1));
  double r_act = n_e / n_c;
  double p0_e = nuisance - delta_NI / (1 + r_act);
  double p0_c = p0_e + delta_NI;



  if (type == "power") {
    p1_e = nuisance - delta / (1 + r_act);
    p1_c = p1_e + delta;
  }

  for (int i = 0; i <= n_c; ++i) {
    for (int j = 0; j <= n_e; ++j) {
      Rcpp::checkUserInterrupt();
      p_c = i / n_c;
      p_e = j / n_e;
      pt = p_rml(p_c, p_e, r_act, delta_NI);
      pt_c = pt[0];
      pt_e = pt[1];

      se = sqrt((r_act * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_e);
      ts = (p_e - p_c + delta_NI) / se;
      ind = ts > krit;

      if (ind & (type == "size")) {
        sum_prob = Rf_choose(n_c, i) * (Rf_choose(n_e, j) * pow(p0_c, i)) *
          pow(1 - p0_c, n_c - i) * pow(p0_e, j) * pow(1 - p0_e, n_e - j);
        reject_prob += sum_prob;
      } else if (ind & (type == "power")) {
        sum_prob = Rf_choose(n_c, i) * (Rf_choose(n_e, j) * pow(p1_c, i)) *
          pow(1 - p1_c, n_c - i) * pow(p1_e, j) * pow(1 - p1_e, n_e - j);
        reject_prob += sum_prob;
      }
    }
  }
  return reject_prob;
}

// [[Rcpp::export]]
double fm_recalc_reject(S4 design, double n1, double nuisance,
                        String type, NumericMatrix nmat) {
  double delta_NI = design.slot("delta_NI");
  double delta = design.slot("delta");
  double r = design.slot("r");
  double alpha = design.slot("alpha");
  double p1_e, p1_c, p_c, p_e, pt_c, pt_e, se, ts, sum_prob;
  double n_new, n2, n_c2, n_e2, n_cdot, n_edot, r_act2;
  NumericVector pt;
  bool ind;

  double reject_prob = 0;
  double krit = R::qnorm(1 - alpha, 0, 1, TRUE, FALSE);
  double p0_e = nuisance - delta_NI / (1 + r);
  double p0_c = p0_e + delta_NI;
  double n_c1 = ceil(n1 / (r + 1));
  double n_e1 = ceil(n1 * r / (r + 1));
  double r_act1 = n_e1 / n_c1;


  if (type == "power") {
    p1_e = nuisance - delta / (1 + r_act1);
    p1_c = p1_e + delta;
  }

  for (int i = 0; i <= n_c1; ++i) {
    for (int j = 0; j <= n_e1; ++j) {
      Rcpp::checkUserInterrupt();
      p_c = i / n_c1;
      p_e = j / n_e1;
      n_new = nmat((i + j * (n_c1 + 1)), 3);

      if (((n_new > 0) & (n_new <= (n_c1 + n_e1))) | (n_new == -99)) {
        pt = p_rml(p_c, p_e, r_act1, delta_NI);
        pt_c = pt[0];
        pt_e = pt[1];

        se = sqrt((r_act1 * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_e1);
        ts = (p_e - p_c + delta_NI) / se;
        ind = ts > krit;

        if (ind & (type == "size")) {
          sum_prob = Rf_choose(n_c1, i) * (Rf_choose(n_e1, j) *
            pow(p0_c, i)) * pow(1 - p0_c, n_c1 - i) * pow(p0_e, j) *
            pow(1 - p0_e, n_e1 - j);
          reject_prob += sum_prob;
        } else if (ind & (type == "power")) {
          sum_prob = Rf_choose(n_c1, i) * (Rf_choose(n_e1, j) *
            pow(p1_c, i)) * pow(1 - p1_c, n_c1 - i) * pow(p1_e, j) *
            pow(1 - p1_e, n_e1 - j);
          reject_prob += sum_prob;
        }
      } else {
        n2 = n_new - (n_c1 + n_e1);
        n_c2 = ceil(n2 / (r + 1));
        n_e2 = ceil(n2 * r / (r + 1));
        r_act2 = (n_e1 + n_e2) / (n_c1 + n_c2);

        for (int k = 0; k <= n_c2; ++k) {
          for (int l = 0; l <= n_e2; ++l) {
            n_cdot = n_c1 + n_c2;
            n_edot = n_e1 + n_e2;
            p_c = (i + k) / n_cdot;
            p_e = (j + l) / n_edot;
            pt = p_rml(p_c, p_e, r_act2, delta_NI);
            pt_c = pt[0];
            pt_e = pt[1];
            se = sqrt((r_act2 * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_edot);
            ts = (p_e - p_c + delta_NI) / se;
            ind = ts > krit;

            if (ind & (type == "size")) {
              sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
                Rf_choose(n_c2, k) * Rf_choose(n_e2, l) *
                pow(p0_c, i + k) * pow(1 - p0_c, n_cdot - i - k) *
                pow(p0_e, j + l) * pow(1 - p0_e, n_edot - j - l);
              reject_prob += sum_prob;
            } else if (ind & (type == "power")) {
              sum_prob = Rf_choose(n_c1, i) * Rf_choose(n_e1, j) *
                Rf_choose(n_c2, k) * Rf_choose(n_e2, l) *
                pow(p1_c, i + k) * pow(1 - p1_c, n_cdot - i - k) *
                pow(p1_e, j + l) * pow(1 - p1_e, n_edot - j - l);
              reject_prob += sum_prob;
            }
          }
        }
      }
    }
  }
  return reject_prob;
}
