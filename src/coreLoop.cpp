#include "utils.h"

/* Using Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying
 * such matrices would imply a prohibitive RAM overconsumption.
 */

// for locus_core function
// [[Rcpp::export]]
void coreLoop(const MapMat X,
              const MapMat Y,
              MapArr2D gam_vb,
              const MapArr1D log_om_vb,
              const MapArr1D log_1_min_om_vb,
              const double log_sig2_inv_vb,
              const MapArr1D log_tau_vb,
              MapMat m1_beta,
              MapMat mat_x_m1,
              MapArr2D mu_beta_vb,
              const MapArr1D sig2_beta_vb,
              const MapArr1D tau_vb) {

  const Arr1D c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_om_vb(j) - log_om_vb(j) -
      mu_beta_vb.row(j).transpose().square() / (2 * sig2_beta_vb) + c));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_z_core and locus_mix_core function
// [[Rcpp::export]]
void coreZLoop(const MapMat X,
               const MapMat Y,
               MapArr2D gam_vb,
               const MapArr1D log_om_vb,
               const MapArr1D log_1_min_om_vb,
               const double log_sig2_inv_vb,
               const MapArr1D log_tau_vb,
               MapMat m1_beta,
               MapMat mat_x_m1,
               MapMat mat_z_mu,
               MapArr2D mu_beta_vb,
               const MapArr1D sig2_beta_vb,
               const MapArr1D tau_vb) {

  const Arr1D c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_om_vb(j) - log_om_vb(j) -
      mu_beta_vb.row(j).transpose().square() / (2 * sig2_beta_vb) + c));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_logit_core function
// [[Rcpp::export]]
void coreLogitLoop(const MapMat X,
                   const MapArr2D Y,
                   MapArr2D gam_vb,
                   const MapArr1D log_om_vb,
                   const MapArr1D log_1_min_om_vb,
                   const double log_sig2_inv_vb,
                   MapMat m1_beta,
                   MapArr2D mat_x_m1,
                   MapArr2D mat_z_mu,
                   MapArr2D mu_beta_vb,
                   const MapArr2D psi_vb,
                   const MapArr2D sig2_beta_vb) {

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.matrix().noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb.row(j) * (X.col(j).transpose() * (Y - 2 * psi_vb * (mat_x_m1 + mat_z_mu)).matrix()).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_om_vb(j) - log_om_vb(j) -
      log_sig2_inv_vb / 2 - mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.row(j)) -
      log(sig2_beta_vb.row(j)) / 2));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.matrix().noalias() += X.col(j) * m1_beta.row(j);

  }

}



// for locus_probit_core function
// [[Rcpp::export]]
void coreProbitLoop(const MapMat X,
                    const MapMat W,
                    MapArr2D gam_vb,
                    const MapArr1D log_om_vb,
                    const MapArr1D log_1_min_om_vb,
                    const double log_sig2_inv_vb,
                    MapMat m1_beta,
                    MapMat mat_x_m1,
                    MapMat mat_z_mu,
                    MapArr2D mu_beta_vb,
                    const double sig2_beta_vb) {

  const double c = -(log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * ((W - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_om_vb(j) - log_om_vb(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb) + c));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_horseshoe_core function
// [[Rcpp::export]]
void coreHorseShoeLoop(const MapMat X,
              const MapMat Y,
              MapMat mat_x_m1,
              MapMat m1_beta,
              MapArr2D mu_beta_vb,
              const MapArr1D sig2_beta_vb,
              const MapArr1D tau_vb) {

  //const Arr1D c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1).transpose() * X.col(j)).array();

    m1_beta.row(j) = mu_beta_vb.row(j);
    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }



}

