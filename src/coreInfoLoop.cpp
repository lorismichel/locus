#include "utils.h"

/* Using Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying
 * such matrices would imply a prohibitive RAM overconsumption.
 */

// for locus_info_core function
// [[Rcpp::export]]
void coreInfoLoop(const MapMat X,
                  const MapMat Y,
                  MapArr2D gam_vb,
                  const MapArr2D log_Phi_mat_v_mu,
                  const MapArr2D log_1_min_Phi_mat_v_mu,
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

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) + c.transpose()));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_z_info_core and locus_mix_info_core function
// [[Rcpp::export]]
void coreZInfoLoop(const MapMat X,
                   const MapMat Y,
                   MapArr2D gam_vb,
                   const MapArr2D log_Phi_mat_v_mu,
                   const MapArr2D log_1_min_Phi_mat_v_mu,
                   const double log_sig2_inv_vb,
                   const MapArr1D log_tau_vb,
                   MapMat m1_beta,
                   MapMat mat_x_m1,
                   const MapMat mat_z_mu,
                   MapArr2D mu_beta_vb,
                   const MapArr1D sig2_beta_vb,
                   const MapArr1D tau_vb) {

  const Arr1D c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) + c.transpose()));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_logit_info_core function
// [[Rcpp::export]]
void coreLogitInfoLoop(const MapMat X,
                       const MapArr2D Y,
                       MapArr2D gam_vb,
                       const MapArr2D log_Phi_mat_v_mu,
                       const MapArr2D log_1_min_Phi_mat_v_mu,
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

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      log_sig2_inv_vb / 2 - mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.row(j)) -
      log(sig2_beta_vb.row(j)) / 2));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.matrix().noalias() += X.col(j) * m1_beta.row(j);

  }

}



// for locus_probit_info_core function
// [[Rcpp::export]]
void coreProbitInfoLoop(const MapMat X,
                        const MapMat Wy,
                        MapArr2D gam_vb,
                        const MapArr2D log_Phi_mat_v_mu,
                        const MapArr2D log_1_min_Phi_mat_v_mu,
                        const double log_sig2_inv_vb,
                        MapMat m1_beta,
                        MapMat mat_x_m1,
                        MapMat mat_z_mu,
                        MapArr2D mu_beta_vb,
                        const double sig2_beta_vb) {

  const double c = -(log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * ((Wy - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb) + c));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}

