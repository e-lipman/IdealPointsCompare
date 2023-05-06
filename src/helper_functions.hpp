#ifndef helper_functions
#define helper_functions

#include "parameters.hpp"

void print_datetime();

double rpg(double shape, double scale);

template<class N, class M>
void prebind_col_to_matrix(N &col, M &mat, M &out){
  if ( (out.nrow()!=mat.nrow()) | (out.ncol()!=mat.ncol()+1) ){
    throw std::invalid_argument("Dimension error in 'prebind_col_to_matrix'");
  }
  out(_,0) = col;
  for (int j=0; j<mat.ncol(); j++){
    out(_, j+1) = mat(_, j);
  }
}

double compute_lpost_zeta_hier(int I, int p,
                               Rcpp::LogicalVector &zeta,
                               arma::uvec &eps,
                               arma::mat &x,
                               arma::vec &eta,
                               double intercept,
                               Rcpp::NumericVector &theta,
                               double int_prior_scale);
  
double compute_llik(int I, int J, int gam0_max,
                    const Rcpp::LogicalMatrix &y,
                    const Rcpp::LogicalMatrix &missing,
                    
                    Rcpp::NumericVector &mu,
                    Rcpp::NumericVector &alf,
                    Rcpp::NumericMatrix &beta);

double compute_lpost(int I, int J, int gam0_max,
                     const Rcpp::LogicalMatrix &y,
                     const Rcpp::LogicalMatrix &missing,
                     
                     Rcpp::NumericVector &mu,
                     Rcpp::NumericVector &alf,
                     Rcpp::NumericMatrix &beta,
                     Rcpp::LogicalVector &zeta,
                     Rcpp::NumericVector &theta,
                     
                     double rho_mu, double kap2_mu,
                     double kap2_alf, double w_alf,
                     double eta, double sig2,
                     
                     ParametersHier *pp_x,
                     double int_prior_scale);

double get_scale_boonstra(int int_prior_fam,
                          int I, double q);
#endif