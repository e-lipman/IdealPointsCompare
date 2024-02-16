#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "configs.hpp"
#include "parameters.hpp"
#include "helloPG/RNG.h" // from https://github.com/jgscott/helloPG
#include "helloPG/PolyaGamma.h"

static double const log2pi = std::log(2.0 * M_PI);

using namespace Rcpp;

void print_datetime(){
  time_t now;
  time(&now);
  
  struct tm *local = localtime(&now);
  int hours = local->tm_hour;         // get hours since midnight (0-23)
  int minutes = local->tm_min;        // get minutes passed after the hour (0-59)
  int seconds = local->tm_sec;        // get seconds passed after a minute (0-59)
  int day = local->tm_mday;            // get day of month (1 to 31)
  int month = local->tm_mon + 1;      // get month of year (0 to 11)
  int year = local->tm_year + 1900;   // get year since 1900
  
  Rprintf("%04d-%02d-%02d %02d:%02d:%02d: ", 
          year, month, day, hours, minutes, seconds);
}

double rpg(double shape, double scale) {
  RNG r;
  PolyaGamma pg;
#ifdef USE_R
  GetRNGstate();
#endif
  double result = pg.draw(shape, scale, r);
#ifdef USE_R
  PutRNGstate();
#endif
  return result;
}

/* helper functions for unnormalized posterior */
// https://gallery.rcpp.org/articles/dmvnorm_arma/

static void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

static arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

/* compute unnormalized posterior */

double compute_lpost_zeta_hier(int I, int p,
                                Rcpp::LogicalVector &zeta,
                                arma::uvec &eps,
                                arma::mat &x,
                                arma::vec &eta,
                                double intercept,
                                Rcpp::NumericVector &theta,
                                double int_prior_scale){
  
  double zeta_prior=0, eps_prior=0, eta_prior=0, int_prior=0;
  
  // zeta ~  logit(x[,eps]*eta)
  for (int i=0; i<I; i++){
    zeta_prior+=R::dbinom(zeta[i], 1, theta[i], true);
  }
  
  // eps ~ 1/[(p+1)]*(p choose sum(eps))]
  eps_prior -= log(Rf_choose(p,sum(eps)))-log(p+1); // total covs, selected covs
  
  // eta | eps ~ N(0, 4I(xtx)^{-1})
  arma::mat x1 = x.cols(arma::find(eps == 1));
  arma::mat Sig = 4*I*inv(x1.t()*x1);
  eta_prior += dmvnrm_arma_fast(x, arma::rowvec(p).fill(0), Sig, true)[0];
  
  // intercept prior
  if (int_prior_fam==0){
    int_prior = R::dnorm(intercept, 0, int_prior_scale, true);
  } else {
    int_prior = R::dlogis(intercept, 0, int_prior_scale, true);
  }
  
  return(zeta_prior + eps_prior + eta_prior+int_prior);
}

double compute_llik(int I, int J, int gam0_max,
                    const Rcpp::LogicalMatrix &y,
                    const Rcpp::LogicalMatrix &missing,
                    
                    Rcpp::NumericVector &mu,
                    Rcpp::NumericVector &alf,
                    Rcpp::NumericMatrix &beta){
  // y_ij ~ logit(mu_j + alf_j*beta_ig)
  double llik = 0;
  for (int i=0; i<I; i++){
    for (int j=0; j<J; j++){
      if (missing(i,j)==0){
        int gam = j<=gam0_max ? 0 : 1;
        double psi = mu[j] + alf[j]*beta(i,gam);
        llik += R::dbinom(y(i,j), 1, exp(psi)/(1+exp(psi)), true);
      }
    }
  }
  return(llik);
}

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
                     double int_prior_scale){
  
  double llik=0, lprior_mu=0, lprior_alf=0, lprior_beta=0, lprior_zeta=0;
  
  /* likelihood */
  llik = compute_llik(I, J, gam0_max, y, missing, mu, alf, beta);

  /* mu prior */
  // mu_j ~ N(rho_mu, kap2_mu)
  for (int j=0; j<J; j++){
    lprior_mu += R::dnorm(mu[j], rho_mu, sqrt(kap2_mu), true);
  }
  lprior_mu += R::dnorm(rho_mu, rho_mu_mean, sqrt(rho_mu_var), true);
  lprior_mu += R::dgamma(1/kap2_mu, kap2_mu_a, 1/kap2_mu_b, true);

  /* alf prior */
  // alf_j = w_alf*I(0) + (1-w_alf)*N(0, kap2_alf)
  int count_alf_nonzero=0;
  for (int j=0; j<J; j++){
    if (alf[j]!=0){
      count_alf_nonzero++;
      lprior_alf += R::dnorm(alf[j], 0, sqrt(kap2_alf), true);
    }
  }
  lprior_alf += R::dbinom(J-count_alf_nonzero, J, w_alf, true);
  lprior_alf += R::dgamma(1/kap2_alf, kap2_alf_a, 1/kap2_alf_b, true);
  lprior_alf += R::dbeta(w_alf, w_a, w_b, true);

  /* beta prior */
  // zeta_i=0: beta_i0=beta_i1 ~ N(eta, sig2) 
  // zeta_i=1: beta_i0, beta_i1 ~ N(eta, sig2) 
  for (int i=0; i<I; i++){
    lprior_beta += R::dnorm(beta(i,0), eta, sqrt(sig2), true);
    if (zeta[i]==0){
      lprior_beta += R::dnorm(beta(i,1), eta, sqrt(sig2), true);
    }
  }
  lprior_beta += R::dnorm(eta, eta_mean, sqrt(eta_var), true);
  lprior_beta += R::dgamma(1/sig2, sig2_a, 1/sig2_b, true);
  
  /* zeta prior */
  if (pp_x == NULL){ // non-hierarchical model
    int count_zeta_one=0;
    for (int i=0; i<I; i++){count_zeta_one += zeta[i];}
    lprior_zeta += R::dbinom(count_zeta_one, I, theta[0], true);
    lprior_zeta += R::dbeta(theta[0], theta_a, theta_b, true);
  } else {
    lprior_zeta += compute_lpost_zeta_hier(I, pp_x->p_, zeta,  
                                           pp_x->eps, pp_x->x_, pp_x->eta, 
                                           pp_x->intercept, theta, int_prior_scale);
  }
  
  return(llik + lprior_mu + lprior_mu + lprior_beta + lprior_zeta);
}

// intercept prior scale Boonstra
double get_scale_boonstra(int int_prior_fam,
                          int I, double q){
  double I2 = 2*I;
  double s = exp(-1/I2);
  double s2 = log(s) - log(1-s);
  double scale=0;
  if (int_prior_fam==0){
    scale = s2/R::qnorm(1-q/2,0,1,true,false);
  } else if (int_prior_fam==1){
    scale = s2/R::qlogis(1-q/2,0,1,true,false);
  }
  return(scale);
}
