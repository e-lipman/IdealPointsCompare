#ifndef parameters
#define parameters

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "configs.hpp"
#include "gibbs_updates.hpp"
using namespace Rcpp;

class Parameters
{
public: 
  Parameters(int I, int J, int gam0_max,
             LogicalMatrix y_raw, LogicalMatrix missing): 
  
  I_(I), J_(J), gam0_max_(gam0_max), missing_(I,J),
  
  // hyperparams
  rho_mu(0), kap2_mu(1), kap2_alf(1), w_alf(.5), 
  rho_beta(0), kap2_beta(1), theta (I),
  
  // params
  mu (J), alf(J), beta(I,2), zeta(I,1),
  
  // augmented data
  y(I,J), v(I,J), z(I,J)
  
  // sample initial values
  {
    for (int j=0; j<J; j++){ // sample from prior
      mu[j] = R::rnorm(rho_mu, sqrt(kap2_mu));
      alf[j] = R::rbinom(1, 1-w_alf) * R::rnorm(0, sqrt(kap2_alf));
    }
    // sample from full conditionals
    y = y_raw;
    missing_ = missing;
    update_yvz(I, J, gam0_max, y, v, z, missing, mu, alf, beta);
  }
  
  Parameters(int I):
  I_(I),
  
  J_(0), gam0_max_(0),  
  rho_mu(0), kap2_mu(0), kap2_alf(0), w_alf(0),  
  rho_beta(0), kap2_beta(0),
  
  theta(I)
  { 
  }
  
  virtual ~Parameters()=default;
  void call_update_zeta_beta(){
    update_zeta_beta(I_, J_, gam0_max_, rho_beta, kap2_beta, false,
                     theta, v, z, mu, alf, beta, zeta);
  }
  
  int I_;
  int J_;
  int gam0_max_;
  LogicalMatrix missing_;
  
  double rho_mu;
  double kap2_mu;  
  double kap2_alf;
  double w_alf; 
  double rho_beta;
  double kap2_beta;
  NumericVector theta;
  
  NumericVector mu;
  NumericVector alf;
  NumericMatrix beta;
  LogicalVector zeta;
  
  LogicalMatrix y;
  NumericMatrix v;
  NumericMatrix z;
};

class ParametersNohier:public Parameters
{
public:
  ParametersNohier(int I, int J, int gam0_max,
                   LogicalMatrix y_raw, LogicalMatrix missing):
  Parameters(I, J, gam0_max, y_raw, missing)
  {
    double theta_val = R::rbeta(1, 1);
    theta = NumericVector (I_, theta_val);
    call_update_zeta_beta();
  }
  virtual ~ParametersNohier()=default;
};


class ParametersHier:public Parameters
{
public:
  ParametersHier(int I, int J, arma::mat &x, int gam0_max, int p,
                 LogicalMatrix y_raw, LogicalMatrix missing,
                 int eps_init_method):
  Parameters(I, J, gam0_max, y_raw, missing),
  p_(p), x_(x), eps(p), vx(I), zx(I), eta(p), intercept(0), linpred(I),
  eps_init(eps_init_method)
  {
    initialize_params();
    call_update_zeta_beta();
  }
  ParametersHier(int I,int p, arma::mat x, 
                 int eps_init_method):
  Parameters(I), p_(p), x_(x), eps(p), vx(I), zx(I), 
  eta(p), intercept(0), linpred(I), 
  accept_count(0), accept_count2(0), accept_count3(0),
  eps_init(eps_init_method), g(0)
  {
    initialize_params();
  }
  
  virtual~ParametersHier()=default;
  
  void initialize_params(){
    g = I_;
    if (eps_init == 0){ // null model
      eps.fill(0);
    } else if (eps_init == 1){ // full model
      eps.fill(1);
    } else if (eps_init == 2){ // random model
      for (int k=0; k<p_; k++){eps[k] = R::rbinom(1,0.5);}
    } else {
      throw std::invalid_argument("'eps_init_method' must be in 0,1,2");
    }
    //Rcout << eps.as_row() << std::endl;
    
    eta.fill(0);
    intercept=R::rnorm(0,1);
    for (int k=0; k<p_; k++){
      if (eps[k]==1){
        eta[k] = R::rnorm(0,1);
      }
    }
    linpred = x_*eta + intercept;
    for (int i=0; i<I_; i++){theta[i] = exp(linpred[i])/(1+exp(linpred[i]));}
  }
  
  int p_;
  arma::mat x_;
  arma::uvec eps;
  arma::vec vx;
  arma::vec zx;
  arma::vec eta;
  double intercept;
  arma::vec linpred;
  int accept_count;
  int accept_count2;
  int accept_count3;
  int eps_init;
  double g;
};

#endif