#include <string.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "helper_functions.hpp"

using namespace Rcpp;

static arma::uvec gen_eps_proposal(int p, arma::uvec eps){
  // eps should not include intercept
  
  const Rcpp::NumericVector p1 = {0.9, 0.1};
  const Rcpp::NumericVector p2 = {0.6, 0.2, 0.15, 0.05};
  
  arma::uvec eps_prop = eps;
  
  // generate set of flips
  int num_flips = 0;
  unsigned int move_type = Rcpp::sample(2, 1, false, p1)[0];
  
  if (move_type==1){
    // flip 1, 2, 3 or 4 variables
    num_flips = Rcpp::sample(p2.length(), 1, false, p2)[0];
    if (num_flips>p){num_flips=p;}
    Rcpp::IntegerVector flips = Rcpp::sample(p, num_flips, false)-1; //0:(p-1)
    
    for (int k=0; k<num_flips; k++){
      arma::uvec idx = arma::uvec(1).fill(flips[k]); 
      arma::uvec vals = arma::uvec(1).fill(1-eps_prop[idx[0]]); 
      eps_prop.elem(idx) = vals;
    }
  } else {
    // Add and remove one variable
    if (sum(eps)>0 && sum(eps)<eps.n_elem){
      Rcpp::IntegerVector idx_0 = as<IntegerVector>(wrap(find(eps == 0)));
      arma::uvec flip0 = arma::uvec(1).fill(Rcpp::sample(idx_0, 1, false)[0]);
      eps_prop.elem(flip0) = arma::uvec(1).fill(1);
    
      Rcpp::IntegerVector idx_1 = as<IntegerVector>(wrap(find(eps == 1)));
      arma::uvec flip1 = arma::uvec(1).fill(Rcpp::sample(idx_1, 1, false)[0]); 
      eps_prop.elem(flip1) = arma::uvec(1).fill(0);
    }
  }

  return(eps_prop);
}

static double get_lBF0(int I, int p,
                     arma::uvec eps, 
                     arma::mat &x, 
                     arma::mat &vx, 
                     arma::vec &zx, 
                     arma::vec &mu,
                     arma::mat &Sig,
                     double intercept,
                     double g){
  
  if (p==0){
    return(0);
  } else {
    double gg = 4*g; // scale for inv(Xtx)
    arma::mat vmat = arma::diagmat(vx);
    arma::mat x1 = x.cols(arma::find(eps == 1));
    
    arma::mat xx = x1.t() * x1;
    arma::mat xvx = x1.t() * vmat * x1;
    arma::mat xx_inv = arma::inv(xx);
    arma::mat xvgx_inv = inv(xvx + (1/gg)*xx);
    
    arma::vec za = zx - intercept;
    arma::vec zavx = x1.t() * vmat * za;
    
    // calculate Bayes factorial(agains null model)
    double lterm1 = 
      -0.5*real(arma::log_det(arma::eye(p,p) + gg*xvx*xx_inv));
    arma::vec lterm2 = 0.5*zavx.t() * xvgx_inv * zavx;
    double lBF = lterm1 + lterm2[0];
      
    // calculate conditional mean and var for eta under eps
    Sig = xvgx_inv;
      
    mu = Sig * x1.t()*vmat*za;
      
    return(lBF);  
  }
}

void update_eps_eta(int I, int p,
                    arma::uvec &eps,
                    arma::mat &x,
                    arma::vec &vx,
                    arma::vec &zx,
                    arma::vec &eta,
                    double intercept,
                    int *accept_count,
                    double *g,
                    int g_hyper){
  
  
  // proposal for eps (non-intercept covariates)
  arma::uvec eps_prop;
  eps_prop = gen_eps_proposal(p, eps);
  int p1 = sum(eps);
  int p2 = sum(eps_prop);
  
  // propose g
  double g_prop,thresh=0;
  if (g_hyper!=0){
    double thresh = (I-p2)/(p2+1);
    double prop_sd = I/2;
    g_prop = *g + R::rnorm(0, prop_sd);
    g_prop = g_prop>thresh ? g_prop : thresh + (thresh-g_prop); // reflect
  } else {
    g_prop = I;
  }
  
  // Get Bayes factor against null model + conditional pars
  arma::vec mu, mu_prop;
  arma::mat Sig, Sig_prop;

  double lBF1 = get_lBF0(I, p1, eps, x, vx, zx, mu, Sig,
                         intercept, *g);  
  double lBF2 = get_lBF0(I, p2, eps_prop, x, vx, zx, mu_prop, Sig_prop,
                         intercept, g_prop);

  // accept or reject proposal
  double lp_accept, g_prior_ratio;
  if (g_hyper!=0){
    g_prior_ratio = 
      (3/2)*(log(*g+1)-log(g_prop+1)) +
      (1/2)*(log(p1+1)-log(p2+1));
  }
  else {
    g_prior_ratio=0;
  }
  lp_accept = (lBF2 - log(Rf_choose(p,p2)))-(lBF1 - log(Rf_choose(p,p1))) +
    g_prior_ratio;
  
  bool accept = log(R::runif(0,1)) < lp_accept;
  if (accept){
    eps = eps_prop;
    *g = g_prop;
    *accept_count=*accept_count + 1;
  }
  
  // update eta (sample from mvn normal)
  int p_new = accept ? p2 : p1;
  arma::vec eta_new(p_new);
  if (accept){
    eta_new =  mu_prop + arma::chol(Sig_prop)*arma::randn(p_new, 1);;
  } else {
    eta_new =  mu + arma::chol(Sig)*arma::randn(p_new, 1);;
  }
  //Rcout << "p: " << p_new << " g: " << *g << std::endl;
  
  eta.fill(0);
  eta.elem(arma::find(eps==1)) = eta_new;
}

/////////////////////////////////////////////////////////

void update_g(int I, 
              arma::uvec &eps,
              double intercept,
              arma::vec &linpred,
              double *g,
              int *accept_count3){
  
  double p = sum(eps);
  arma::vec x_eta = linpred - (intercept * arma::ones(I));
  arma::vec x_eta2 = x_eta.t() * x_eta;
  
  // propose new g
  double thresh = (I-p)/(p+1);
  double prop_sd = I/2;
  double g_prop = *g + R::rnorm(0, prop_sd);
  g_prop = g_prop>thresh ? g_prop : thresh + (thresh-g_prop); 
  
  // calc acceptance prob
  double prior_ratio, lik_ratio;
  prior_ratio = (3/2)*(log(*g+1)-log(g_prop+1));
  lik_ratio = (p/2)*(log(*g)-log(g_prop)) -
    ((1/ *g)-(1/g_prop))*(x_eta2[0]/8);
  
  // accept and update
  bool accept = log(R::runif(0,1)) < (prior_ratio + lik_ratio);
  if (accept){
    *g = g_prop;
    *accept_count3=*accept_count3 + 1;
  }
  //Rcout << "g: " << *g << std::endl;
}

/////////////////////////////////////////////////////////

void update_intercept(int I,
                      arma::uvec eps,
                      arma::vec &vx,
                      arma::vec &zx,
                      double *intercept,
                      arma::vec &linpred,
                      double int_prior_scale,
                      int *accept_count2){
  
  // amount to multiply proposal scale to match normal
  double mult = int_prior_fam == 0 ? 1 : 1.6;
  arma::vec x_eta = linpred - (*intercept * arma::ones(I));
  double sig2 = 1/(sum(vx) + 1/pow(mult*int_prior_scale,2));
  double mu = sig2 * ( dot(vx, zx - x_eta) );
  //Rcout << " mu: " << mu << " sig2" << sig2 << std::endl;
  
  double prop;
  if (int_prior_fam==0){ //normal
    prop = R::rnorm(mu, sqrt(sig2));
    *intercept = prop;
  } else if (int_prior_fam==1){ //logistic
    // calc acceptance prob
    double scale = sqrt(sig2); //sqrt(sig2/3);
    prop = mu + scale*R::rt(3);
    double curr = *intercept;
    double q_prop, f_prop, lik_prop=0, q_curr, f_curr, lik_curr=0;
    q_prop = R::dt((prop-mu)/scale, 3, true);
    q_curr = R::dt((curr-mu)/scale, 3, true);
    f_prop = R::dlogis(prop, 0, int_prior_scale, true);
    f_curr = R::dlogis(curr, 0, int_prior_scale, true);
    for (int i=0; i<I; i++){
      lik_prop += R::dnorm(zx[i], x_eta[i]+prop, 1/sqrt(vx[i]), true);
      lik_curr += R::dnorm(zx[i], x_eta[i]+curr, 1/sqrt(vx[i]), true);
    }
    
    // accept or reject
    double lp_accept = lik_prop + f_prop + q_curr - lik_curr - f_curr - q_prop;
    bool accept = log(R::runif(0,1)) < lp_accept;
    
    if (accept){
      *intercept=prop;
      *accept_count2=*accept_count2+1;
    }
    //Rcout << exp(lp_accept) << " " << accept << std::endl;
  } else {
    Rcout << "Invalid prior family for intercept" << std::endl;
    exit(1);
  }
  
}

void update_vx_zx(int I, 
                  Rcpp::LogicalVector &zeta, 
                  arma::vec &vx, 
                  arma::vec &zx,
                  arma::vec&linpred){
  for (int i=0; i<I; i++){
    vx[i] = rpg(1, linpred[i]);
    zx[i] = (zeta[i]-0.5)/vx[i];
  }
}

void update_theta(int I,
                  Rcpp::NumericVector &theta,
                  arma::uvec &eps,
                  arma::mat &x,
                  arma::vec &eta,
                  double intercept,
                  arma::vec &linpred){
  
  arma::uvec idx = arma::find(eps == 1);
  linpred = x.cols(idx) * eta.rows(idx) + intercept;
  
  for (int i=0; i<I; i++){
    theta[i] = exp(linpred[i])/(1+exp(linpred[i]));
  }
}

/////////////////////////////////////////////////////////

void run_covreg_onestep(int I,
                         Rcpp::LogicalVector &zeta,
                         Rcpp::NumericVector &theta,
                         arma::uvec &eps,
                         double *g,
                         arma::mat &x,
                         arma::vec &vx,
                         arma::vec &zx,
                         arma::vec &eta,
                         double *intercept,
                         arma::vec &linpred,
                         int *accept_count,
                         int *accept_count2,
                         int *accept_count3,
                         double int_prior_scale){
  
  int p = eta.n_elem;
  
  // update vx
  update_vx_zx(I, zeta, vx, zx, linpred);
  
  // update intercept
  update_intercept(I, eps, vx, zx, intercept, linpred, 
                   int_prior_scale, accept_count2);
  
  // update g, epsilon and eta
  update_eps_eta(I, p, eps, x, vx, zx, eta, *intercept, accept_count,
                   g, g_hyper);
  
  if (g_hyper!=0){
    update_g(I, eps, *intercept, linpred, g, accept_count3);
  }
  
  // update theta and linpred
  update_theta(I, theta, eps, x, eta, *intercept, linpred);
}