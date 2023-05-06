#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "helper_functions.hpp"

using namespace Rcpp;

/* a) sample missing data and PG random varibles */

void update_yvz(int I,int J, int gam0_max,
                    Rcpp::LogicalMatrix &y,
                    Rcpp::NumericMatrix &v,
                    Rcpp::NumericMatrix &z,
                    const Rcpp::LogicalMatrix &missing, 
                    Rcpp::NumericVector &mu, 
                    Rcpp::NumericVector &alf, 
                    Rcpp::NumericMatrix &beta){
  
    double psi, p;
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            psi = mu[j] + alf[j]*beta(i, j<=gam0_max ? 0 : 1);
            if(missing(i, j) == 1){
                p = exp(psi)/(1+exp(psi));
                y(i, j) = R::rbinom(1,p);}
            v(i,j) = rpg(1,psi);
            z(i,j) = (y(i,j)-.5)/v(i,j);
            //Rcout << "y: " << y(i,j) << " v: " << v(i,j) << " z: " << z(i,j) <<  std::endl;
        }
    }
}

/* b) sample policy parameters */
void update_mu(int I, int J, int gam0_max,
                double rho_mu, double kap2_mu, 
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta){
    double rho_hat, kap2_hat, kap2_hat_inv;
    for (int j=0; j<J; j++){
        int gam = j<=gam0_max ? 0 : 1;
        rho_hat = rho_mu/kap2_mu;
        kap2_hat_inv = 1/kap2_mu;
        for (int i=0; i<I; i++){
            rho_hat += v(i,j)*(z(i,j)-alf[j]*beta(i,gam));
            kap2_hat_inv += v(i,j);
        }
        kap2_hat = 1/kap2_hat_inv;
        rho_hat = rho_hat * kap2_hat;
        mu[j] = R::rnorm(rho_hat, sqrt(kap2_hat));
        //Rcout << "MU | rho_hat: " << rho_hat << " kap2_hat: " << kap2_hat << " mu : " << mu(j) <<  std::endl;
    }
}

void update_alf(int I, int J, int gam0_max,
                double kap2_alf, double w_alf,
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta){
    double rho_hat, kap2_hat, kap2_hat_inv, w_hat, p0, p1;
    for (int j=0; j<J; j++){
        int gam = j<=gam0_max ? 0 : 1;
        rho_hat = 0;
        kap2_hat_inv = 1/kap2_alf;
        for (int i=0; i<I; i++){
            rho_hat += v(i,j)*(z(i,j)-mu[j])*beta(i, gam);
            kap2_hat_inv += v(i,j)*pow(beta(i,gam),2);
        }
        kap2_hat = 1/kap2_hat_inv;
        rho_hat = rho_hat * kap2_hat;

        p1 = w_alf*sqrt(kap2_alf);
        p0 = (1-w_alf)*sqrt(kap2_hat) *
                exp(0.5*pow(rho_hat,2)/kap2_hat);
        w_hat = p1/(p1+p0);
        if (R::rbinom(1,w_hat)==1){
            alf[j]=0;
        } else {
            alf[j]=R::rnorm(rho_hat, sqrt(kap2_hat));
        }

        //Rcout << "ALF | rho_hat: " << rho_hat << " kap2_hat: " << kap2_hat <<  " w_hat: " << w_hat << " alf : " << alf(j) <<  std::endl;
    }
}

void update_zeta_beta(int I, int J, int gam0_max,
                double eta, double sig2, bool cut, 
                Rcpp::NumericVector &theta,
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta,
                Rcpp::LogicalVector &zeta){
    
    double rhohat, rhohat0, rhohat1;
    double sig2hat, sig2hat0, sig2hat1;
    double sig2hat_inv, sig2hat0_inv, sig2hat1_inv;
    double rho_adder, sig2hat_inv_adder;
    double lp1=0, p1;

    for (int i=0; i<I; i++){
        sig2hat_inv = 1/sig2, sig2hat0_inv = 1/sig2, sig2hat1_inv = 1/sig2;
        rhohat = eta/sig2, rhohat0 = eta/sig2, rhohat1 = eta/sig2;
        for (int j=0; j<J; j++){
            rho_adder = v(i,j)*alf[j]*(z(i,j)-mu[j]);
            sig2hat_inv_adder = v(i,j)*pow(alf(j),2);
            
            rhohat += rho_adder;
            sig2hat_inv += sig2hat_inv_adder;

            if (j<=gam0_max){
                rhohat0 += rho_adder;
                sig2hat0_inv += sig2hat_inv_adder;
            } else {
                rhohat1 += rho_adder;
                sig2hat1_inv += sig2hat_inv_adder;
            }
        }

        sig2hat = 1/sig2hat_inv;
        sig2hat0 = 1/sig2hat0_inv;
        sig2hat1 = 1/sig2hat1_inv;
        
        rhohat = sig2hat* rhohat;
        rhohat0 = sig2hat0* rhohat0;
        rhohat1 = sig2hat1* rhohat1;
        
        lp1 =  
          0.5*( log(sig2 * sig2hat) - log(sig2hat0 * sig2hat1) ) +
          0.5*( pow(eta,2)/sig2 + pow(rhohat,2)/sig2hat - 
                pow(rhohat0,2)/sig2hat0 - pow(rhohat1,2)/sig2hat1);
        lp1+=log(theta[i]) - log(1-theta[i]);
        
        if (abs(lp1)>50){
          p1 = lp1>0 ? 1 : 0;
        } else {
          p1 = exp(lp1);
        }
        
        // sample zeta
        zeta[i] = R::rbinom(1, p1/(1+p1));
        
        //sample beta
        if (zeta[i]==1){
          beta(i,0) = R::rnorm(rhohat, sqrt(sig2hat));
          beta(i,1) = beta(i,0);
        } else {
          beta(i,0) = R::rnorm(rhohat0, sqrt(sig2hat0));
          beta(i,1) = R::rnorm(rhohat1, sqrt(sig2hat1));
        }
    }
}

/* c) enforce identifiability */

void enforce_identifiability(int I, int J, 
                             int demLeader, int repLeader, int gam0_max, 
                             Rcpp::NumericVector &mu,  
                             Rcpp::NumericVector &alf,  
                             Rcpp::NumericMatrix &beta){
    
    // pin beta0 anchors to +-1
    double beta_D = beta(demLeader,0);
    double beta_R = beta(repLeader,0);

    for (int i=0; i<I; i++){
        for (int gam=0; gam<=1; gam++){
            beta(i,gam) = (2*beta(i,gam)-beta_R-beta_D)/(beta_R-beta_D);
        }
    }
    for (int j=0; j<J; j++){
        mu[j] = mu[j] + 0.5*alf[j]*(beta_R+beta_D); 
        alf[j] = 0.5*alf[j]*(beta_R-beta_D);
    }
    
    //flip beta1 so oriented the same way as beta0
    
    bool flip1 = beta(demLeader,1)>beta(repLeader,1);
    if (flip1){
      for (int i=0; i<I; i++){
        beta(i,1)=-beta(i,1);
      }
      for (int j=gam0_max+1; j<J; j++){
        alf[j] = -1*alf[j];
      }
    }
}

/* d) sample hyperparameters */

static double update_normal_normal(Rcpp::NumericVector &x, 
                            double sig2, double mu_prior, double sig2_prior){
    
    int n = x.length();
    double sum_x = 0;
    for (int i=0; i<n; i++){ 
        sum_x += x[i];
    }
     
    double sig2_hat = 1 / (1/sig2_prior + n/sig2);
    double mu_hat = sig2_hat * (mu_prior/sig2_prior + sum_x/sig2);
    return(R::rnorm(mu_hat, sqrt(sig2_hat)));
}

static double update_normal_igam(Rcpp::NumericVector &x, 
                          double mu, double a_prior, double b_prior){
    int n = x.length();
    double sum_sq_dev = 0;
    for (int i=0; i<n; i++){ sum_sq_dev += pow(x[i]-mu, 2); }

    double a_hat = a_prior + n/2;
    double b_hat = b_prior + sum_sq_dev/2;
    return(1/(R::rgamma(a_hat, 1/b_hat))); // note: R::rgamma uses scale param
}

static double update_binom_beta(Rcpp::LogicalVector &x, 
                         double a_prior, double b_prior){
    int n = x.length(), k = 0;
    for (int i=0; i<n; i++){ if (x[i]==1){ k++; }; }
    
    double a_hat = a_prior + k;
    double b_hat = b_prior + (n-k);
    return(R::rbeta(a_hat, b_hat));
}

void update_mu_hyperparams(Rcpp::NumericVector &mu, 
                      double *rho_mu, double *kap2_mu,
                      double rho_mu_mean, double rho_mu_var,
                      double kap2_mu_a, double kap2_mu_b){

    *rho_mu = update_normal_normal(mu, *kap2_mu, rho_mu_mean, rho_mu_var);
    *kap2_mu = update_normal_igam(mu, *rho_mu, kap2_mu_a, kap2_mu_b);  
}

void update_alf_hyperparams(int J, Rcpp::NumericVector &alf, 
                      double *kap2_alf, double *w_alf,
                      double kap2_alf_a, double kap2_alf_b,
                      double w_a, double w_b){
    
    int count_nonzero=0;
    LogicalVector is_zero (J);
    for (int j=0; j<J; j++){ 
      is_zero[j] = alf[j]==0 ? 1 : 0;
      count_nonzero += 1-is_zero[j];
    }
    
    int idx_nonzero = 0;
    NumericVector alf_nonzero (count_nonzero);
    for (int j=0; j<J; j++){
      if (is_zero[j]==0){
        alf_nonzero[idx_nonzero]=alf[j];
        idx_nonzero++;
      }
    }
    
    *kap2_alf = update_normal_igam(alf_nonzero, 0, kap2_alf_a, kap2_alf_b);
    *w_alf = update_binom_beta(is_zero, w_a, w_b);
}

void update_theta_bb(Rcpp::LogicalVector &zeta, 
                  Rcpp::NumericVector &theta,  
                  double theta_a, double theta_b){
    
    double theta_val = update_binom_beta(zeta, theta_a, theta_b);                  
    theta = NumericVector(theta.length(), theta_val);
}

void update_beta_hyperparams(int I,
                      Rcpp::NumericMatrix &beta,
                      Rcpp::LogicalVector &zeta,
                      double *eta, double *sig2,
                      double eta_mean, double eta_var,
                      double sig2_a, double sig2_b){

    int count_zeta0 = 0;
    for (int i=0; i<I; i++){ 
        if (zeta[i]==0){ count_zeta0++; }
    }

    Rcpp::NumericVector beta_unique (I + count_zeta0);
    for (int i=0; i<I; i++){ beta_unique[i]=beta(i,0); }
    int idx = I;
    for (int i=0; i<I; i++){
        if (zeta[i]==0){
            beta_unique[idx] = beta(i,1);
            idx++;
        }
    }
    *eta = update_normal_normal(beta_unique, *sig2, eta_mean, eta_var);
    *sig2 = update_normal_igam(beta_unique, *eta, sig2_a, sig2_b);  
}

