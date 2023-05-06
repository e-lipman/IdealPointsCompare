#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "configs.hpp"
#include "parameters.hpp"
#include "gibbs_updates.hpp"
#include "covariate_regression.hpp"
#include "helper_functions.hpp"

using namespace Rcpp;

//[[Rcpp::export]]
List run_sampler(const LogicalMatrix y_raw, const LogicalMatrix missing, 
                int demLeader, int repLeader, int gam0_max,
                int iter, int burnin, int thin, int progress,
                bool cut,
                Nullable<NumericMatrix> x_ = R_NilValue,
                int eps_init_method = 0) {
    
    /*get I and J*/
    int I = y_raw.nrow();
    int J = y_raw.ncol();
    bool xflag = x_.isNotNull();
    Rcout << "I: " << I << " J: " << J << std::endl;
    
    if (!xflag){
      Rcout << "theta prior params: a=" << theta_a << " b=" << theta_b << std::endl;
    }
  
    /*check arguments*/
    if(y_raw.nrow()!=missing.nrow() ||
       y_raw.ncol()!=missing.ncol()){
      throw std::invalid_argument("'y_raw' and 'missing' must have same dimensions");
    }
    if(demLeader>=I || repLeader>=I || gam0_max>=J-1){
      throw std::invalid_argument("invalid index for 'demLeader', 'repLeader', or 'gam0_max'");
    }
    if(cut & xflag){
      throw std::invalid_argument("can't run both stages of cut model simultaneously with this function");
    }
    
    
    /*initialize class containing parameters*/
    arma::mat x;
    int p = 0;
    Parameters *pp = NULL;
    
    double int_prior_scale=0;
    if (xflag){
      x = as<arma::mat>(NumericMatrix (x_));
      p = x.n_cols; 
      pp = new ParametersHier(I, J, x, gam0_max, p, y_raw, missing, 
                              eps_init_method);
      int_prior_scale = int_prior_scale_raw;
      if (int_prior_scale_raw==0){
        int_prior_scale=get_scale_boonstra(int_prior_fam, I, 0.01);
      }
      Rcout << "Running with " << p << " covariates" << std::endl;
      Rcout << "Intercept prior: family=" << int_prior_fam  << ", scale=" << int_prior_scale << std::endl;
    } else {
      pp = new ParametersNohier(I, J, gam0_max, y_raw, missing);
      Rcout << "Running without covariates" << std::endl;
    }
    
    ParametersHier *pp_x=dynamic_cast<ParametersHier*>(pp);
    if(pp_x == NULL && xflag==1){ exit(1); }
    
    /*params to store*/
    NumericMatrix beta0Out(iter, I);
    NumericMatrix beta1Out(iter, I);
    LogicalMatrix zetaOut(iter, I);
    NumericMatrix thetaOut(iter, xflag ? I : 1);
    NumericVector lpost(iter);
    NumericVector gOut(iter);
    
    NumericMatrix etaOut (iter, p);
    LogicalMatrix epsOut (iter, p);
    NumericVector interceptOut (iter);
    
    /*Gibbs sampler*/
    for (int k=0; k < burnin + iter; k++){
      if (k % progress == 0){
        print_datetime();
        Rcout << "iter " << k << std::endl;
      }
      if ((k==burnin) & xflag){ // reset counters
        pp_x->accept_count=0;
        pp_x->accept_count2=0;
        pp_x->accept_count3=0;
      }
      
      int rep_iter = k<burnin ? 1 : thin;
      for (int t=0; t<rep_iter; t++){
        R_CheckUserInterrupt();
        
        /* a) sample missing data and PG random varibles
              i) impute missing data (y)
              ii) sample PG latent variable (v)
              iii) compute transformed outcome (z) */
        update_yvz(I, J, gam0_max, 
                   pp->y, pp->v, pp->z, missing, pp->mu, pp->alf, pp->beta);
        
        /* b) sample policy parameters */
        update_mu(I, J, gam0_max, 
                  pp->rho_mu, pp->kap2_mu, pp->v, pp->z, 
                  pp->mu, pp->alf, pp->beta);
        update_alf(I, J, gam0_max, 
                  pp->kap2_alf, pp->w_alf, pp->v, pp->z, 
                  pp->mu, pp->alf, pp->beta);
        
        update_zeta_beta(I, J, gam0_max, 
                  pp->rho_beta, pp->kap2_beta, cut, 
                  pp->theta, pp->v, pp->z, 
                  pp->mu, pp->alf, pp->beta, pp->zeta);
        
        /* c) enforce identifiability */
        enforce_identifiability(I, J, demLeader, repLeader, gam0_max,
                                pp->mu, pp->alf, pp->beta);

        /* d) sample hyperparameters */
        update_mu_hyperparams(pp->mu, &pp->rho_mu, &pp->kap2_mu, 
                               rho_mu_mean, rho_mu_var, kap2_mu_a, kap2_mu_b);
        update_alf_hyperparams(J, pp->alf, &pp->kap2_alf, &pp->w_alf,
                              kap2_alf_a, kap2_alf_b, w_a, w_b);
        if (xflag){
          run_covreg_onestep(I, pp->zeta, pp->theta, 
                              pp_x->eps, &pp_x->g, x, pp_x->vx, pp_x->zx, 
                              pp_x->eta, &pp_x->intercept, pp_x->linpred,
                              &pp_x->accept_count, 
                              &pp_x->accept_count2, 
                              &pp_x->accept_count3, 
                              int_prior_scale);
        } else {
          update_theta_bb(pp->zeta, pp->theta, theta_a, theta_b); 
        }
        update_beta_hyperparams(I, pp->beta, pp->zeta, 
                                &pp->rho_beta, &pp->kap2_beta,
                              eta_mean, eta_var, sig2_a, sig2_b);
      }
      
      /* update outputs */
      if (k>=burnin){
          beta0Out(k-burnin,_) = pp->beta(_,0);
          beta1Out(k-burnin,_) = pp->beta(_,1);
          zetaOut(k-burnin,_) = pp->zeta;
          
          if (xflag){
            thetaOut(k-burnin,_) = pp->theta;
            interceptOut[k-burnin] = pp_x->intercept;
            for (int idx=0; idx<p; idx++){
              etaOut(k-burnin,idx) = pp_x->eta[idx];
              epsOut(k-burnin,idx) = pp_x->eps[idx];
            }
            gOut[k-burnin] = pp_x->g;
          } else {
            thetaOut(k-burnin,0) = pp->theta[0];
          }
          
          lpost[k-burnin] = compute_lpost(I, J, gam0_max,    
                                          y_raw, 
                                          missing, 
                                          pp->mu, pp->alf, pp->beta, 
                                          pp->zeta, pp->theta,    
                                          pp->rho_mu, pp->kap2_mu,   
                                          pp->kap2_alf, pp->w_alf,   
                                          pp->rho_beta, pp->kap2_beta,
                                          pp_x, int_prior_scale);
          //llik[k-burnin] = compute_llik(I, J, gam0_max, y_raw, missing, 
          //                                pp->mu, pp->alf, pp->beta);
      }
    }
    
    if (xflag){
      double accept_rate = (double)  pp_x->accept_count/(iter*thin);
      double accept_rate2 = (double)  pp_x->accept_count2/(iter*thin);
      double accept_rate3 = (double)  pp_x->accept_count3/(iter*thin);
      
      return List::create(Named("beta0") = beta0Out, 
                          Named("beta1") = beta1Out,
                          Named("zeta") = zetaOut,
                          //Named("theta") = thetaOut,
                          Named("lpost") = lpost,
                          Named("eta") = etaOut,
                          Named("intercept") = interceptOut,
                          //Named("eps") = epsOut,
                          Named("g") = gOut,
                          Named("accept_eps") = accept_rate,
                          Named("accept_int") = accept_rate2,
                          Named("accept_g") = accept_rate3);      
    } else if (cut) {
      return List::create(Named("beta0") = beta0Out, 
                          Named("beta1") = beta1Out,
                          Named("zeta") = zetaOut,
                          Named("theta") = thetaOut);
    } else {
      return List::create(Named("beta0") = beta0Out, 
                          Named("beta1") = beta1Out,
                          Named("zeta") = zetaOut,
                          Named("theta") = thetaOut,
                          Named("lpost") = lpost);
    }

}

//[[Rcpp::export]]
List run_covreg_only(Rcpp::LogicalVector &zeta,
                     arma::mat &x,
                     int iter, int burnin, int thin, int progress,
                     int eps_init_method = 0, 
                     bool quiet=false){
  
  int I = zeta.length();
    
  int p = x.n_cols;
  if (not quiet){
    Rcout << "I: " << I << " p: " << p - 1 << std::endl;
  }
  
  double int_prior_scale = int_prior_scale_raw;
  if (int_prior_scale_raw==0){
    int_prior_scale=get_scale_boonstra(int_prior_fam, I, 0.01);
  }
  Rcout << "Intercept prior: family=" << int_prior_fam  << ", scale=" << int_prior_scale << std::endl;
  
  // initialize
  ParametersHier *pp;
  pp = new ParametersHier(I, p, x, eps_init_method);
  
  NumericMatrix thetaOut (iter, I);
  NumericMatrix etaOut (iter, p);
  LogicalMatrix epsOut (iter, p);
  NumericVector interceptOut (iter);
  NumericVector lpost(iter);
  NumericVector gOut(iter);
  
  // run mcmc
  for (int k=0; k < burnin + iter; k++){
    
    if (not quiet and (k % progress == 0)){
      print_datetime();
      Rcout << "iter " << k << std::endl;
    }
    
    int rep_iter = k<burnin ? 1 : thin;
    for (int t=0; t<rep_iter; t++){
    
    run_covreg_onestep(I, zeta, pp->theta, 
                       pp->eps, &pp->g, x, pp->vx, pp->zx, 
                       pp->eta, &pp->intercept, pp->linpred, 
                       &pp->accept_count, 
                       &pp->accept_count2, 
                       &pp->accept_count3, 
                       int_prior_scale);
    }
    if (k>=burnin){
      lpost[k-burnin] = compute_lpost_zeta_hier(I, p, zeta, pp-> eps, 
                                                x, pp->eta, 
                                                pp->intercept, pp->theta,
                                                int_prior_scale);
      gOut[k-burnin] = pp->g;
      thetaOut(k-burnin,_) = pp->theta;
      interceptOut[k-burnin] = pp->intercept;
      for (int idx=0; idx<p; idx++){
        etaOut(k-burnin,idx) = pp->eta[idx];
        epsOut(k-burnin,idx) = pp->eps[idx];
      }
    }
  }
  return List::create(Named("theta") = thetaOut,
                      Named("eta") = etaOut,
                      Named("intercept") = interceptOut,
                      Named("eps") = epsOut,
                      Named("lpost") = lpost,
                      Named("gOut") = gOut);  
}






// //[[Rcpp::export]]
// List run_covreg_only_multi(Rcpp::LogicalMatrix &zeta, 
//                            Rcpp::NumericVector &intercept,
//                            arma::mat &x, 
//                            int steps, int progress, 
//                            int eps_init_method = 0){
//   
//   int I = zeta.ncol();
//   int iter = zeta.nrow();
//   int p = x.n_cols;
//   Rcout << "I: " << I << " iter: " << iter  <<" p: " << p << std::endl;
//   
//   // initialize
//   NumericMatrix thetaOut (iter, I);
//   NumericMatrix etaOut (iter, p);
//   LogicalMatrix epsOut (iter, p);
//   
//   NumericVector theta (I);
//   LogicalVector eps (p);
//   NumericVector eta (p);
//   
//   for (int k=0; k<iter; k++){
//     if (k % progress == 0){
//       print_datetime();
//       Rcout << "iter " << k << std::endl;
//     }
//     
//     LogicalVector zeta_k = zeta(k,_);
//     List res_list = run_covreg_only(zeta_k, x,
//                                1, steps-1, 1, steps, 
//                                eps_init_method, 
//                                true, true, intercept[k]);
//     
//     theta = res_list["theta"];
//     eps = res_list["eps"];
//     eta = res_list["eta"];
//     thetaOut(k,_) = theta;
//     epsOut(k,_) = eps;
//     etaOut(k,_) = eta;
//   }
//   
//   return List::create(Named("theta") = thetaOut,
//                       Named("eta") = etaOut,
//                       Named("eps") = epsOut);
// }
