#ifndef gibbs_updates
#define gibbs_updates

/* a) sample missing data and PG random varibles */

void update_yvz(int I,int J, int gam0_max,
                    Rcpp::LogicalMatrix &y,
                    Rcpp::NumericMatrix &v,
                    Rcpp::NumericMatrix &z,
                    const Rcpp::LogicalMatrix &missing, 
                    Rcpp::NumericVector &mu, 
                    Rcpp::NumericVector &alf, 
                    Rcpp::NumericMatrix &beta);

/* b) sample policy parameters */

void update_mu(int I, int J, int gam0_max,
                double rho_mu, double kap2_mu, 
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta);

void update_alf(int I, int J, int gam0_max,
                double kap2_alf, double w_alf,
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta);


void update_zeta_beta(int I, int J, int gam0_max,
                double eta, double sig2, bool cut,
                Rcpp::NumericVector &theta,
                Rcpp::NumericMatrix &v,
                Rcpp::NumericMatrix &z,
                Rcpp::NumericVector &mu, 
                Rcpp::NumericVector &alf, 
                Rcpp::NumericMatrix &beta,
                Rcpp::LogicalVector &zeta);

/* c) enforce identifiability */

void enforce_identifiability(int I, int J, 
                             int demLeader, int repLeader, int gam0_max, 
                             Rcpp::NumericVector &mu,  
                             Rcpp::NumericVector &alf,  
                             Rcpp::NumericMatrix &beta);

/* d) sample hyperparameters */

void update_mu_hyperparams(Rcpp::NumericVector &mu, 
                      double *rho_mu, double *kap2_mu,
                      double rho_mu_mean, double rho_mu_var,
                      double kap2_mu_a, double kap2_mu_b);

void update_alf_hyperparams(int J, Rcpp::NumericVector &alf, 
                      double *kap2_alf, double *w_alf,
                      double kap2_alf_a, double kap2_alf_b,
                      double w_a, double w_b);

void update_theta_bb(Rcpp::LogicalVector &zeta, 
                     Rcpp::NumericVector &theta,  
                     double theta_a, double theta_b);

void update_beta_hyperparams(int I,
                      Rcpp::NumericMatrix &beta,
                      Rcpp::LogicalVector &zeta,
                      double *eta, double *sig2,
                      double eta_mean, double eta_var,
                      double sig2_a, double sig2_b);

#endif