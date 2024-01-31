#ifndef covariate_regression
#define covariate_regression

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
                        double int_prior_scale,
                        bool fix_int);

#endif