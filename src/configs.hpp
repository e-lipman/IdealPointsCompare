#ifndef configs
#define configs

/* mu hyperpriors */
const double rho_mu_mean = 0;
const double rho_mu_var = 1;
const double kap2_mu_a = 2;
const double kap2_mu_b = 1;

/* alpha hyperpriors */
const double kap2_alf_a = 2;
const double kap2_alf_b = 1;
const double w_a = 1;
const double w_b = 1;

/* beta hyperpriors */
const double eta_mean = 0;
const double eta_var = 1;
const double sig2_a = 1.5;       
const double sig2_b = 0.5;

/* zeta hyperpriors */

//stage 1
const double theta_a = 1;
const double theta_b = 1; 

//stage 2
const double int_prior_scale_raw = 1; // (sd) normal 1.6 gives approx uniform
const int int_prior_fam = 1; // 0 = normal, 1 = logistic

const int g_hyper = 0; // 0 = unit info, 1=robust

#endif