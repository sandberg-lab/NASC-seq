
data {
  int<lower=0> reads;
  int<lower=0> content[reads];
  int<lower=0> conversions[reads];
  real<lower=0, upper=1> p_c;
  real<lower=0, upper=1> p_e;
}

parameters {
  real<lower=-3, upper=3> alpha_logged;
  real<lower=-3, upper=3> beta_logged;
  real<lower=0, upper=1> pi_g;
}

transformed parameters{
    real alpha;
    real beta;
    alpha = exp(alpha_logged);
    beta = exp(beta_logged);
}


model {
    pi_g ~ beta(alpha, beta);
    for (i in 1:reads){
        target += log_sum_exp(binomial_lpmf(conversions[i] | content[i], p_c) + bernoulli_lpmf(1|pi_g),binomial_lpmf(conversions[i] | content[i],p_e) + bernoulli_lpmf(0|pi_g));
    }
}
