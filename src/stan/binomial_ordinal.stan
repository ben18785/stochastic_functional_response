functions {
  real rogers_II(real N0, real h, real a) {
    real term1 = h * N0 * a;
    real term2 = lambert_w0(exp(-a + term1) * term1);
    return (term1 - term2) / (h * a);
  }
}


data {
  // data to fit model
  int T; // # of rows in dataset
  array[T] int prey_density;
  array[T] int prey_eaten;
  int<lower=1> n_levels_h;
  int<lower=1> n_levels_a;
  array[T] int<lower=1> levels_h;
  array[T] int<lower=1> levels_a;
  
  // posterior predictive data
  int T_sim;
  array[T_sim] int prey_density_sim;
  int n_levels_sim;
  array[n_levels_sim, 2] int levels_sim;
}


parameters {
  ordered[n_levels_h] rev_logit_h;
  positive_ordered[n_levels_a] a;
}

transformed parameters {
  vector[T] functional_response;
  vector[n_levels_h] logit_h;
  for (i in 1:n_levels_h) {
    logit_h[i] = rev_logit_h[n_levels_h - i + 1];
  }
  vector[n_levels_h] h = inv_logit(logit_h);
  
  for(i in 1:T) {
    real h_temp = h[levels_h[i]];
    real a_temp = a[levels_a[i]];
    functional_response[i] = rogers_II(prey_density[i], h_temp, a_temp);
  }
}

model {
  for(i in 1:T) {
    prey_eaten[i] ~ binomial(prey_density[i], functional_response[i]/prey_density[i]);
  }
  
  // priors
  a ~ cauchy(0, 1);
  rev_logit_h ~ logistic(0, 1); // implies h ~ uniform(0, 1)
}

generated quantities {
  matrix[T_sim, n_levels_sim] functional_response_sim;
  matrix[T_sim, n_levels_sim] prey_eaten_sim;
  
  for(j in 1:n_levels_sim) {
    real h_temp = h[levels_sim[j, 1]];
    real a_temp = a[levels_sim[j, 2]];
    for(i in 1:T_sim) {
      functional_response_sim[i, j] = rogers_II(prey_density_sim[i], h_temp, a_temp);
      prey_eaten_sim[i, j] = binomial_rng(prey_density_sim[i], functional_response_sim[i, j] / prey_density_sim[i]);
    }
  }
}

