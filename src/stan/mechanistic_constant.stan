functions {

  matrix construct_matrix(int n0, real a, real r) {

    vector[n0 + 1] a_ns;
    for(i in 0:n0) {
      a_ns[i + 1] = -(n0 - i) * a;
    }
    vector[n0 + 1] rs = rep_vector(r, n0 + 1);
    rs[1] = 0;

    matrix[n0 + 1, n0 + 1] top_left = diag_matrix(a_ns);
    matrix[n0 + 1, n0 + 1] top_right = diag_matrix(rs);
    matrix[n0 + 1, n0 + 1] bottom_right = diag_matrix(-rs);
    matrix[n0 + 1, n0 + 1] bottom_left = rep_matrix(0, n0 + 1, n0 + 1);
    int k = 1;
    for(i in 2:(n0 + 1)) {
      for(j in 1:(n0 + 1)) {
        if(i - j == 1) {
          bottom_left[i, j] = -a_ns[k];
          k += 1;
        }
      }
    }

    matrix[n0 + 1, 2 * (n0 +1)] A_top = append_col(top_left, top_right);
    matrix[n0 + 1, 2 * (n0 +1)] A_bottom = append_col(bottom_left, bottom_right);

    matrix[2 * (n0 + 1), 2 * (n0 +1)] A = append_row(A_top, A_bottom);

    return A;
  }

  vector solve_system_logp(real t, int n0, real a, real r,  vector p_zero) {
    int dim_matrix = 2 * (n0 + 1);
    matrix[dim_matrix, dim_matrix] A = construct_matrix(n0, a, r);
    matrix[dim_matrix, dim_matrix] exp_A = matrix_exp(A * t);
    vector[dim_matrix] p_both_states = exp_A * p_zero;
    int dim_p = n0 + 1;
    vector[dim_p] p;
    for(i in 1:dim_p) {
      p[i] = p_both_states[i] + p_both_states[i + dim_p]; // because system is stacked
    }

    return(log(p));
  }

  real overall_logp(array[] int n_prey_remaining, vector log_p_v) {
    real log_p = 0.0;
    int n = num_elements(n_prey_remaining);
    int n_lambda = num_elements(log_p_v);
    for(i in 1:n)
      log_p += log_p_v[n_lambda - n_prey_remaining[i]];
    return(log_p);
  }

  real stochastic_reaction_lpmf(
    array[] int n_prey_remaining, int n_prey_initial, real t_obs,
    real a, real r) {

      // construct helper vectors
      vector[n_prey_initial] n_possible_prey_remaining;
      vector[2 * (n_prey_initial + 1)] p_zero = rep_vector(0, 2 * (n_prey_initial + 1));
      p_zero[1] = 1;

      // calculate probability distribution
      vector[n_prey_initial + 1] log_p = solve_system_logp(t_obs, n_prey_initial, a, r, p_zero);

      // calculate log p(data|parameters)
      return(overall_logp(n_prey_remaining, log_p));
  }
}

data {
  int n; // number of data points
  array[n] int n_prey_remaining;
  real t_obs;

  // covariates
  int n_covariate_blocks;
  int n_covariates_a;
  int n_covariates_r;
  matrix[n_covariate_blocks, n_covariates_a] X_a;
  matrix[n_covariate_blocks, n_covariates_r] X_r;
  array[n_covariate_blocks] int n_prey_initial_distinct;
  array[n_covariate_blocks] int len_blocks;

  int T_sim;
  array[T_sim] int prey_density_sim;
  int n_covariate_blocks_simple_sim;
  matrix[n_covariate_blocks_simple_sim, n_covariates_a] X_a_sim;
  matrix[n_covariate_blocks_simple_sim, n_covariates_r] X_r_sim;
}

parameters {
  vector[n_covariates_a] beta_a;
  vector[n_covariates_r] beta_r;
}

transformed parameters {
  vector<lower=0>[n_covariate_blocks] a;
  vector<lower=0>[n_covariate_blocks] r;
  for(i in 1:n_covariate_blocks) {
    a[i] = exp(X_a[i] * beta_a);
    r[i] = exp(X_r[i] * beta_r);
  }
}

model {
  int pos = 1;
  for(i in 1:n_covariate_blocks) {
      segment(n_prey_remaining, pos, len_blocks[i]) ~ stochastic_reaction(
        n_prey_initial_distinct[i], t_obs, a[i], r[i]);
      pos += len_blocks[i];
  }

  beta_a ~ normal(0, 2.5);
  beta_r ~ normal(0, 2.5);
}

generated quantities {
  vector[n] log_likelihood;
  matrix[T_sim, n_covariate_blocks_simple_sim] prey_eaten_sim;
  vector[n_covariate_blocks_simple_sim] a_sim_unique;
  vector[n_covariate_blocks_simple_sim] r_sim_unique;
  vector[n_covariate_blocks_simple_sim] h_sim_unique;
  
  for(j in 1:n_covariate_blocks_simple_sim) {
    
    real a_sim = exp(X_a_sim[j] * beta_a);
    real r_sim = exp(X_r_sim[j] * beta_r);
    a_sim_unique[j] = a_sim;
    r_sim_unique[j] = r_sim;
    h_sim_unique[j] = 1 / r_sim;
    
    
    for(i in 1:T_sim) {
      int n_prey_initial_tmp = prey_density_sim[i];
      vector[n_prey_initial_tmp] n_possible_prey_remaining;
      vector[2 * (n_prey_initial_tmp + 1)] p_zero = rep_vector(0, 2 * (n_prey_initial_tmp + 1));
      p_zero[1] = 1;
  
      // calculate probability distribution
      vector[n_prey_initial_tmp + 1] log_p = solve_system_logp(t_obs, n_prey_initial_tmp, a_sim, r_sim, p_zero);
  
      array[n_prey_initial_tmp + 1] int n_possible_prey_remaining_longer;
      for(k in 1:(n_prey_initial_tmp + 1))
        n_possible_prey_remaining_longer[k] = n_prey_initial_tmp - k + 1;
 
      int idx = categorical_rng(softmax(log_p));
      int n_prey_remaining_sim = n_possible_prey_remaining_longer[idx];
      prey_eaten_sim[i, j] = n_prey_initial_tmp - n_prey_remaining_sim;
    }
  }
  
  int k = 1;
  for(i in 1:n_covariate_blocks) {

      // construct helper vectors
      int n_prey_initial_tmp = n_prey_initial_distinct[i];
      vector[n_prey_initial_tmp] n_possible_prey_remaining;
      vector[2 * (n_prey_initial_tmp + 1)] p_zero = rep_vector(0, 2 * (n_prey_initial_tmp + 1));
      p_zero[1] = 1;
      
      // calculate probability distribution
      vector[n_prey_initial_tmp + 1] log_p = solve_system_logp(t_obs, n_prey_initial_tmp, a[i], r[i], p_zero);

      for(j in 1:len_blocks[i]) {
        log_likelihood[k] = overall_logp({n_prey_remaining[k]}, log_p);
        k += 1;
      }

  }

}
