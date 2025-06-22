# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan", "loo")
)
options(mc.cores=4)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("src/r/")

# Replace the target list below with your own:
list(
  
  # datasets from frair r package: https://github.com/dpritchard/frair/tree/master/data
  # "size" here refers to the prey size (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12784)
  tar_target(data_bythotrephes, {
    load("data/raw/bythotrephes.RData")
    df <- bythotrephes
    df %>% 
      rename(
        n_prey_initial=density,
        n_eaten=eaten,
        n_alive=alive
        )
  }),
  tar_target(data_gammarus, {
    load("data/raw/gammarus.RData")
    df <- gammarus %>% 
      rename(
        n_prey_initial=density,
        n_eaten=eaten,
        n_alive=alive,
        species=spp
      ) %>% 
      arrange(species)
    species_names <- unique(df$species)
    df_1 <- df %>% 
      filter(species==species_names[1])
    df_2 <- df %>% 
      filter(species==species_names[2])
    list(
      celticus=df_1,
      pulex=df_2
    )
  }),
  tar_target(files_data_frair, {
    write_csv(data_bythotrephes, "data/processed/bythotrephes.csv")
    write_csv(data_gammarus$celticus, "data/processed/celticus.csv")
    write_csv(data_gammarus$pulex, "data/processed/pulex.csv")
  }),
  
  # data from: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13039
  tar_target(data_archer, {
    
    # dataset contains rows without predators present to measure natural mortality
    # will need a model that contains both predator-dependent and predator-independent
    # mortality
    read.csv("data/raw/data_archer_manual.csv") %>% 
      rename(
        n_prey_initial=N0,
        n_predators=PredNo
      ) %>% 
      mutate(
        n_alive=n_prey_initial-Ndead
      ) %>% 
      rename(n_dead=Ndead) %>% 
      select(n_prey_initial, n_dead, n_alive, n_predators)
  }),
  tar_target(data_sentis, {
    df <- read.csv("data/raw/data_sentis.csv") %>% 
      rename(
        n_prey_initial=Daphnia_density,
        n_prey_eaten=Daphnia_eaten
      ) %>% 
      mutate(
        n_prey_alive=n_prey_initial-n_prey_eaten
      )
    df
  }),
  tar_target(data_fw, {
    df <- read.csv("data/raw/rawFRdata.csv") %>% 
      rename(
        species_predator=predator,
        n_prey_initial=prey.density,
        n_eaten=eaten,
        n_alive=alive
      ) %>% 
      select(n_prey_initial, n_eaten, n_alive,
             species_predator, refuge) %>% 
      arrange(species_predator, refuge)
    list(
      gdc=df %>% filter(species_predator=="Gdc"),
      gp=df %>% filter(species_predator=="Gp")
    )
  }),
  tar_target(data_vulcic, {
    df <- read.csv("data/raw/data_vucic-pestic_manual.csv")
    df %>% 
      rename(
        n_prey_initial=N0,
        n_prey_eaten=Neaten
      ) %>% 
      mutate(
        n_prey_alive=n_prey_initial-n_prey_eaten
      )
  }),
  tar_target(files_data_odes, {
    write_csv(data_archer, "data/processed/archer.csv")
    write_csv(data_sentis, "data/processed/sentis.csv")
    write_csv(data_fw$gdc, "data/processed/gdc.csv")
    write_csv(data_fw$gp, "data/processed/gp.csv")
    write_csv(data_vulcic, "data/processed/vulcic.csv")
  }),
  
  
  # bythotrephes
  
  ## Fit basic binomial type-II model
  tar_target(stan_model_binomial_raw, "src/stan/binomial.stan",
             format = "file"),
  tar_target(stan_model_binomial, stan_model(stan_model_binomial_raw)),
  tar_target(stan_data_bythotrephes, {
    
    df <- data_bythotrephes %>% 
      mutate(
        levels=as.numeric(size)
      )
    levels_h <- df$levels
    levels_a <- df$levels
    prepare_stan_data_covariates(
       data_bythotrephes,
       levels_h=levels_h,
       levels_a=levels_a)
    }),
  tar_target(fit_bythotrephes_binomial,
             sampling(stan_model_binomial,
                      data=stan_data_bythotrephes,
                      iter=200,
                      chains=4)),
  tar_target(plot_fit_bythotrephes_binomial,
             plot_posterior_predictive_stan(
               fit_bythotrephes_binomial,
               data_bythotrephes,
               "size")),
  
  ## fit mechanistic model
  tar_target(stan_model_mechanistic_raw, "src/stan/mechanistic_constant.stan",
             format = "file"),
  tar_target(stan_model_mechanistic, stan_model(stan_model_mechanistic_raw)),
  tar_target(stan_data_mechanistic_bythotrephes, {
    
    df <- data_bythotrephes %>% 
      arrange(size, n_prey_initial) %>% 
      group_by(size, n_prey_initial) %>% 
      mutate(group_index = cur_group_id()) %>% 
      ungroup()
    df_sum <- df %>% 
      group_by(group_index) %>% 
      summarise(
        n=n()
      )
    
    df_covar <- df %>% 
      select(size, n_prey_initial) %>% 
      unique()
    
    X <- model.matrix(~ size - 1, data = df_covar)
    
    prey_density_sim <- seq(1, max(df$n_prey_initial))
    X_a_sim <- unique(X)
    
    list(
      n = nrow(df),
      n_prey_remaining = df$n_alive,
      t_obs = 1,
      n_covariate_blocks = nrow(X),
      n_covariates_a = ncol(X),
      n_covariates_r = ncol(X),
      X_a = X,
      X_r = X,
      n_prey_initial_distinct = df_covar$n_prey_initial,
      len_blocks = df_sum$n,
      prey_density_sim = prey_density_sim,
      T_sim = length(prey_density_sim),
      X_a_sim = X_a_sim,
      X_r_sim = X_a_sim,
      n_covariate_blocks_simple_sim = nrow(X_a_sim)
    )
  }),
  tar_target(fit_bythotrephes_mechanistic,
             sampling(stan_model_mechanistic,
                      data=stan_data_mechanistic_bythotrephes,
                      iter=200,
                      chains=4)),
  tar_target(plot_fit_bythotrephes_mechanistic,
             plot_posterior_predictive_stan(
               fit_bythotrephes_mechanistic,
               data_bythotrephes,
               "size")),
  tar_target(loo_bythotrephes, {
    
    log_like_binomial <- loo::extract_log_lik(fit_bythotrephes_binomial, "log_likelihood")
    log_like_mechanistic <- loo::extract_log_lik(fit_bythotrephes_mechanistic, "log_likelihood")
    
    loo_binomial <- loo::loo(log_like_binomial)
    loo_mechanistic <- loo::loo(log_like_mechanistic)
    loo::loo_compare(loo_binomial, loo_mechanistic)
  }),
  
  # Gammarus celticus (native species of parasite preying on dipteran larvae)
  # original data here (with more untapped datasets + more info on data e.g. replicates, experimental arms):
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12292 (which indicates some predator-independent prey mortality albeit at low levels)
  ## binomial
  tar_target(stan_data_gammarus_celticus, {
    
    df <- data_gammarus$celticus %>% 
      mutate(
        levels=1
      )
    levels_h <- df$levels
    levels_a <- df$levels
    prepare_stan_data_covariates(
      data_gammarus$celticus,
      levels_h=levels_h,
      levels_a=levels_a)
  }),
  tar_target(fit_gammarus_celticus_binomial,
             sampling(stan_model_binomial,
                      data=stan_data_gammarus_celticus,
                      iter=200,
                      chains=4)),
  tar_target(plot_fit_gammarus_celticus_binomial,
             plot_posterior_predictive_stan(
               fit_gammarus_celticus_binomial,
               data_gammarus$celticus)),
  
  ## mechanistic model
  tar_target(stan_data_mechanistic_gammarus_celticus, {
    
    df <- data_gammarus$celticus %>% 
      arrange(n_prey_initial) %>% 
      group_by(n_prey_initial) %>% 
      mutate(group_index = cur_group_id()) %>% 
      ungroup()
    df_sum <- df %>% 
      group_by(group_index) %>% 
      summarise(
        n=n()
      )
    
    df_covar <- df %>% 
      select(n_prey_initial) %>% 
      unique()
    
    X <- rep(1, nrow(df_covar)) %>% 
      as.matrix()
    
    prey_density_sim <- seq(1, max(df$n_prey_initial))
    X_a_sim <- unique(X)
    
    list(
      n = nrow(df),
      n_prey_remaining = df$n_alive,
      t_obs = 1,
      n_covariate_blocks = nrow(X),
      n_covariates_a = ncol(X),
      n_covariates_r = ncol(X),
      X_a = X,
      X_r = X,
      n_prey_initial_distinct = df_covar$n_prey_initial,
      len_blocks = df_sum$n,
      prey_density_sim = prey_density_sim,
      T_sim = length(prey_density_sim),
      X_a_sim = X_a_sim,
      X_r_sim = X_a_sim,
      n_covariate_blocks_simple_sim = nrow(X_a_sim)
    )
  }),
  tar_target(fit_gammarus_celticus_mechanistic,
             sampling(stan_model_mechanistic,
                      data=stan_data_mechanistic_gammarus_celticus,
                      iter=200,
                      chains=4)),
  tar_target(plot_gammarus_celticus_mechanistic,
             plot_posterior_predictive_stan(
               fit_gammarus_celticus_mechanistic,
               data_gammarus$celticus)),
  tar_target(loo_gammarus_celticus, {
    
    log_like_binomial <- loo::extract_log_lik(fit_gammarus_celticus_binomial, "log_likelihood")
    log_like_mechanistic <- loo::extract_log_lik(fit_gammarus_celticus_mechanistic, "log_likelihood")
    
    loo_binomial <- loo::loo(log_like_binomial)
    loo_mechanistic <- loo::loo(log_like_mechanistic)
    loo::loo_compare(loo_binomial, loo_mechanistic)
  }),
  
  # Gammarus pulex (invasive species of parasite preying on dipteran larvae)
  ## binomial
  tar_target(stan_data_gammarus_pulex, {
    
    df <- data_gammarus$pulex %>% 
      mutate(
        levels=1
      )
    levels_h <- df$levels
    levels_a <- df$levels
    prepare_stan_data_covariates(
      data_gammarus$pulex,
      levels_h=levels_h,
      levels_a=levels_a)
  }),
  tar_target(fit_gammarus_pulex_binomial,
             sampling(stan_model_binomial,
                      data=stan_data_gammarus_pulex,
                      iter=200,
                      chains=4)),
  tar_target(plot_fit_gammarus_pulex_binomial,
             plot_posterior_predictive_stan(
               fit_gammarus_pulex_binomial,
               data_gammarus$pulex)),
  
  ## mechanistic model
  tar_target(stan_data_mechanistic_gammarus_pulex, {
    
    df <- data_gammarus$pulex %>% 
      arrange(n_prey_initial) %>% 
      group_by(n_prey_initial) %>% 
      mutate(group_index = cur_group_id()) %>% 
      ungroup()
    df_sum <- df %>% 
      group_by(group_index) %>% 
      summarise(
        n=n()
      )
    
    df_covar <- df %>% 
      select(n_prey_initial) %>% 
      unique()
    
    X <- rep(1, nrow(df_covar)) %>% 
      as.matrix()
    
    prey_density_sim <- seq(1, max(df$n_prey_initial))
    X_a_sim <- unique(X)
    
    list(
      n = nrow(df),
      n_prey_remaining = df$n_alive,
      t_obs = 1,
      n_covariate_blocks = nrow(X),
      n_covariates_a = ncol(X),
      n_covariates_r = ncol(X),
      X_a = X,
      X_r = X,
      n_prey_initial_distinct = df_covar$n_prey_initial,
      len_blocks = df_sum$n,
      prey_density_sim = prey_density_sim,
      T_sim = length(prey_density_sim),
      X_a_sim = X_a_sim,
      X_r_sim = X_a_sim,
      n_covariate_blocks_simple_sim = nrow(X_a_sim)
    )
  }),
  tar_target(fit_gammarus_pulex_mechanistic,
             sampling(stan_model_mechanistic,
                      data=stan_data_mechanistic_gammarus_pulex,
                      iter=200,
                      chains=4)),
  tar_target(plot_gammarus_pulex_mechanistic,
             plot_posterior_predictive_stan(
               fit_gammarus_pulex_mechanistic,
               data_gammarus$pulex)),
  tar_target(loo_gammarus_pulex, {
    
    log_like_binomial <- loo::extract_log_lik(fit_gammarus_pulex_binomial, "log_likelihood")
    log_like_mechanistic <- loo::extract_log_lik(fit_gammarus_pulex_mechanistic, "log_likelihood")
    
    loo_binomial <- loo::loo(log_like_binomial)
    loo_mechanistic <- loo::loo(log_like_mechanistic)
    loo::loo_compare(loo_binomial, loo_mechanistic)
  })
)
