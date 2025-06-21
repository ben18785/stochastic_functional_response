# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("src/r/")

# Replace the target list below with your own:
list(
  
  # datasets from frair r package: https://github.com/dpritchard/frair/tree/master/data
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
  
  
  # Fit basic binomial type-II model
  tar_target(stan_model_binomial_raw, "src/stan/binomial_ordinal.stan",
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
  tar_target(fit_bythotrephes_ordinal,
             sampling(stan_model_binomial,
                      data=stan_data_bythotrephes,
                      iter=200,
                      chains=4)),
  tar_target(plot_fit_bythotrephes_ordinal,
             plot_posterior_predictive_stan(
               fit_bythotrephes_ordinal,
               stan_data_bythotrephes,
               data_bythotrephes,
               "size"))
)
