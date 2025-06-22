
prepare_stan_data_covariates <- function(df, levels_h=NULL, levels_a=NULL, T_sim=100) {

  if(is.null(levels_h))
    levels_h = rep(1, nrow(df))
  if(is.null(levels_a))
    levels_a = rep(1, nrow(df))
  
  n_levels_h = max(levels_h)
  n_levels_a = max(levels_a)
  
  levels_sim <- tibble(
    levels_h=1:n_levels_h,
    levels_a=1:n_levels_a
  ) %>% 
    as.matrix()
  
  prey_density_sim = seq(1, max(df$n_prey_initial))
  T_sim <- length(prey_density_sim)
  
  data_stan <- list(
    T=nrow(df),
    prey_eaten=df$n_eaten,
    prey_density=df$n_prey_initial,
    n_levels_h=n_levels_h,
    n_levels_a=n_levels_a,
    levels_a=levels_a,
    levels_h=levels_h,
    T_sim=T_sim,
    prey_density_sim=prey_density_sim,
    n_levels_sim=nrow(levels_sim),
    levels_sim=levels_sim
  )
  
  data_stan
}

copy_column_str <- function(df, from_col_name, to_col_name) {
  from_sym <- sym(from_col_name)
  to_sym <- sym(to_col_name)
  
  df %>%
    mutate(!!to_sym := !!from_sym)
}

extract_quantile_longer <- function(quantile_level, functional_response_sim) {
  apply(functional_response_sim, c(2, 3), function(x) quantile(x, quantile_level)) %>% 
    as.data.frame() %>% 
    mutate(n_prey_initial=1:nrow(.)) %>% 
    pivot_longer(-n_prey_initial) %>% 
    mutate(name=as.numeric(as.factor(name))) %>% 
    rename(covariate_level=name)
}

plot_posterior_predictive_stan <- function(fit, df, covar_name=NULL) {
  
  functional_response_sim = rstan::extract(fit, "prey_eaten_sim")[[1]]
  lower <- extract_quantile_longer(0.025, functional_response_sim) %>% 
    rename(lower=value)
  upper <- extract_quantile_longer(0.975, functional_response_sim) %>% 
    rename(upper=value)
  middle <- extract_quantile_longer(0.5, functional_response_sim) %>% 
    rename(middle=value)
  lower_1 <- extract_quantile_longer(0.25, functional_response_sim) %>% 
    rename(lower_1=value)
  upper_1 <- extract_quantile_longer(0.75, functional_response_sim) %>% 
    rename(upper_1=value)
  df_qs <- lower %>% 
    left_join(upper) %>% 
    left_join(middle) %>% 
    left_join(lower_1) %>% 
    left_join(upper_1)
  
  if(!is.null(covar_name)) {
    df_1 <- copy_column_str(df, covar_name, "covariate") %>% 
      mutate(covariate_level=as.numeric(covariate))
    
    lookup_levels <- df_1 %>% 
      select(covariate_level, covariate) %>% 
      unique()
    
    df_qs <- df_qs %>% 
      left_join(lookup_levels)
  } else {
    df_1 <- df %>% 
      mutate(covariate="")
  }
  
  g <- ggplot(df_1, aes(x=n_prey_initial)) +
    geom_ribbon(data=df_qs, aes(ymin=lower, ymax=upper),
                fill="blue", alpha=0.5) +
    geom_ribbon(data=df_qs, aes(ymin=lower_1, ymax=upper_1),
                fill="blue", alpha=0.5) +
    geom_line(data=df_qs, aes(y=middle), colour="blue") +
    geom_jitter(aes(y=n_eaten), height = 0, width = 0.25) +
    geom_smooth(aes(y=n_eaten), se=FALSE, colour="black",
                linetype = "dashed",
                linewidth = 0.5) +
    xlab("# initial prey") +
    ylab("# prey eaten")
  
  if(!is.null(covar_name))
    g <- g + facet_wrap(~covariate)
  
  g
}
