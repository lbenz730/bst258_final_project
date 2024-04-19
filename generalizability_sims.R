### Simulation Study of Generalizability Estimators

library(tidyverse)
library(furrr)
library(ranger)
library(ggridges)
source('helpers.R')
source('estimators.R')

### Function to generate generalizability data 
### Use data generation process similar to (Colnet, 2024)
### https://github.com/BenedicteColnet/combine-rct-rwd-review/blob/main/simulations/estimators_and_simulations-wo-cw.R
###
### n_trial = desired trial size (approximate)
### n_obs = desired size of observational study
### selection_model = list of coefficients for the selection model (into the trial)
### outcome_model = list of coefficients for the outcome model
### mult = multiplier to get the final trial sample size to ~ n_trial
generate_data <- function(n_trial, n_obs, selection_model, outcome_model, mult) {
  ### Generate Covariates from MVN
  n_subjects <- n_trial * mult  ### Expected Incidence
  df <- 
    tibble('X1' = rnorm(n_subjects, 1, 1),
           'X2' = rnorm(n_subjects, 1, 1),
           'X3' = rnorm(n_subjects, 1, 1),
           'X4' = rnorm(n_subjects, 1, 1))
  
  ### Sample trial probability and treatment probability (randomized w/in trial)
  df_trial <- 
    df %>% 
    mutate('prob_S' = expit(compute_model(., selection_model))) %>% 
    mutate('S' = rbinom(n_subjects, 1, prob_S))  %>% 
    filter(S == 1) %>% 
    mutate('A' = rbinom(n = nrow(.), size = 1, p = 0.5))
  
  ### Fresh Sample for observational study
  df_obs <- 
    tibble('X1' = rnorm(n_obs, 1, 1),
           'X2' = rnorm(n_obs, 1, 1),
           'X3' = rnorm(n_obs, 1, 1),
           'X4' = rnorm(n_obs, 1, 1),
           'S' = 0,
           'A' = NA)
  
  
  
  ### Stack Trial and Observational and generate outcomes
  df_final <- 
    bind_rows(df_trial, df_obs) %>% 
    select(-prob_S) %>% 
    mutate('Y' = compute_model(., outcome_model) + rnorm(n = nrow(.), mean = 0, sd = 1))
  
  return(df_final)
}

### Function to compute all estimators for a given dataset
compute_estimators <- function(df) {
  ### IPSW Estimators
  ipsw_glm <- ipsw(df, model_type = 'GLM', model_formula = S ~ X1 + X2 + X3 + X4)
  ipsw_rf <- ipsw(df, model_type = 'RF', model_formula = S ~ X1 + X2 + X3 + X4)
  
  ### G-Formula Estimators
  gformula_lm <- g_formula(df, model_type = 'LM', Y ~ A:X1 + X2 + X3 + X4)
  gformula_rf <- g_formula(df, model_type = 'RF', Y ~ X1 + X2 + X3 + X4 + AX1 + AX2 + AX3 + AX4)
  
  ### AIPSW Estimators
  aipsw_lm <- 
    aipsw(df, model_type = 'GLM/LM', 
          outcome_formula = Y ~ A:X1 + X2 + X3 + X4, 
          selection_formula = S ~ X1 + X2 + X3 + X4)
  
  aipsw_rf <- 
    aipsw(df, model_type = 'RF', 
          outcome_formula = Y ~ X1 + X2 + X3 + X4 + AX1 + AX2 + AX3 + AX4, 
          selection_formula = S ~ X1 + X2 + X3 + X4)
  
  df_results <- bind_rows(ipsw_glm, ipsw_rf, gformula_lm, gformula_rf, aipsw_lm, aipsw_rf)
  
  return(df_results)
}


### Run Simulation ###

### Set up parallelization and sim parameters
plan(multisession(workers = 32))
n_sims <- 1000

### Low S incidence
### a) Models Are Correctly Specified
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 50,
                  selection_model = list('(Intercept)' = -2.5,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 1391))

results_1a <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 380),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 2%',
         'model_spec' = 'Outcome/S Propensity Correct')

### Low S incidence
### b) Outcome Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 50,
                  selection_model = list('(Intercept)' = -2.5,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 131))

results_1b <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 30),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 2%',
         'model_spec' = 'Outcome Misspecified')

### Low S incidence
### c) Propensity Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 50,
                  selection_model = list('(Intercept)' = 0.5,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 11))

results_1c <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 3),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 2%',
         'model_spec' = 'S Propensity Misspecified')

### Low S incidence
### Propensity Model Incorrect and outcome model incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 50,
                  selection_model = list('(Intercept)' = 0.5,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 991))

results_1d <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 91943103),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 2%',
         'model_spec' = 'Outcome/S Propensity Misspecified')



### Medium S incidence
### a) Models Are Correctly Specified
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 5,
                  selection_model = list('(Intercept)' = 0,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 1391))

results_2a <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 380),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 20%',
         'model_spec' = 'Outcome/S Propensity Correct')

### Medium S incidence
### b) Outcome Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 5,
                  selection_model = list('(Intercept)' = 0,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 131))

results_2b <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 30),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 20%',
         'model_spec' = 'Outcome Misspecified')

### Medium S incidence
### c) Propensity Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 5,
                  selection_model = list('(Intercept)' = 3.6,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 11))

results_2c <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 3),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 20%',
         'model_spec' = 'S Propensity Misspecified')

### Medium S incidence
### Propensity Model Incorrect and outcome model incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 5,
                  selection_model = list('(Intercept)' = 3.6,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 991))

results_2d <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 91943103),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 20%',
         'model_spec' = 'Outcome/S Propensity Misspecified')


### High S incidence
### a) Models Are Correctly Specified
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 2,
                  selection_model = list('(Intercept)' = 1.7,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 1391))

results_3a <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 380),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 50%',
         'model_spec' = 'Outcome/S Propensity Correct')

### High S incidence
### b) Outcome Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 2,
                  selection_model = list('(Intercept)' = 1.7,
                                         'X1' = -0.5, 
                                         'X2' = -0.3, 
                                         'X3' = -0.5, 
                                         'X4' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 131))

results_3b <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 30),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 50%',
         'model_spec' = 'Outcome Misspecified')

### High S incidence
### c) Propensity Model Incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 2,
                  selection_model = list('(Intercept)' = 6.5,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 11))

results_3c <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 3),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 50%',
         'model_spec' = 'S Propensity Misspecified')

### High S incidence
### Propensity Model Incorrect and outcome model incorrect
dfs <- 
  future_map(1:n_sims, ~{
    generate_data(n_trial = 1000,
                  n_obs = 5000,
                  mult = 2,
                  selection_model = list('(Intercept)' = 6.5,
                                         'I(exp(X1))' = -0.5, 
                                         'I(exp(X2))' = -0.3, 
                                         'I(exp(X3))' = -0.5, 
                                         'I(exp(X4))' = -0.4),
                  outcome_model = list('(Intercept)' = -100,
                                       'A:X1:X2' = 27.4,
                                       'X2' = 13.7,
                                       'X3' = 13.7,
                                       'X4' = 13.7))
  }, .options = furrr_options(seed = 991))

results_3d <- 
  future_map_dfr(dfs, compute_estimators, 
                 .options = furrr_options(seed = 91943103),
                 .id = 'dataset_id') %>% 
  mutate('s_incidence' = 'Baseline S Incidence ~ 50%',
         'model_spec' = 'Outcome/S Propensity Misspecified')



### Plot Results
df_results <- 
  bind_rows(results_1a, results_1b, results_1c, results_1d,
            results_2a, results_2b, results_2c, results_2d,
            results_3a, results_3b, results_3c, results_3d) %>% 
  mutate('model_spec' = factor(model_spec, levels = c('Outcome/S Propensity Correct',
                                                      'Outcome Misspecified',
                                                      'S Propensity Misspecified',
                                                      'Outcome/S Propensity Misspecified')))
write_csv(df_results, 'sim_results/generalizability.csv')

ggplot(df_results %>% filter(tau_hat >= 0, tau_hat <= 40), aes(x = tau_hat, y = method)) + 
  facet_grid(s_incidence ~ model_spec) + 
  scale_x_continuous(limits = c(5, 35)) +
  geom_vline(xintercept = 27.4, lty = 2, lwd = 1.2, col = 'grey') + 
  geom_density_ridges(aes(fill = method), 
                      alpha = 0.5,
                      scale = 1,
                      rel_min_height = 0.02,
                      quantile_lines = T, 
                      quantiles = 0.5, 
                      panel_scaling = F) + 
  theme(legend.position = 'none') +
  labs(x = 'Estimated ATE',
       y = 'Method',
       title = 'Simulation Study of Generalizability Estimators')

ggsave('figures/generalizability.png', height = 12, width = 16)


