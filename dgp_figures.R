library(tidyverse)
library(ggridges)
library(ranger)
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

set.seed(19217)
df1 <- 
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
                                     'X4' = 13.7)) %>% 
  mutate('kappa' = 'Baseline S Incidence ~ 2%')

df2 <- 
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
                                     'X4' = 13.7)) %>% 
  mutate('kappa' = 'Baseline S Incidence ~ 20%')

df3 <- 
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
                                     'X4' = 13.7)) %>% 
  mutate('kappa' = 'Baseline S Incidence ~ 50%')


bind_rows(df1, df2, df3) %>% 
  pivot_longer(cols = starts_with('X'),
               names_to = 'covariate',
               values_to = 'value') %>% 
  mutate('pop' = ifelse(S == 1, 'Trial', 'Non-Trial')) %>% 
  ggplot(aes(x = value, y = kappa)) + 
  facet_wrap(~covariate) +
  geom_vline(xintercept = 1, lwd = 1.2, lty = 2, col = 'grey') + 
  geom_density_ridges(aes(fill = pop),
                      scale = 1,
                      quantile_lines = T,
                      panel_scaling = F,
                      quantiles = 0.5, 
                      alpha = 0.4, 
                      rel_min_height = 0.01) + 
  theme_bw() +
  scale_x_continuous(limits = c(-2.5, 4.5)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size = 10),
        legend.position = "bottom") + 
  labs(x = 'Covariate Value',
       y = 'Baseline Trial S Incidence',
       fill = 'Population',
       color = 'Population',
       title = 'Distribution of Simulated Covariates')

ggsave('figures/covariate_dist.png', height = 9/1.2, width = 16/1.2)




### Propensity/Outcome figure
df_preds <- 
  bind_rows(
    g_formula_preds(df1, 'LM', Y ~ A:X1 + X2 + X3 + X4) %>% mutate('kappa' = 'Baseline S Incidence ~ 2%'),
    g_formula_preds(df2, 'LM', Y ~ A:X1 + X2 + X3 + X4) %>% mutate('kappa' = 'Baseline S Incidence ~ 20%'),
    g_formula_preds(df3, 'LM', Y ~ A:X1 + X2 + X3 + X4) %>% mutate('kappa' = 'Baseline S Incidence ~ 50%'),
    g_formula_preds(df1, 'RF', Y ~ X1 + X2 + X3 + X4 + AX1 + AX2 + AX3 + AX4) %>% mutate('kappa' = 'Baseline S Incidence ~ 2%'),
    g_formula_preds(df2, 'RF', Y ~ X1 + X2 + X3 + X4 + AX1 + AX2 + AX3 + AX4) %>% mutate('kappa' = 'Baseline S Incidence ~ 20%'),
    g_formula_preds(df3, 'RF', Y ~ X1 + X2 + X3 + X4 + AX1 + AX2 + AX3 + AX4) %>% mutate('kappa' = 'Baseline S Incidence ~ 50%'),
    df1 %>% mutate('mu1' = ifelse(A == 1, Y, Y + 27.4*X1),
                   'mu0' = ifelse(A == 0, Y, Y - 27.4*X1),
                   'method' = 'True Distibution w/in S = 1'),
    df2 %>% mutate('mu1' = ifelse(A == 1, Y, Y + 27.4*X1),
                   'mu0' = ifelse(A == 0, Y, Y - 27.4*X1),
                   'method' = 'True Distibution w/in S = 1'),
    df3 %>% mutate('mu1' = ifelse(A == 1, Y, Y + 27.4*X1),
                   'mu0' = ifelse(A == 0, Y, Y - 27.4*X1),
                   'method' = 'True Distibution w/in S = 1'),
    
    df1 %>% filter(S == 0) %>%  
      mutate('mu1' = -100 + 27.4*X1 + 13.7 * (X1 + X2 + X3),
             'mu0' = -100 + 13.7 * (X1 + X2 + X3),
             'method' = 'True Distibution w/in S = 0'),
    df2 %>% filter(S == 0) %>% 
      mutate('mu1' = -100 + 27.4*X1 + 13.7 * (X1 + X2 + X3),
             'mu0' = -100 + 13.7 * (X1 + X2 + X3),
             'method' = 'True Distibution w/in S = 0'),
    df3 %>% filter(S == 0) %>% 
      mutate('mu1' = -100 + 27.4*X1 + 13.7 * (X1 + X2 + X3),
             'mu0' = -100 + 13.7 * (X1 + X2 + X3),
             'method' = 'True Distibution w/in S = 0')
    
    
  )

df_preds %>% 
  pivot_longer(cols = contains('mu'),
               names_to = 'mu', 
               values_to = 'mu_hat') %>% 
  ggplot(aes(x = mu_hat)) + 
  facet_grid(kappa~mu) + 
  geom_density(aes(fill = method), alpha = 0.5, rel_min_height = 0.01, scale = 0.7) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size = 10),
        legend.position = "bottom") + 
  labs(x = 'Predicted Value',
       y = 'Density',
       fill = 'Outcome Model',
       title = 'Distribution of Outcome Regression Predictions')

ggsave('figures/outcome_dist.png', height = 12/1.2, width = 16/1.2)