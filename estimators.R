### Compute the ATE in the RCT via difference in means
rct_ate <- function(df) {
  df_tau <- 
    df %>% 
    filter(S == 1) %>% 
    summarise('method' = 'RCT ATE', 
              'tau_hat' = mean(Y[A == 1]) - mean(Y[A == 0]))
  
  return(df_tau)
}

### IPSW
ipsw <- function(df, model_type, model_formula) {
  ### Compute Selection Propensity
  if(model_type == 'GLM') {
    selection_model <- 
      glm(model_formula, 
          family = 'binomial',
          data = df)
    
    df$prob_s <- selection_model$fitted.values
    
  } else if(model_type == 'RF') {
    selection_model <- 
      ranger(model_formula, 
             probability = T,
             max.depth = 8,
             data = df)
    
    
    df$prob_s <- selection_model$predictions[,'1']
    
  }
  
  df_trial <- 
    df %>% 
    filter(S == 1)
  
  # ### Weight Truncation
  # if(model_type == 'Random Forest') {
  #   df_trial$prob_s <- truncate_weights(df_trial$prob_s, 0.01, 0.99)
  # }
  
  
  ### Compute Estimate
  df_tau <-
    df_trial %>% 
    summarise('method' = paste0('IPW (', model_type, ')'),
              'tau_hat' = sum(Y * (1 - prob_s)/prob_s * (A/0.5 - (1-A)/0.5))/sum(df$S == 0))
  
  return(df_tau)
}

g_formula <- function(df, model_type, model_formula) {
  df_trial <- 
    df %>% 
    filter(S == 1)
  
  df_obs <- 
    df %>% 
    filter(S == 0)
  
  if(model_type == 'LM') {
    outcome_model <- lm(model_formula, data = df_trial)
    mu0 <- predict(outcome_model, newdata = mutate(df_obs, 'A' = 0))
    mu1 <- predict(outcome_model, newdata = mutate(df_obs, 'A' = 1))
  } else if(model_type == 'RF') {
    outcome_model <- 
      ranger(model_formula, 
             data = interact(df_trial),
             max.depth = 6)
    
    
    
    mu0 <- predict(outcome_model, data = interact(mutate(df_obs, 'A' = 0)))$predictions
    mu1 <- predict(outcome_model, data = interact(mutate(df_obs, 'A' = 1)))$predictions 
    
  }
  
  df_tau <-
    tibble('method' = paste0('G-Formula (', model_type, ')'),
           'tau_hat' = mean(mu1 - mu0))
  
  return(df_tau)
}

aipsw <- function(df, selection_formula, outcome_formula, model_type) {
  
  if(model_type == 'GLM/LM') {
    selection_model <- 
      glm(selection_formula, 
          family = 'binomial',
          data = df)
    
    df$prob_s <- selection_model$fitted.values
    
  } else if(model_type == 'RF') {
    selection_model <- 
      ranger(selection_formula, 
             probability = T,
             max.depth = 8,
             data = df)
    
    
    df$prob_s <- selection_model$predictions[,'1']
    
  }
  
  df <- 
    df %>% 
    mutate('inv_odds_s' = (1 - prob_s)/prob_s) 
  
  ### Outcome Model 
  df_trial <- 
    df %>% 
    filter(S == 1)
  
  if(model_type == 'GLM/LM') {
    outcome_model <- lm(outcome_formula, data = df_trial)
    df$mu0 <- predict(outcome_model, newdata = mutate(df, 'A' = 0))
    df$mu1 <- predict(outcome_model, newdata = mutate(df, 'A' = 1))
  } else if(model_type == 'RF') {
    outcome_model <- 
      ranger(outcome_formula, 
             data = interact(df_trial), 
             max.depth = 6)
    
    
    df$mu0 <- predict(outcome_model, data = interact(mutate(df, 'A' = 0)))$predictions
    df$mu1 <- predict(outcome_model, data = interact(mutate(df, 'A' = 1)))$predictions 
  }
  
  
  
  
  tau_g <- mean(df$mu1[df$S == 0] - df$mu0[df$S == 0])
  
  tau_aug <- 
    df %>% 
    filter(S == 1) %>% 
    summarise('tau_aug' = sum(inv_odds_s * (A/0.5 * (Y-mu1) - (1-A)/0.5 * (Y-mu0)))/sum(df$S == 0)) %>% 
    pull(tau_aug)
  
  df_tau <- 
    tibble('method' = paste0('AIPSW (', model_type, ')'),
           'tau_hat' = tau_g + tau_aug )
  
  
  return(df_tau)
  
}


g_formula_preds <- function(df, model_type, model_formula) {
  df_trial <- 
    df %>% 
    filter(S == 1)
  
  df_obs <- 
    df %>% 
    filter(S == 0)
  
  if(model_type == 'LM') {
    outcome_model <- lm(model_formula, data = df_trial)
    mu0 <- predict(outcome_model, newdata = mutate(df_obs, 'A' = 0))
    mu1 <- predict(outcome_model, newdata = mutate(df_obs, 'A' = 1))
  } else if(model_type == 'RF') {
    outcome_model <- 
      ranger(model_formula, 
             data = interact(df_trial),
             max.depth = 6)
    
    
    
    mu0 <- predict(outcome_model, data = interact(mutate(df_obs, 'A' = 0)))$predictions
    mu1 <- predict(outcome_model, data = interact(mutate(df_obs, 'A' = 1)))$predictions 
    
  }
  
  df_tau <-
    tibble('method' = model_type,
           'mu0' = mu0,
           'mu1' = mu1)
  
  return(df_tau)
}