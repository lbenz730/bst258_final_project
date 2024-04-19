### Collection of useful helper functions

### Inverse Logit Function
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

### Logit Function
logit <- function(x) {
  return(log(x/(1-x)))
}

### Function to compute the fitted value of a model given a list of coefficients
compute_model <- function(df, beta) {
  ### Initialize fitted values to be 0
  preds <- rep(0, nrow(df))
  
  ### Iterate over coefficients
  for(i in 1:length(beta)) {
    coeff <- names(beta)[i]
    
    if(coeff == '(Intercept)') { ### Intercept
      preds <- preds + beta[[coeff]]
    } else if(grepl('\\[', coeff) & grepl(':', coeff)) { ### Level of categorical variable w/ interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeffs[1]))
      variable <- gsub('\\[.*$', '', coeffs[1])
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + df[[ coeffs[2] ]][ix] * beta[[ coeff ]]
    } else if(grepl('\\[', coeff)) { ### Level of categorical variable
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeff))
      variable <- gsub('\\[.*$', '', coeff)
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + beta[[coeff]]
    } else if(grepl(':', coeff)) { ### Interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      if(length(coeffs) == 2) {
        preds <- preds + df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * beta[[coeff]]
      } else {
        preds <- preds + df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * df[[ coeffs[3] ]] *  beta[[coeff]]
      }
    } else if(grepl('\\^2\\)$', coeff)) { ### quadratic 
      variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
      preds <- preds + df[[ variable ]]^2 * beta[[coeff]]
    } else if(grepl('^I\\(exp\\(', coeff)) { ### exponential
      variable <- gsub('\\)\\)$', '', gsub('^I\\(exp\\(', '', coeff))
      preds <- preds + exp(df[[ variable ]]) * beta[[coeff]]
    } else { ### Linear Term 
      preds <- preds + df[[coeff]] * beta[[coeff]]
    }
  }
  
  return(preds)
}


### Weight Truncation Function
truncate_weights <- function(w, q1, q2) {
  lower_bound <- quantile(w, q1)
  upper_bound <- quantile(w, q2)
  w[w <= lower_bound] <- lower_bound
  w[w >= upper_bound] <- upper_bound
  return(w)
}

interact <- function(df) {
  df  <- 
    df %>% 
    mutate('AX1' = A * X1,
           'AX2' = A * X2,
           'AX3' = A * X3,
           'AX4' = A * X4) 
  
  return(df)
}

theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 12),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))
