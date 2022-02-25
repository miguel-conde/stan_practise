

compute_lppd <- function(ll) {
  
  # ll: log-likelihood matrix from a stan model
  # log-pointwise-predictive-density:
  #
  #        sum_i {log(Pr_i)}
  # Pr_i = Mean likelihood of each observation in the training sample
  
  # dfmean <- ll %>% 
  #   as_tibble() %>% 
  #   exp() %>%               # Get likelihoods
  #   summarise_all(mean) %>% # Mean likelihood of each observation in the training sample
  #   gather(key, lppd_i) %>% 
  #   select(lppd_i) %>%
  #   log()
  
  mu <- mean(ll)
  
  lppd_i <- mu + log(colMeans(exp(as_tibble(ll - mu))))
  sd_i <- ll %>% apply(2, sd)
  
  lppd <-  lppd_i %>% sum()
  
  out <- list(dfmean = tibble(lppd_i, sd_i),
              lppd = lppd,
              lppd_se = sqrt(var(lppd_i) * length(lppd_i)))
  
  return(out)
}

compute_pWAIC <- function(ll) {
  # ll: log-likelihood matrix from a stan model
  #
  # Effective number of parameters:
  #
  #        sum_i {V_i}
  # V_i = variance in log-likelihood for observation in the training sample
  
  dfvar <- ll %>% 
    as_tibble() %>% 
    summarise_all(var) %>% # variance in log-likelihood for observation in the training sample
    gather(key, V_i) %>% 
    select(V_i) 
  
  pWAIC = dfvar %>% 
    sum()
  
  out <- list(dfvar = dfvar,
              pWAIC = pWAIC,
              pWAIC_se = sqrt(var(dfvar$V_i) * nrow(dfvar)))
  
  return(out)
}

compute_WAIC <- function(ll) {
  
  lppd <- compute_lppd(ll)
  pWAIC  <- compute_pWAIC(ll)
  
  WAIC <- -2 * (lppd$lppd - pWAIC$pWAIC)
  
  pointwise <- lppd$dfmean %>%
    mutate(waic_vec   = -2 * (lppd$dfmean$lppd_i - pWAIC$dfvar$V_i))
  
  waic_se <- pointwise %>%
    summarise(waic_se = (var(waic_vec) * nrow(lppd$dfmean)) %>% sqrt())
  
  out <- list(WAIC = WAIC,
              WAIC_se = as.numeric(waic_se),
              lppd = lppd$lppd,
              lppd_se = lppd$lppd_se,
              pWAIC = pWAIC$pWAIC,
              pWAIC_se = pWAIC$pWAIC_se,
              pointwise = pointwise)
  
  return(out)
}