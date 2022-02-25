library(rstan)
library(bayesplot)
library(tidybayes)
library(tidyverse)

sales_ts <- readRDS(here::here("STAN", "data", "sales_ts.Rds"))

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

# LLT ---------------------------------------------------------------------

Y = sales_ts$IT_existing %>% scale() %>% as.numeric()

N = length(Y)

fit_llt_1 <- stan(here::here("STAN", "stan_ts", "llt_1.stan"), 
                data = list(N = N, X = Y),
                iter = 2000, 
                chains = 4)

plot(Y, type = "l")

as.data.frame(fit_llt_1) %>% 
  select(starts_with("u")) %>% 
  colMeans() %>% 
  lines(col = "blue")

fit_llt_2 <- stan(here::here("STAN", "stan_ts", "llt_2.stan"), 
                  data = list(T = N, y = Y),
                  iter = 2000, 
                  chains = 4)

as.data.frame(fit_llt_2) %>% 
  select(starts_with("u")) %>% 
  colMeans() %>% 
  lines(col = "red", lwd = 2)






fit_llt<- stan(here::here("STAN", "stan_ts", "llt.stan"), 
                  data = list(N = N, y = Y),
                  iter = 2000, 
                  chains = 4)

as.data.frame(fit_llt) %>% 
  select(starts_with("level")) %>% 
  colMeans() %>% 
  lines(col = "green")





library(loo)

logLikelihood_1 <- extract_log_lik(fit_llt_1, "logLikelihood")
WAIC_1 <- waic(logLikelihood_1)

logLikelihood_2 <- extract_log_lik(fit_llt_1, "logLikelihood")
WAIC_2 <- waic(logLikelihood_2)

logLikelihood_mine <- extract_log_lik(fit_llt, "logLikelihood")
WAIC_mine <- waic(logLikelihood_mine)

compare(WAIC_1, 
        WAIC_2, 
        WAIC_mine)
loo_compare(WAIC_1, 
            WAIC_2, 
            WAIC_mine)

LOO_1 <- loo(logLikelihood_1)
LOO_2 <- loo(logLikelihood_2)
LOO_mine <- loo(logLikelihood_mine)

compare(LOO_1, 
        LOO_2, 
        LOO_mine)
loo_compare(LOO_1, 
            LOO_2, 
            LOO_mine)
