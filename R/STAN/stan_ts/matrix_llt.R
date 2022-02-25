library(rstan)
library(bayesplot)
library(tidybayes)
library(tidyverse)

sales_ts <- readRDS(here::here("STAN", "data", "sales_ts.Rds")) %>% 
  gather(ctry_segment, value, -date) %>% 
  separate(ctry_segment, into = c("country", "segment")) %>% 
  mutate_at(vars(country, segment), ~ factor(.))

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

y <- sales_ts %>% 
  filter(country == "IT", segment == "existing") %>% 
  pull(value)

fit_matrix_llt <- stan(here::here("STAN", "stan_ts", "matrix_llt.stan"), 
                     data = list(y = y,
                                  N = length(y)),
                     iter = 2000, 
                     chains = 4)

plot(y, type = "l")

rstan::extract(fit_matrix_llt, "x")[[1]] %>% 
  .[,,1] %>% 
  colMeans() %>% 
  lines(col = "blue") 

as.data.frame(fit_matrix_llt) %>% 
  select(starts_with("mu_y[")) %>% 
  colMeans() %>% 
  lines(col = "blue")
