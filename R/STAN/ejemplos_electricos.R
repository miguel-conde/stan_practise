library(tidyverse)

library(readr)

library(units)

library(lubridate)

# GENERACIÓN (MWh) PROGRAMADA PBF SOLAR FOTOVOLTAICA
# DESDE EL 01-09-2010 A LAS 00:00 HASTA EL 31-08-2021 A LAS 23:50 AGRUPADOS POR HORA
# https://www.esios.ree.es/es/analisis/14?vis=1&start_date=01-09-2010T00%3A00&end_date=31-08-2021T23%3A50&compare_start_date=31-08-2010T00%3A00&groupby=hour

historico_generacion_fotovoltaica <- 
  read_delim("data/historico_generacion_fotovoltaica.csv", 
             ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(datetime = datetime + hours(1)) %>% 
  mutate(value = as_units(value, "megawatthour"))

# GENERACIÓN (MWh) PROGRAMADA PBF EOLICA
# DESDE EL 01-09-2010 A LAS 00:00 HASTA EL 31-08-2021 A LAS 23:50 AGRUPADOS POR HORA
# https://www.esios.ree.es/es/analisis/10073?vis=1&start_date=01-09-2010T00%3A00&end_date=31-08-2021T23%3A50&compare_start_date=31-08-2010T00%3A00&groupby=hour
historico_generacion_eolica <- 
  read_delim("data/historico_generacion_eolica.csv", 
             ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(datetime = datetime + hours(1)) %>% 
  mutate(value = as_units(value, "megawatthour"))

# POTENCIAS INSTALADAS (MW) por TECNOLOGÍA | Del 09/2019 al 08/2021
# https://www.ree.es/es/datos/generacion/potencia-instalada
historico_potencias_instaladas <-
  read_csv("data/historico_potencias_instaladas.csv", 
           locale = locale(decimal_mark = ",", grouping_mark = ".", 
                           encoding = "WINDOWS-1252"), skip = 4) %>% 
  drop_na() %>% 
  rename(tecnologia = `...1`) %>% 
  mutate_at(vars(2), ~ as.numeric(.)) %>% 
  gather(date, potencia, -tecnologia) %>% 
  spread(tecnologia, potencia) %>% 
  janitor::clean_names() %>% 
  mutate(date = lubridate::dmy(paste0("01/", date))) %>% 
  mutate_if(is.numeric, ~ as_units(., "megawatt")) %>% 
  mutate(datetime = as.POSIXct(date))






res <- historico_generacion_fotovoltaica %>% 
  select(datetime, generation_pv = value) %>% 
  full_join(historico_generacion_eolica %>%  
              select(datetime, generation_eolic = value),
            by = "datetime") %>% 
  full_join(historico_potencias_instaladas %>% 
              mutate(datetime = as.POSIXct(date)) %>% 
              select(datetime, 
                     installed_power_eolic = eolica, 
                     installed_power_pv = solar_fotovoltaica),
            by = "datetime") %>% 
  fill(starts_with("installed_power_")) %>% 
  drop_na() %>% 
  mutate(load_factor_pv = generation_pv / installed_power_pv,
         load_factor_eolic = generation_eolic / installed_power_eolic) %>% 
  # Gropup by month, day, hour
  group_by(mth = month(datetime), d = day(datetime), hr = hour(datetime)) %>% 
  summarise_at(vars(-datetime), ~ mean(.), .groups = "drop") %>% 
  ungroup() %>% 
  mutate(datetime = ymd_hms(paste0("2020/", mth, "/", d, " ", hr, ":00:00")),
         .before = mth)

res2 <- res %>% ungroup() %>% 
  mutate(datetime = ymd_hms(paste0("2020/", mth, "/", d, " ", hr, ":00:00"))) %>%
  select(datetime, starts_with("load"))


res2 %>% mutate_at(vars(starts_with("load_")), ~ 100 * .  ) %>% 
  Conento::hchart_lineas()

# MOCK DATA PV ------------------------------------------------------------

library(forecast)
auto_arima_res2_pv <- auto.arima(ts(res2$load_factor_pv))

synthetic_serie <- fitted(auto_arima_res2_pv) + 
  sample(residuals(auto_arima_res2_pv), 
         length(residuals(auto_arima_res2_pv)), 
         replace = TRUE) 
synthetic_serie[synthetic_serie < 0] <- 0

idx_low_load_factor <- res$load_factor_pv < 
  quantile(res$load_factor_pv, 
           probs = seq(0, 1, by = .05))["55%"]
synthetic_serie[idx_low_load_factor] <- 0
synthetic_serie <- synthetic_serie/max(synthetic_serie) * as.numeric(max(res2$load_factor_pv))

plot(synthetic_serie %>% tail(200), type = "l")
plot(synthetic_serie, type = "l")

synthetic_serie_pv <- synthetic_serie


# MOCK DATA EOLIC ---------------------------------------------------------


auto_arima_res2_eo <- auto.arima(ts(res2$load_factor_eolic))

synthetic_serie <- fitted(auto_arima_res2_eo) + 
  sample(residuals(auto_arima_res2_eo), 
         length(residuals(auto_arima_res2_eo)), 
         replace = FALSE) 
synthetic_serie[synthetic_serie < 0] <- 0

idx_low_load_factor <- res$load_factor_pv < 
  quantile(res$load_factor_pv, 
           probs = seq(0, 1, by = .05))["0%"]
synthetic_serie[idx_low_load_factor] <- 0
synthetic_serie <- synthetic_serie/max(synthetic_serie) * as.numeric(max(res2$load_factor_eolic))

plot(synthetic_serie %>% tail(200), type = "l")
plot(synthetic_serie, type = "l")

synthetic_serie_eolic <- synthetic_serie



# MOCK DATA MANAGED RENWABLES ---------------------------------------------



probe <- res %>% 
  mutate(installed_powermngd_renwables = installed_power_eolic + 
           installed_power_pv)

probe <- probe %>% mutate(generation_power_mngd_renwables = 
                            generation_eolic + generation_pv)
probe <- probe %>% mutate(load_factor_mngd_renwables = 
                            generation_power_mngd_renwables / installed_powermngd_renwables )

probe <- probe %>% mutate(datetime = ymd_hms(paste0("2020/", mth, "/", d, " ", hr, ":00:00")),
                          .before = mth)
probe_2 <- probe %>% select(datetime, starts_with("load"))

probe_2 <- probe %>% ungroup() %>% select(datetime, starts_with("load"))

probe_2

probe_2 %>% mutate_at(vars(starts_with("load_")), ~ 100 * .  ) %>% 
  Conento::hchart_lineas()

probe_2$load_factor_mngd_renwables %>% mean()

probe_2$load_factor_mngd_renwables %>% hist()



auto_arima_res2_renew <- auto.arima(ts(probe_2$load_factor_mngd_renwables))

synthetic_serie <- fitted(auto_arima_res2_renew) + 
  sample(residuals(auto_arima_res2_renew), 
         length(residuals(auto_arima_res2_renew)), 
         replace = FALSE) 
synthetic_serie[synthetic_serie < 0] <- 0

idx_low_load_factor <- res$load_factor_pv < 
  quantile(res$load_factor_pv, 
           probs = seq(0, 1, by = .05))["0%"]
synthetic_serie[idx_low_load_factor] <- 0
synthetic_serie <- synthetic_serie/max(synthetic_serie) * as.numeric(max(res2$load_factor_eolic))

plot(synthetic_serie %>% tail(200), type = "l")
plot(synthetic_serie, type = "l")

synthetic_serie_renew <- synthetic_serie

# BETA - STAN -------------------------------------------------------------

library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


stan_code <- "
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  alpha ~ exponential(1);
  beta ~ exponential(1);
  y ~ beta(alpha, beta);
}
"

## DATA

Y <- as.numeric(probe_2$load_factor_mngd_renwables)

## PRIORS
curve(dexp(x, 1), 0, 10)

# Prior predictive simulation
sample_alpha <- rexp(1e4, 1/7)
sample_beta <- rexp(1e4, 1/18)
prior_ps_y <- rbeta(1e4, sample_alpha, sample_beta)

plot(density(prior_ps_y))

## MODEL

# Fit the model

fit <- stan(model_code = stan_code, 
            iter = 2000, chains = 4, 
            data = list(N = length(Y), y = Y),
            verbose = TRUE)

alpha <- extract(fit, "alpha")[[1]] 
qplot(alpha)

beta <- extract(fit, "beta")[[1]] 
qplot(beta)

pairs(fit)

library(shinystan) 

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)

# Residuals
stan_code_resid <- "
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  mu ~ normal(0, .1);
  sigma ~ exponential(.1);
  y ~ cauchy(mu, sigma);
}

generated quantities {
  vector[N] lSimData;
  int aMax_indicator;
  int aMin_indicator;
  
  // Generate posterior predictive samples
  for (i in 1:N) {
    lSimData[i] = cauchy_rng(mu, sigma);
  }
  
  // Compare with real data
  aMax_indicator = max(lSimData) > max(y);
  aMin_indicator = min(lSimData) < max(y);
}
"

Y <- as.numeric(residuals(auto_arima_res2_renew))

fit_resid <- stan(model_code = stan_code_resid, 
                  iter = 500, chains = 4, 
                  data = list(N = length(Y), y = Y),
                  verbose = TRUE)

print(fit_resid, pars = c("mu", "sigma"))

mu <- extract(fit_resid, "mu")[[1]] 
qplot(mu)

sigma <- extract(fit_resid, "sigma")[[1]] 
qplot(sigma)

pairs(fit_resid)

library(shinystan) 

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)

## PV
Y_pv <- as.numeric(probe_2$load_factor_pv)

fit_pv <- stan(model_code = stan_code, 
               iter = 2000, chains = 4, 
               data = list(N = length(Y_pv), y = Y_pv),
               verbose = TRUE) 
fit_pv

# Residuals
Y <- as.numeric(residuals(auto_arima_res2_pv))

fit_pv_resid <- stan(model_code = stan_code_resid, 
                     iter = 2000, chains = 4, 
                     data = list(N = length(Y), y = Y),
                     verbose = TRUE)

print(fit_pv_resid, pars = c("mu", "sigma"))

mu <- extract(fit_pv_resid, "mu")[[1]] 
qplot(mu)

sigma <- extract(fit_pv_resid, "sigma")[[1]] 
qplot(sigma)

pairs(fit_pv_resid)

library(shinystan) 

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)

## EOLIC

Y_eo <- as.numeric(probe_2$load_factor_eolic)

fit_eo <- stan(model_code = stan_code, 
               iter = 2000, chains = 4, 
               data = list(N = length(Y_eo), y = Y_eo),
               verbose = TRUE) 
fit_eo

library(rethinking)

precis(fit_eo)

# Posterior Predictive Sampling
rows_samples_eo <- 
  lapply(rstan::extract(fit_eo), function(x) sample(x, 1e5, replace = TRUE))
rbeta(1e5, rows_samples_eo$alpha, rows_samples_eo$beta) %>% 
  density() %>% 
  plot()
