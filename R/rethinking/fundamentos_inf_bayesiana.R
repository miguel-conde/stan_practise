library(tidyverse)

library(rstan)

library(rethinking)


## Vamos a hacer un ejercicio sencillo de inferencia bayesiana paso a paso,
## desde la generación de datos fake a la creación del modelo, primero "a mano"
## y lugo con rstan y rethinking


# 1 - MODELO 1 ------------------------------------------------------------


# 1.1 - FAKE DATA ---------------------------------------------------------

## Tenemos 10 bolas (7 azules y 3 verdes) en 1 bolsa y sacamos 9 con reposición. 
## Vamos a hacer que
## tenemos .
## Nuestros datos fake serán una extracción aleatoria:

N_azul <- 7
N_verde <- 3
N_bolas <- N_azul + N_verde

p <- N_azul / N_bolas

n_obs <- 1
n_bernouillis_x_obs <- 9

# n_obs <- número de observaciones
# size <- número de experimentos de Bernouilli en cada observación
# p <- probabilidad de éxito de cada experimento de Bernouilli
# dato_num_azules <- es el dato: número de éxitos en toda la secuencia de 
# experimentos de Bernouilli
dato_num_azules <- rbinom(n = n_obs, size = n_bernouillis_x_obs, prob = p) 

## Evidentemente, haremos como que no conocemos como se han generado los datos

# 1.2 - Modelo ------------------------------------------------------------

# Nos dicen que hab generado 1 observación consistente en n_bernouillis_x_obs
# extracciones de bolas de una bolsa en la que hay bolas azules y verdes.
# Nos piden que estimemos la proporción de bolas azules sabiendo que en esa 
# observación han salido num_azules.
#
# Partimos de nuestros a prioris acerca de como se han generado los datos: tiene
# que ser:
#
# - Likelihood: será una binomial de parámetro p
# - Prior del parámetro: uniforme de 0 a 1
#

# MÉTODO GRID
ps_grid <- seq(0, 1, length.out = 100) # La grid de p's

prior <- dunif(ps_grid, 0, 1)
likelihood <- dbinom(x = dato_num_azules, 
                     size = n_bernouillis_x_obs, 
                     prob = ps_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

plot(ps_grid, posterior, type = "o")
ps_grid[which.max(posterior)]

# Método HMCMH - rstan

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


stan_code <- "
data {
  int<lower=0> N;
  int y;
}

parameters {
  real<lower=0, upper=1> p;
}

model {
  p ~ uniform(0, 1); 
  y ~ binomial(N, p);
}
"

fit <- stan(model_code = stan_code, 
            iter = 1000, chains = 1, 
            data = list(N = n_bernouillis_x_obs, y = dato_num_azules),
            verbose = TRUE)

fit
precis(fit)

stan_posterior <- rstan::extract(fit)$p

hist(stan_posterior)


# 2 - MODELO 2 ------------------------------------------------------------

n_obs2 <- 100
n_bernouillis_x_obs2 <- 9

dato_num_azules2 <- rbinom(n = n_obs2, size = n_bernouillis_x_obs2, prob = p) 

## MÉTODO GRID
ps_grid2 <- seq(0, 1, length.out = 1000) # La grid de p's

posterior2 <- sapply(ps_grid2,
               function(p) {
                 log_likelihood <- sum(dbinom(x = dato_num_azules2,
                                              size = n_bernouillis_x_obs2,
                                              prob = p,
                                              log = TRUE))
                 likelihood <- exp(log_likelihood)
                 post <- likelihood * dunif(p, 0, 1)
                 post
               })
posterior2 <- posterior2 / sum(posterior2)

plot(ps_grid2, posterior2, type = "o")
ps_grid2[which.max(posterior2)]

# Método HMCMH - rstan

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


## MÉTODO HMCMC - stan

stan_code2 <- "
data {
  int<lower=0> N;
  int<lower=0> N_obs;
  int y[N_obs];
}

parameters {
  real<lower=0, upper=1> p;
}

model {
  p ~ uniform(0, 1); 
  y ~ binomial(N, p);
}
"

fit2 <- stan(model_code = stan_code2, 
            iter = 1000, chains = 1, 
            data = list(N = n_bernouillis_x_obs2, 
                        N_obs = n_obs2,
                        y = dato_num_azules2),
            verbose = TRUE)

fit2
precis(fit2)

stan_posterior2 <- rstan::extract(fit2)$p

hist(stan_posterior2)

# 3 - MODELO 3 ------------------------------------------------------------

# Con 2 parámetros ahora

## MÉTODO GRID
ps <- seq(0, 1, length.out = 100) # La grid de p's
Ns <- 0:30
pars_grid <- expand_grid(ps, Ns)

posterior3 <- sapply(1:nrow(pars_grid),
                    function(idx) {
                      log_likelihood <- sum( dbinom(x = dato_num_azules2,
                                                    size = pars_grid$Ns[idx],
                                                    prob = pars_grid$ps[idx],
                                                    log = TRUE))
                      likelihood <- exp(log_likelihood)
                      post <- likelihood * dunif(p, 0, 1)
                      post
                    })
posterior3 <- posterior3 / sum(posterior)

plot(pars_grid$ps, posterior3, type = "o")
plot(pars_grid$Ns, posterior3, type = "o")
pars_grid[which.max(posterior3),]


# Método HMCMH - rstan

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


## MÉTODO HMCMC - stan



stan_code3 <- "
data {
  int<lower=0> N_obs;
  int y[N_obs];
}

parameters {
  real<lower=0, upper=1> p;
  real<lower = 0> N;
}

transformed parameters {
int<lower=0> int_N;

int_N = round(N);
}

model {
  p ~ uniform(0, 1); 
  N ~ uniform(1, 20);
  y ~ binomial(int_N, p);
}
"

fit3 <- stan(model_code = stan_code3, 
             iter = 1000, chains = 1, 
             data = list(N_obs = n_obs,
                         y = dato_num_azules),
             verbose = TRUE)

fit3
precis(fit3)

stan_posterior3 <- rstan::extract(fit3)$p

hist(stan_posterior3)



# 4 - SAMPLING ------------------------------------------------------------

# 4.1 - Muestreo de la posterior ------------------------------------------

## Del modelo 2
## 10000 muestras de la posterior:
post_samples_2 <- sample(x = ps_grid2, 
                         size = 10000, 
                         prob = posterior2,
                         replace = TRUE)
plot(post_samples_2)
density(post_samples_2) %>% plot()

# 4.2 - Aplicaciones del muestreo de la posterior -------------------------


# 4.2.1 - Probabilidad de que un parámetro esté acotado entre límites -----

# Probabilidad posterior de que p < 0.75
sum(posterior2[ps_grid2 < 0.75]) # Exacta, casi nunca la podremos calcular
mean(post_samples_2 < 0.75)      # De la posterior, casi igual
mean(rstan::extract(fit2)$p < 0.75)

# Probabilidad posterior de que p < 0.7 y p > 0.6
sum(posterior2[ps_grid2 < 0.7 & ps_grid2 > 0.6]) # Exacta, casi nunca la podremos calcular
mean(post_samples_2 < 0.7 & post_samples_2 > 0.6)      # De la posterior, casi igual
mean(rstan::extract(fit2)$p < 0.7 & 
       rstan::extract(fit2)$p > 0.6)


# 4.2.2 - Intervalos de masa definida - Intervalos de credibilidad --------

# Entre qué 2 valores de p se encuentra el 89% de la probabilidad
quantile(post_samples_2, probs = c((1-.89)/2, 1-(1-.89)/2))
quantile(rstan::extract(fit2)$p, probs = c((1-.89)/2, 1-(1-.89)/2))

# Este tipo de intervalos asignan la misma masa de probailiadd a cada cola.
# Se llaman "PERCENTILE INTERVALS", PI.
quantile(post_samples_2, probs = c(0.1, 1 - 0.1))
quantile(post_samples_2, probs = c(0.2, 1 - 0.2))
quantile(post_samples_2, probs = c(0.3, 1 - 0.3))
quantile(post_samples_2, probs = c(0.25, 1 - 0.25))

# Funcionan bien si la posterior es bastante simétrica. Si está bastante sesgada,
# mjor usar HIGHEST POSTERIOR DENSITY INTERVALS (HDPI) = el intervalo más estrecho
# que contiene la masa de probabilidad especificada
rethinking::HPDI(post_samples_2, prob = 0.5)


# 4.2.3 - Estimaciones puntuales ------------------------------------------

## Sirven como resúmenes de la posterior

## 1 - Estimación MAP (Maximum A Posteriori) - MODA
ps_grid2[which.max(posterior2)] # Exacta, no siempre podremos calcularla
rethinking::chainmode(post_samples_2, adj = 0.01)
rethinking::chainmode(rstan::extract(fit2)$p, adj = 0.01)

## 2 - Media
sum(ps_grid2 * posterior2)
mean(post_samples_2)
mean(rstan::extract(fit2)$p)

## 3 - Mediana
quantile(post_samples_2, prob = .5)
  quantile(rstan::extract(fit2)$p, prob = .5)

## 4 - Function loss

# Minimizamos una medida de la incertidumbre sobre el valor real

## L1 - Coincide con la MEDIANA
loss <-sapply(ps_grid2, function(d) sum(posterior2 * abs(d - ps_grid2)))
ps_grid2[which.min(loss)]

## L2 - Coincide con la media
loss <-sapply(ps_grid2, function(d) sum(posterior2 * (d - ps_grid2)^2))
ps_grid2[which.min(loss)]


# 4.3 - Muestreo para simular predicciones --------------------------------


# 4.3.1 - Dummy data ------------------------------------------------------

# Es lo que hemos hecho para los datos de los modelos de arriba

# Comprobación:
# Probabilidades de cada caso:
dbinom(0:n_bernouillis_x_obs2, size = n_bernouillis_x_obs2, prob = p)

# Simulamos dummy data
dummy_data <- rbinom(n = 1e6, size = n_bernouillis_x_obs2, prob = p)
table(dummy_data) / 1e6


# 4.3.2 - Model Checking --------------------------------------------------


# 4.3.2.1 - Posterior Predictive Distribution -----------------------------

## PPS (Posterior Predictive Sampling)
w <- rbinom(1e5, size = n_bernouillis_x_obs2, prob = post_samples_2)

# Con Stan
stan_code2_pps <- "
data {
  int<lower=0> N;
  int<lower=0> N_obs;
  int y[N_obs];
}

parameters {
  real<lower=0, upper=1> p;
}

model {
  p ~ uniform(0, 1); 
  y ~ binomial(N, p);
}

generated quantities {
  int y_pps[N_obs];
  
  for (n in 1:N_obs) {
    y_pps[n] = binomial_rng(N, p);
  }
    
}
"

fit2_pps <- stan(model_code = stan_code2_pps, 
             iter = 1000, chains = 1, 
             data = list(N = n_bernouillis_x_obs2, 
                         N_obs = n_obs2,
                         y = dato_num_azules2),
             verbose = TRUE)

# Salen 500 (1 por cada post-warmup draws) distribuciones de cada simulación predictiva

extract(fit2_pps)$y_pps %>% apply(1, function(x) table(x)/100) %>% bind_rows()


# EJERCICIOS --------------------------------------------------------------


# Easy --------------------------------------------------------------------

# The Easy problems use the samples from the posterior distribution for 
# the globe tossing example.
# This code will give you a specific set of samples, so that you can check your 
# answers exactly.
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size = 9, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

set.seed(100)
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

# Use the values in samples to answer the questions that follow.

plot(p_grid, posterior, type = "l")
hist(samples, freq = FALSE)
lines(density(samples), col = "red")


# 3E1 ---------------------------------------------------------------------

# How much posterior probability lies below p = 0.2?

# Respuesta exacta de la posterior distribution:
sum(posterior[p_grid < 0.2])

# Respuesta de la posterior sample:
mean(samples < 0.2) 


# 3E2 ---------------------------------------------------------------------

# How much posterior probability lies above p = 0.8?

# Respuesta exacta de la posterior distribution:
sum(posterior[p_grid > 0.8])

# Respuesta de la posterior sample:
mean(samples > 0.8) 

# 3E3 ---------------------------------------------------------------------

# How much posterior probability lies between p = 0.2 and p = 0.8?

# Respuesta exacta de la posterior distribution:
sum(posterior[p_grid < 0.8]) - sum(posterior[p_grid < 0.2])
sum(posterior[p_grid > 0.2 & p_grid < 0.8])

# Respuesta de la posterior sample:
mean(samples > 0.2 & samples < 0.8) 

# 3E4 ---------------------------------------------------------------------

# 20% of the posterior probability lies below which value of p?

# Respuesta de la posterior sample:
quantile(samples, 0.2)

# Comprobación con la posterior:
sum(posterior[p_grid < quantile(samples, 0.2)])

# Y con la mustra d ela posterior
sum(samples < quantile(samples, 0.2)) / 1e4

# 3E5 ---------------------------------------------------------------------

# 20% of the posterior probability lies above which value of p?

# Respuesta de la posterior sample:
quantile(samples, 0.8,)

# 3E6 ---------------------------------------------------------------------

# Which values of p contain the narrowest interval equal to 66 % of the 
# posterior probability?
rethinking::HPDI(samples, prob = 0.66) # TODO

# 3E7 ---------------------------------------------------------------------

# Which values of p contain 66 % of the posterior probability, assuming equal 
# posterior probability both below and above the interval?
rethinking::PI(samples, prob = 0.66) # TODO
quantile(samples, c((1-.66)/2, 1 - (1-.66)/2))

# Medium ------------------------------------------------------------------


# 3M1 ---------------------------------------------------------------------

# Suppose the globe tossing data had turned out to be 8 water in 15 tosses.
# Construct the posterior distribution, using grid approximation. Use the same 
# flat prior as before.

p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(8, size = 15, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

plot(p_grid, posterior, type = "l")

# 3M2 ---------------------------------------------------------------------

# Draw 10,000 samples from the grid approximation from above.Then use the 
# samples to calculate the 90 % HPDI for p.

set.seed(100)
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)


hist(samples, freq = FALSE)
lines(density(samples), col = "red")

rethinking::HPDI(samples, prob = 0.9) # TODO

# 3M3 ---------------------------------------------------------------------

# Construct a posterior predictive check for this model and data.This means 
# simulate the distribution of samples, averaging over the posterior 
# uncertainty in p. What is the probability of observing 8 water in 15 tosses?

pps <- rbinom(1e4, size = 15, prob = samples)

mean(pps == 8)

# PAra cada posible resultado, de la posterior simulation se obtienen estas
# probabilidades:
table(pps) / 1e4

# Comprobemos: la probabilidad de la posterior es:
p <- p_grid[which.max(posterior)]

# Luego las probabilidades de los diferentes resultados son:
dbinom(0:15, 15, p)

# Comparemos gráficamente:
plot(as.vector(table(pps) / 1e4), dbinom(0:15, 15, p),
     xlab = "PPS", ylab = "POSTERIOR")
abline(a = 0, b = 1, lty = 2)

# O:
plot(0:15, as.vector(table(pps) / 1e4), type = "o",
     ylim = c(0, 0.20),
     xlab = "Number of Ws", ylab = "Probability")
lines(0:15, dbinom(0:15, 15, p), col = "red", type = "o")
legend("topleft", legend = c("Posterior", "PPS"), 
       # type = "o",
       lty = 1,
       col = c(1, "red"),
       bty = "n")

# 3M4 ---------------------------------------------------------------------

# Using the posterior distribution constructed from the new (8/15) data, now 
# calculate the probability of observing 6 water in 9 tosses.
pps <- rbinom(1e4, size = 9, prob = samples)

mean(pps == 6)

# 3M5 ---------------------------------------------------------------------

# Start over at 3M1, but now use a prior that is zero below p = 0.5 and a 
# constant above p = 0.5. This corresponds to prior information that a majority 
# of the Earth’s surface is water. Repeat each problem above and compare the 
# inferences (using both priors) to the true value p = 0.7.

p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- c(rep(0, length(p_grid[p_grid < 0.5])),
           rep(2, length(p_grid[p_grid >= 0.5])))
likelihood <- dbinom(8, size = 15, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

plot(p_grid, posterior, type = "l")

# 3M5.2
set.seed(100)
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)


hist(samples, freq = FALSE)
lines(density(samples), col = "red")

rethinking::HPDI(samples, prob = 0.9) # TODO

# 3M5.3
pps <- rbinom(1e4, size = 15, prob = samples)

mean(pps == 8)

# 3M5.4
pps <- rbinom(1e4, size = 9, prob = samples)

mean(pps == 6)

# 3M6 ---------------------------------------------------------------------

# Suppose you want to estimate the Earth’s proportion of water very precisely. 
# Specifically, you want the 99% percentile interval of the posterior 
# distribution of p to be only 0.05 wide. This means the distance between the 
# upper and lower bound of the interval should be 0.05. How many times will you 
# have to toss the globe to do this?

xxx <- function(n_tosses, p = 0.7, p_i = .99) {
  
  p_grid <- seq(from = 0, to = 1, length.out = 1000)
  prior <- dunif(p_grid, 0, 1)
  likelihood <- dbinom(n_tosses * p, size = n_tosses, prob = p_grid)
  posterior <- likelihood * prior
  posterior <- posterior / sum(posterior)
  
  samples <- sample(p_grid, size = 1e4, prob = posterior, replace = TRUE)
  
  out <- quantile(samples, c((1-p_i)/2, 1 - (1-p_i)/2))
  
  return(as.numeric(diff(out)))
}

xxx(100)
xxx(1000)
xxx(1500)
xxx(1600)
xxx(1800)
xxx(1900)
xxx(2200)

# Hard --------------------------------------------------------------------

# The Hard problems here all use the data below. These data indicate the gender 
# (male = 1, female = 0) of officially reported first and second born children 
# in 100 two-child families. So for example, the first family in the data 
# reported a boy (1) and then a girl (0). The second family reported a girl (0) 
# and then a boy (1). The third family reported two girls. You can load these 
# tow vectors into R’s memory by typing:

data(homeworkch3)

birth1 
birth2 

# Use these vectors as data. So for example to compute the total number of boys 
# born across all of these births, you could use:

sum(birth1) + sum(birth2)


# 3H1 ---------------------------------------------------------------------

# Using grid approximation, compute the posterior distribution for the 
# probability of a birth being a boy. Assume a uniform prior probability. Which 
# parameter value maximizes the posterior probability?

n_boys <- sum(birth1) + sum(birth2)
n_births <- length(birth1) + length(birth2)

p_grid <- seq(0, 1, length.out = 1000)

prior <- dunif(p_grid)
likelihood <- dbinom(x = n_boys, size = n_births, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

plot(p_grid, posterior, type = "l",
     xlab = "p",
     ylab = "p Posterior Distribution")

p_grid[which.max(posterior)]

# 3H2 ---------------------------------------------------------------------

# Using the sample function, draw 10,000 random parameter values from the 
# posterior distribution you calculated above. Use these sample to estimate the 
# 50%, 89%, and 97% highest posterior density intervals.

samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

rethinking::HPDI(samples, c(.5, .89, .97))

# 3H3 ---------------------------------------------------------------------

# Use rbinom to simulate 10,000 replicates of 200 births. You should end up with 
# 10,000 numbers, each one a count of boys out of 200 births. Compare the 
# distribution of predicted numbers of boys to the actual count in the data 
# (111 boys out of 200 births). There are many good ways to visualize the 
# simulations, but the dens command (part of the rethinking package) is probably 
# the easiest way in this case. Does it look like the model fits the data well? 
# That is, does the distribution of predictions include the actual observation 
# as a central, likely outcome?

pps <- rbinom(1e4, size = n_births, prob = samples)

mean(pps)

rethinking::dens(pps, show.HPDI = .89)
abline(v = n_boys, lty = 2, col = "red")
abline(v = mean(pps), lty = 2)

# Based on this posterior predictive distribution, the model appears to fit well, 
# with the observed value of 111 in the middle of the distribution.

# 3H4 ---------------------------------------------------------------------

# Now compare 10,000 counts of boys from 100 simulated first borns only the 
# number of boys in the first births, birth1. How does the model look in this 
# light?

pps <- rbinom(1e4, size = length(birth1), prob = samples)

mean(pps)

rethinking::dens(pps, show.HPDI = .89)
abline(v = sum(birth1), lty = 2, col = "red")
abline(v = mean(pps), lty = 2)

# Based on first births only, the model appears to fit less well. Specifically, 
# the model appears to be overestimating the number of first births who are boys. 
# However, it does not appear to be a large discrepancy, as the observed value 
# is still within the middle 66% interval.

# 3H5 ---------------------------------------------------------------------

# The model assumes that sex of first and second births are independent. To 
# check this assumption, focus now on second births that followed female first 
# borns. Compare 10,000 simulated conts of boys to only those second births that 
# followed girls. To do this correctly, you need to count the number of first 
# borns who were girls and simulate that many births, 10,000 times. Compare the 
# counts of boys in your simulations to the actual observed count of boys 
# following girls. How does the model look in this light? Any guesses what is 
# going on in these data?
pps <- rbinom(1e4, size = sum(birth1 == 0), prob = samples)

mean(pps)

obs <- sum(birth2[birth1 == 0])

rethinking::dens(pps, show.HPDI = .89)
abline(v = obs, lty = 2, col = "red")
abline(v = mean(pps), lty = 2)

# The model is severely under estimating the number of second-born boys when the 
# first born child is a girl. Thus, our assumption that births are independent 
# appears to be violated.
