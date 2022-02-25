library(tidyverse)

library(rstan)

library(rethinking)


# DATA --------------------------------------------------------------------


data(Howell1)

d <- Howell1
d2 <- Howell1 %>% filter(age >= 18)


# HEIGHT MODEL ------------------------------------------------------------


stan_code <- "
data{
  int<lower=0> N;
  vector[N] h;
}

parameters {
  real mu;
  real<lower=0,upper=50> sigma;
}

model {
  mu ~ normal(180, 20); // 95% mu in 140 - 220
  sigma ~ uniform(0, 50); // 95 % h in mu +/- 0 - 100
  h ~ normal(mu, sigma); 
}

"

# 1.Prior exploration - Prior predictive ----------------------------------

# Explore priors
curve(dnorm(x, 180, 20), from = 100, to = 240, main = "Prior mu")
curve(dunif(x, 0, 50), from = -100, to = 100, main = "Prior sigma")

# Prior predictive simulation
prior_sample_mu <- rnorm(1e4, 180, 20)
prior_sample_sigma <- runif(1e4, 0, 50)
prior_sample_h <- rnorm(1e4, prior_sample_mu, prior_sample_sigma)

plot(density(prior_sample_h))

# What if greater sigma prior?
prior_sample_h2 <- rnorm(1e4, prior_sample_mu, runif(1e4, 0, 100))
lines(density(prior_sample_h2), col = "red")

sum(prior_sample_h2 < 0) / 1e4   # Not credible
sum(prior_sample_h2 > 272) / 1e4 # Not credible


# 2.Grid approximation ----------------------------------------------------

post_table <- expand.grid(mu = seq(150, 160, length.out = 100),
                 sigma = seq(7, 9, length.out = 100))

# Log-likelihood
post_table$LL <- sapply(1:nrow(post_table),
                        function(i) {
                          sum(
                            dnorm(d2$height, post_table$mu[i], post_table$sigma[i],
                                  log = TRUE)
                          ) 
                        }) 

post_table$log_posterior <- post_table$LL +
  dnorm(post_table$mu, 180, 20, TRUE) +
  dunif(post_table$sigma, 0, 50, TRUE)

# Finally,the obstacle for getting back on the probability scale is that 
# rounding error is always a threat when moving from log-probability to 
# probability. If you use the obvious approach, like exp(post$prod), you’ll
# get a vector full of zeros, which isn’t very helpful.This is a result of 
# R’s rounding very small probabilities to zero.Remember,in large samples, all
# unique samples are unlikely.This is why you have to work with log-probability.
# The code in the box dodges this problem by scaling all of the log-products
# by the maximum log-product.As a result, the values in post$prob are not all 
# zero, but the yals oaren’t exactly probabilities.Instead they are relative 
# posterior probabilities.But that’s good enough for what we wish to do with 
# these values.

post_table$posterior <- exp(post_table$log_posterior - max(post_table$log_posterior))


rethinking::contour_xyz(post_table$mu, post_table$sigma, post_table$posterior)
rethinking::image_xyz(post_table$mu, post_table$sigma, post_table$posterior)


# 2.1 Sampling from the posterior ------------------------------------------

sample_rows <- sample(1:nrow(post_table), size = 1e4, 
                      prob = post_table$posterior, replace = TRUE)
sample_mu <- post_table$mu[sample_rows]
sample_sigma <- post_table$sigma[sample_rows]

plot(sample_mu, sample_sigma,
     cex = 0.5, pch = 16, col = rethinking::col.alpha(rangi2,0.1))

density(sample_mu) %>% plot()
density(sample_sigma) %>% plot()

PI(sample_mu)
PI(sample_sigma)


# 2.2 Posterior predictive Sampling (PPS) ----------------------------------

pps_sample <- rnorm(1e4, sample_mu, sample_sigma)

sum(pps_sample < 0) / 1e4   # Quite credible
sum(pps_sample > 272) / 1e4 # Quite credible


# 3. Stan model -----------------------------------------------------------

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(model_code = stan_code, 
            iter = 1000, chains = 1, 
            data = list(N = length(d2$height), h = d2$height),
            verbose = TRUE)

fit
precis(fit)


# 3.1 Sampling from the posterior -----------------------------------------

sample_mu <- rstan::extract(fit, "mu")[[1]]
sample_sigma <- rstan::extract(fit, "sigma")[[1]]

pairs(fit)

density(sample_mu) %>% plot()
density(sample_sigma) %>% plot()

PI(sample_mu)
PI(sample_sigma)

# 3.2 - Posterior predictive sampling - PPS -------------------------------

stan_code_pps <- "
data{
  int<lower=0> N;
  vector[N] h;
}

parameters {
  real mu;
  real<lower=0,upper=50> sigma;
}

model {
  mu ~ normal(180, 20); // 95% mu in 140 - 220
  sigma ~ uniform(0, 50); // 95 % h in mu +/- 0 - 100
  h ~ normal(mu, sigma); 
}

generated quantities {
  vector[N] h_sim;
  
  for (i in 1:N) {
    h_sim[i] = normal_rng(mu, sigma);
  }
}

"

fit_pps <- stan(model_code = stan_code_pps, 
            iter = 1000, chains = 1, 
            data = list(N = length(d2$height), h = d2$height),
            verbose = TRUE)

fit_pps
precis(fit_pps)

pps_samples <- extract(fit_pps)$h_sim 

# También:
smaple_mu <- sample(extract(fit_pps)$mu, size = 1e4, replace = TRUE)
sample_sigma <- sample(extract(fit_pps)$sigma, size = 1e4, replace = TRUE)

pps_sample2 <- rnorm(1e4, sample_mu, sample_sigma)

sum(pps_samples[1, ] < 0) / 1e4   # Quite credible
sum(pps_samples[1, ] > 272) / 1e4 # Quite credible


plot(density(d2$height), ylim = c(0, 0.06), lwd = 2)

for (i in sample(1:nrow(pps_samples), size = 20)) {
  lines(density(pps_samples[i, ]), col = "red", lty = 2)
}

lines(density(pps_sample2), col = "blue", lty = 1, lwd = 2)

# Sampling from posterior
precis(extract(fit_pps), hist = FALSE)

# Posterior predictive sampling
precis(pps_samples[1, ], hist = FALSE)
precis(pps_samples[10, ], hist = FALSE)
precis(pps_samples[100, ], hist = FALSE)
precis(pps_sample2, hist = FALSE)


# LINEAR PREDICTION -------------------------------------------------------

ggplot(d2, aes(x = weight, y = height)) + geom_point()


stan_code_linear <- "
data {
  int<lower=0> N;
  
  vector [N] weights;
  vector [N] heights;
  real weights_bar;
}

parameters {
  real<lower=0, upper=50> sigma;
  real<lower=0> beta;
  real<lower=0> alpha;
}

transformed parameters {
  real<lower=0> mu[N];
  
  for ( i in 1:N) {
    mu[i] = beta * (weights[i] - weights_bar) + alpha;
  }
}

model {

  sigma ~ uniform(0, 50);
  beta ~ lognormal(0, 1);
  alpha ~ normal(180, 20); // Como la mu del modelo anterior
                           // ya que alfa es el valor que toma mu
                           // cuando el peso es el peso medio
  heights ~ normal(mu, sigma);

}

generated quantities {
  real mu_sim[N];
  vector[N] heights_sim;
 
  for (i in 1:N) {
    // POSTERIOR PREDICTIVE SAMPLING (PPS)
    
    // mu PPS
    mu_sim[i] = beta * (weights[i] - weights_bar) + alpha;
    
    // heights PPS
    heights_sim[i] = normal_rng(mu[i], sigma);
  }
}

"

# 1. Priors exploration ---------------------------------------------------

## Prior Predictive Distribution
set.seed(2971)
N <- 100
a <- rnorm(N, 180, 20)
b <- rnorm(N, 0, 10)

plot(NULL,
     xlim = range(d2$weight),
     ylim = c(-100, 400),
     xlab = "weight",
     ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd =0.5)
mtext("b ~ dnorm(0, 10)")

x_bar <- mean(d2$weight)

for (i in 1:N) {
  curve(a[i] + b[i] * (x - x_bar),
        from = min(d2$weight),
        to = max(d2$weight),
        add = TRUE,
        col = col.alpha("black", 0.2))
}

# Con esta prior la relación entre peso y altura puede ser positiva o negativa
# e incluso a algunos pesos les corresponden alturas negativas o superiores
# a la del hombre más alto del mundo.
# Mejor usar una log-normal:

b <- rlnorm(N, 0, 1)
density(b) %>% plot()

plot(NULL,
     xlim = range(d2$weight),
     ylim = c(-100, 400),
     xlab = "weight",
     ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd =0.5)
mtext("b ~ dnorm(0, 10)")

x_bar <- mean(d2$weight)

for (i in 1:N) {
  curve(a[i] + b[i] * (x - x_bar),
        from = min(d2$weight),
        to = max(d2$weight),
        add = TRUE,
        col = col.alpha("black", 0.2))
}

# Mucho mejor esta prior


# 2. Ulam model -----------------------------------------------------------

ulam_linear <- ulam(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * (weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2 %>% mutate(xbar = x_bar))

ulam_linear
precis(ulam_linear)

# 3. Stan model -----------------------------------------------------------


fit_linear <- stan(model_code = stan_code_linear, 
            iter = 1000, chains = 1, 
            data = list(N = length(d2$height), 
                        heights = d2$height, 
                        weights_bar = mean(d2$weight),
                        weights = d2$weight),
            verbose = TRUE)

print(fit_linear, pars = c("alpha", "beta", "sigma"))
precis(fit_linear)


# UNDERSTANDING THE POSTERIOR ---------------------------------------------

# I emphasize plotting posterior distributions and posterior preictions,
# instead of attempting to understand a table.


# 1. Tablas de DISTRIBUCIONES MARGINALES ----------------------------------

precis(fit_linear)

# Beta = 0.91 => cada kilo más / menos, 1 cm más / menos
# 89% de la probabilidad muy concentrado OK

vcov(ulam_linear)

cov(data.frame(alpha = extract(fit_linear, "alpha"), 
               beta = extract(fit_linear, "beta"), 
               sigma = extract(fit_linear, "sigma")))

# Very little covariation among the parameters in this case.

pairs(ulam_linear)
pairs(fit_linear, pars = c("alpha", "beta", "sigma"))

# The lack of covariance among the parameters results from centering


# 2. Plotting posterior inference vs data ---------------------------------

# It’s almost always much more useful to plot the posteriorinference against
# the data.Not only does plotting help in interpretng the posterior, but it 
# also provides an informal check on model assumptions. Whent he model’s
# predictions don’t come close to key observations or patterns in the plotted
# data, then you might suspect the model either did not fit correctly or is 
# rather badly specified.

## Versión 1 - supersimple
post_samples <- extract(fit_linear, pars = c("alpha", "beta", "sigma"))
a_map <- mean(post_samples$alpha)
b_map <- mean(post_samples$beta)

plot(height ~ weight, d2, col = rangi2)
curve(a_map + b_map * (x - x_bar), add = TRUE)


# a. Visualizar la incertidumbre ------------------------------------------



for (i in sample(length(post_samples$alpha), size = 100)) {
  curve(post_samples$alpha[i] + post_samples$beta[i] * (x - x_bar), 
        col = col.alpha("black", 0.3), 
        add = TRUE)
}


# b. Regression intervals and contours ------------------------------------

mu_at_50 <- post_samples$alpha + post_samples$beta * (50 - x_bar)
dens(mu_at_50,col = rangi2, lwd = 2, xlab = "mu|weight=50")

PI(mu_at_50, prob = 0.89)

# Para todos los valores de weight - LINK function:
mu <- rethinking::link(ulam_linear)
str(mu)

# Igual que:
extract(fit_linear, "mu_sim")[[1]] %>% str()

# O:
sapply(d2$weight, 
       function(x) post_samples$alpha + post_samples$beta * (x - x_bar)) %>% 
  str()

# Si queremos más:
sapply(d2$weight, 
       function(x) sample(post_samples$alpha, size = 1000, replace = TRUE) + 
         sample(post_samples$beta, size = 1000, replace = TRUE) * (x - x_bar)) %>% 
  str()

# Es decir, cada columna es una muestra de la distribución posterior de cada 
# caso (352 individuos)

## Calculemos lo mismo pero no para los weights en d2, sino para los que 
# queramos:
weight_seq <- seq(25, 70, 1)

mu <- link(ulam_linear, data = data.frame(weight = weight_seq, xbar = x_bar))
mu %>% str()

# O:
sapply(weight_seq, 
       function(x) post_samples$alpha + post_samples$beta * (x - x_bar)) %>% 
  str()

# Y ahora pintamos:
plot(height ~ weight, d2, type = "n")

for(i in 1:nrow(mu)) {
  points(weight_seq, mu[i,], pch = 16, col = col.alpha(rangi2, 0.1))
}

# Cada apilamiento es una distribución posterior gaussiana. Hagamos sus 
# summaries:
mu_mean <- apply(mu, 2, mean)
mu_PIs <- apply(mu, 2, PI, 0.89)
mu_HPDIs <- apply(mu, 2, HPDI, 0.89)

# Y pintémoslo:
plot(height ~ weight, d2, col = col.alpha(rangi2, 0.5))
lines(weight_seq, mu_mean)
shade(mu_PIs, weight_seq, col = col.alpha("gray", 1))


# c. Intervalos de predicción ---------------------------------------------

post_sample_alpha <- extract(fit_linear)$alpha
post_sample_beta <- extract(fit_linear)$beta
post_sample_sigma <- extract(fit_linear)$sigma

sim_heights <- sapply(weight_seq,
                      function(x) {
                        
                        post_sample_mu <- post_sample_beta * (x - x_bar) + 
                          post_sample_alpha
                        
                        rnorm(length(post_sample_mu), 
                              post_sample_mu, 
                              post_sample_sigma)
                      })

sim_heights_means <- apply(sim_heights, 2, mean)
sim_heights_PIs <- apply(sim_heights, 2, PI, 0.89)
sim_heights_HPDIs <- apply(sim_heights, 2, HPDI, 0.89)

# Igual que:
sim(ulam_linear, data = list(weight = weight_seq, 
                             xbar = rep(x_bar, length(weight_seq))))

# O:
extract(fit_linear, "heights_sim")[[1]]

# (aunque esto no es con weight_seq)

# Plot de todo

# 1 - Los datos
plot(height ~ weight, d2, col = col.alpha(rangi2, 0.5))

# 2 - La línea de regresión media (MAP line)
lines(weight_seq, mu_mean)

# 3 - Región HPDI de la MAP line
shade(mu_HPDIs, weight_seq, col = col.alpha("blue", alpha = 1))

# 4 - Región PI de las predicciones o simulaciones
shade(sim_heights_PIs, weight_seq, col = col.alpha("grey", alpha = 1))

# SPLINES -----------------------------------------------------------------

data(Howell1)
d <- Howell1

# Ahora height ~ weight completo, sin filtrar no adultos

plot(height ~ weight, d)

## STANDARDIZE the PREDICTOR VARIABLE
# ( the square or cube of a large number can be truly massive)

d$weight_s <- scale(d$weight) %>% as.numeric()

# 1. Polynomial regression ------------------------------------------------

d$weight_s2 <- d$weight_s^2

# a. Stan -----------------------------------------------------------------

stan_code_parabolic <- "
data {
  int<lower=0> N;
  
  vector [N] weights_s;
  vector [N] weights_s2;
  vector [N] heights;
}

parameters {
  real<lower=0, upper=50> sigma;
  real<lower=0> beta1;
  real beta2;
  real<lower=0> alpha;
}

transformed parameters {
  real<lower=0> mu[N];
  
  for ( i in 1:N) {
    mu[i] = alpha + beta1 * weights_s[i] + beta2 * weights_s2[i];
  }
}

model {

  sigma ~ uniform(0, 50);
  beta1 ~ lognormal(0, 1);
  beta2 ~ normal(0, 1);
  alpha ~ normal(180, 20); // Como la mu del modelo anterior
                           // ya que alfa es el valor que toma mu
                           // cuando el peso es el peso medio
  heights ~ normal(mu, sigma);

}

generated quantities {
  real mu_sim[N];
  vector[N] heights_sim;
 
  for (i in 1:N) {
    // POSTERIOR PREDICTIVE SAMPLING (PPS)
    
    // mu PPS
    mu_sim[i] = alpha + beta1 * weights_s[i] + beta2 * weights_s2[i];
    
    // heights PPS
    heights_sim[i] = normal_rng(mu[i], sigma);
  }
}
"

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

fit_parabolic <- stan(model_code = stan_code_parabolic, 
                      iter = 1000, chains = 1, 
                      data = list(N = length(d$height), 
                                  heights = d$height,
                                  weights_s = d$weight_s,
                                  weights_s2 = d$weight_s2),
                      verbose = TRUE)

print(fit_parabolic, pars = c("beta1", "beta2", "sigma"))
precis(fit_parabolic)

# Plot

weight_s_seq <- seq(-2.2, 2, length.out = 30)

post_samples <- extract(fit_parabolic, 
                        pars = c("sigma", "alpha", "beta1", "beta2"))

## PPS de mu - Línea de regresión MAP
mu <- sapply(weight_s_seq, 
       function(x) post_samples$alpha + 
         post_samples$beta1 * x + 
         post_samples$beta2 * x^2) 

mu_mean <- apply(mu, 2, mean)
mu_PIs <- apply(mu, 2, PI, 0.89)
mu_HPDIs <- apply(mu, 2, HPDI, 0.89)


## PPS de la respuesta - Predicción o simulación

sim_heights <- sapply(1:ncol(mu),
                      function(i) {
                        
                        rnorm(length(mu[i,]), 
                              mu[i,], 
                              post_samples$sigma)
                      })

sim_heights_means <- apply(sim_heights, 1, mean)
sim_heights_PIs <- apply(sim_heights, 1, PI, 0.89)
sim_heights_HPDIs <- apply(sim_heights, 1, HPDI, 0.89)

## # Plot de todo

# 1 - Los datos
plot(height ~ weight_s, d, col = col.alpha(rangi2, 0.5))

# 2 - La línea de regresión media (MAP line)
lines(weight_s_seq, mu_mean)

# 3 - Región HPDI de la MAP line
shade(mu_HPDIs, weight_s_seq, col = col.alpha("blue", alpha = 1))

# 4 - Región PI de las predicciones o simulaciones
shade(sim_heights_PIs, weight_s_seq, col = col.alpha("grey", alpha = 1))


# b. ulam -----------------------------------------------------------------

ulam_parabolic <- ulam(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1 * weight_s + b2 * weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 1),
    sigma ~ dunif(0, 50)
    
  ), data = d
)

precis(ulam_parabolic)

mu <- link(ulam_parabolic, 
           data = list(weight_s <- weight_s_seq,
                       weight_s2 <- weight_s_seq^2)) 

mu_mean <- apply(mu, 2, mean)
mu_PIs <- apply(mu, 2, PI, 0.89)
mu_HPDIs <- apply(mu, 2, HPDI, 0.89)

sim_heights <- sim(ulam_parabolic, 
                   data = list(weight_s <- weight_s_seq,
                               weight_s2 <- weight_s_seq^2))

sim_heights_means <- apply(sim_heights, 2, mean)
sim_heights_PIs <- apply(sim_heights, 2, PI, 0.89)
sim_heights_HPDIs <- apply(sim_heights, 2, HPDI, 0.89)

## # Plot de todo

# 1 - Los datos
plot(height ~ weight_s, d, col = col.alpha(rangi2, 0.5))

# 2 - La línea de regresión media (MAP line)
lines(weight_s_seq, mu_mean)

# 3 - Región HPDI de la MAP line
shade(mu_HPDIs, weight_s_seq, col = col.alpha("blue", alpha = 1))

# 4 - Región PI de las predicciones o simulaciones
shade(sim_heights_PIs, weight_s_seq, col = col.alpha("grey", alpha = 1))

# Con ggplot:

gg_data <- data.frame(
  weight_s_seq = weight_s_seq,
  mu_mean = mu_mean,
  min_mu_HPDIs = mu_HPDIs[1,],
  max_mu_HPDIs = mu_HPDIs[2,],
  min_sim_PIs = sim_heights_PIs[1,],
  max_sim_PIs = sim_heights_PIs[2,],
  min_sim_HPDIs = sim_heights_HPDIs[1,],
  max_sim_HPDIs = sim_heights_HPDIs[2,]
)

# Datos
p <- ggplot() + 
  geom_point(data = d, mapping = aes(x = weight_s, 
                                     y = height, 
                                     color = col.alpha(rangi2, 0.5)))

# Línea de regresión media (MAP line)
p <- p + 
  geom_line(data = gg_data,
            aes(x = weight_s_seq, y = mu_mean)) 

# Región HPDI de la MAP line
p <- p + 
  geom_line(data = gg_data,
            aes(x = weight_s_seq, y = min_mu_HPDIs), linetype = 2) + 
  geom_line(data = gg_data,
            aes(x = weight_s_seq, y = max_mu_HPDIs), linetype = 2) +
  geom_ribbon(data = gg_data,
              aes(x = weight_s_seq, ymin = min_mu_HPDIs, ymax = max_mu_HPDIs),
              fill = "blue", alpha = 0.2)

# Región PI de las predicciones o simulaciones
p <- p  + 
  geom_line(data = gg_data,
            aes(x = weight_s_seq, y = min_sim_PIs), linetype = 2) + 
  geom_line(data = gg_data,
            aes(x = weight_s_seq, y = max_sim_PIs), linetype = 2) +
  geom_ribbon(data = gg_data,
              aes(x = weight_s_seq, ymin = min_sim_PIs, ymax = max_sim_PIs),
              fill = "gray", alpha = 0.2)

p

# 2. Splines --------------------------------------------------------------

library(rethinking)

data(cherry_blossoms)

d <- cherry_blossoms

precis(d, hist = FALSE)

plot(doy ~ year, d)

d2 <- d[complete.cases(d %>% select(doy, year)), ]

num_knots <- 15
knot_list <- quantile(d2$year, probs = seq(0, 1, length.out = num_knots))

library(splines)

# Cubic spline
B <- bs(d2$year,
        knots = knot_list[-c(1, num_knots)],
        degree = 3, intercept = TRUE)

plot(NULL, 
     xlim = range(d2$year), 
     ylim=c(0,1),
     xlab="year",
     ylab="basis")

for (i in 1:ncol(B)) lines(d2$year, B[,i])


# a. Stan -----------------------------------------------------------------

stan_code_spline <- "
data {
  int<lower = 0> N;
  int<lower=0> N_yrs;
  int<lower=0> N_basis;
  vector[N] doy;
  matrix[N_yrs, N_basis] B;
}

parameters {
  real<lower = 0> a;
  vector[N_basis] w;
  real <lower = 0> sigma;
}

transformed parameters {
  row_vector[N_yrs] mu;
  
  for (i in 1:N_yrs) {
    mu[i] = a + B[i, ] * w;
  }
}

model {
  a ~ normal(100, 10);
  w ~ normal(0, 10);
  sigma ~ exponential(1);
    
  doy ~ normal(mu, sigma);
}

generated quantities {
  real mu_sim[N];
  vector[N] doy_sim;

  for (i in 1:N) {
    // POSTERIOR PREDICTIVE SAMPLING (PPS)
    
    // mu PPS
    // mu_sim[i] = alpha + beta1 * weights_s[i] + beta2 * weights_s2[i];
    mu_sim[i] = a + B[i, ] * w;
    
    // heights PPS
    doy_sim[i] = normal_rng(mu[i], sigma);
  }
}

"

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


fit_spline <- stan(model_code = stan_code_spline, 
                   iter = 1000, chains = 1, 
                   data = list(N = length(d2$doy), 
                               doy = d2$doy, 
                               B = B,
                               N_yrs = nrow(B),
                               N_basis = ncol(B)),
                   verbose = TRUE)

print(fit_spline, pars = c("a", "w", "sigma"))
precis(fit_spline)

##
post <- extract.samples(fit_spline)
w <- apply(post$w, 2, mean)

plot(NULL, xlim=range(d2$year), ylim = c(-6, 6),
      xlab = "year", ylab = "basis*weight")

for (i in 1:ncol(B)) lines(d2$year, w[i]*B[,i])

##

post_samples <- extract(fit_spline, pars = c("sigma", "a", "w"))

## PPS de mu - Línea de regresión MAP
mu <- sweep(post_samples$w %*% t(B), 
            MARGIN = 1, post_samples$a, "+")

mu_mean <- apply(mu, 2, mean)
mu_PIs <- apply(mu, 2, PI, 0.89)
mu_HPDIs <- apply(mu, 2, HPDI, 0.89)

## PPS de la respuesta - Predicción o simulación

sim_doys <- sapply(1:nrow(mu),
                   function(i) {
                     
                     rnorm(length(mu[i,]), 
                           mu[i, ], 
                           post_samples$sigma)
                   })

sim_doys_means <- apply(sim_doys, 1, mean)
sim_doys_PIs <- apply(sim_doys, 1, PI, 0.89)
sim_doys_HPDIs <- apply(sim_doys, 1, HPDI, 0.89)

plot(d2$year, d2$doy, col = col.alpha(rangi2,0.3), pch=16)
shade(mu_PIs, d2$year, col = col.alpha("black",0.5))
shade(sim_doys_PIs, d2$year, col = col.alpha("black",0.5))

# a. ulam -----------------------------------------------------------------


ulam_spline <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    # mu <-a + sapply(1:827, function(i) sum(B[i,] * w)),
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data = list(D = d2$doy, B = B),
  start = list(w = rep(0, ncol(B)))
)

# EJERCICIOS --------------------------------------------------------------


# Medium ------------------------------------------------------------------


# 4M1 ---------------------------------------------------------------------

# For the model definition below, simulate observed y values from the prior 
# (not the posterior).
#                      y_i ~ Normal(mu, sigma)
#                      mu ~ Normal(0, 10)
#                      sigma ~ Exponential(1)

sigma <- rexp(1e4, 1)
mu <- rnorm(1e4, 0, 10)
y <- rnorm(1e4, mu, sigma)

dens(y)

# 4M2 ---------------------------------------------------------------------

# Trasnlate the model just above into a quap formula.

flist_4M2 <- alist(
  y ~ dnorm(mu, sigma),
  mu ~ dnorm(0, 10),
  sigma ~ dexp(1)
)

# 4M3 ---------------------------------------------------------------------

# Translate the quap model formula below into a mathematical model definition.

alist(y ~ dnorm(mu, sigma),
      mu <- a + b * x,
      a ~ dnorm(0, 10),
      b ~ dunif(0, 1),
      sigma ~ dexp(1))

# Trivial

# 4M4 ---------------------------------------------------------------------

# A sample of students is measured for height each year for 3 years. After the 
# third year, you want to fit a linear regression predicting height using year 
# as a predictor. Write down the mathematical model definition for this 
# regression, using any variable names and priors you choose. Be prepared to 
# defend your choice of priors.

# height ~ Normal(mu, sigma)
# mu = alpha + beta * (year - bar_year)
# alpha ~ Normal(170, 20)
# beta ~ lLognormal(0, 1)
# sigma ~ Exponetial(1)

alpha_prior <- rnorm(1e4, 170, 20)
beta_prior <- rlnorm(1e4, 0, 1)

plot(NULL, xlim = c(0, 3), ylim = c(0, 220))

for (i in sample(1:1e5, size = 1000)) {
  curve(alpha_prior[i] + beta_prior[i] * (x - 2), add = TRUE)
}

# 4M5 ---------------------------------------------------------------------

# Now suppose I remind you that every student got taller each year. Does this 
# information lead you to change your choice of priors? How?

# 4M6 ---------------------------------------------------------------------

# Now suppose I tell you that the variance among heights for students of the 
# same age is never more than 64cm. How does this lead you to revise your 
# priors?

# 4M7 ---------------------------------------------------------------------

# Refit model m4.3 from the chapter, but omit the mean weight xbar this time. 
# Compare the new model’s posterior to that of the original model. In 
# particular, look at the covariance among the parameters. What is different? 
# Then compare the posterior predictions of both models.

# 4M8 ---------------------------------------------------------------------

# In the chapter, we used 15 knots with the cherry blossom spline. Increase the 
# number of know and observe what happens to the resulting spline. Then adjust 
# also the width of the prior on the weights—change the standard deviation of 
# the prior and watch what happens. What do you think the combination of know 
# number and the prior on the weights controls?


# Hard --------------------------------------------------------------------


# 4H1 ---------------------------------------------------------------------

# The weights listed below were recorded in the !Kung census, but heights were 
# not recorded for these individuals. Provide predicted heights and 89% 
# intervals for each of these individuals. That is, fill in the table, below, 
# using model-based predictions.

#   Individual	weight	expected height	89% interval
#       1	       47.0		
#       2	       43.7		
#       3	       64.8		
#       4	       32.6		
#       5	       54.6		

data(Howell1)

d <- Howell1
d2 <- Howell1 %>% filter(age >= 18)

d_pred <- tibble(weight = c(47.0, 43.7, 64.8, 32.6, 54.6))

fit_linear <- stan(model_code = stan_code_linear, 
                   iter = 1000, chains = 1, 
                   data = list(N = length(d2$height), 
                               heights = d2$height, 
                               weights_bar = mean(d2$weight),
                               weights = d2$weight),
                   verbose = TRUE)

post_betas <- extract(fit_linear, "beta")[[1]]
post_alphas <- extract(fit_linear, "alpha")[[1]]
post_sigmas <- extract(fit_linear, "sigma")[[1]]
weights_bar <- d2$weight %>% mean()

post_mus <- sapply(1:length(post_betas),
                   function(i) {
                     post_alphas[i] + post_betas[i] * (d_pred$weight - weights_bar)
                   })
post_preds <- sapply(1:nrow(post_mus),
                     function(i) {
                       rnorm(length(post_mus[i,]), post_mus[i,], post_sigmas[i])
                     })

apply(post_preds, 2, mean)

apply(post_preds, 2, HPDI) %>% t()

## Predicción directamente desde Stan:
# https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed

stan_code_pred_linear <- "
data {
  int<lower=0> N_pred;
  
  real weights[N_pred];
  real weights_bar;
  
  int<lower = 0> N_samples;
  real<lower=0, upper=50> sigma[N_samples];
  real<lower=0> beta[N_samples];
  real<lower=0> alpha[N_samples];
}

parameters {

}

transformed parameters {

}

model {

}

generated quantities {
  # matrix[N_samples, N_pred] mu_sim;
  # real heights_sim[N_samples, N_pred];
  # 
  # // POSTERIOR PREDICTIVE SAMPLING (PPS)
  # for (n in 1:N_pred) {
  #   for (i in 1: N_samples) {
  #     // mu PPS
  #     mu_sim[i, n] = beta[i] * (weights[n] - weights_bar) + alpha[i];
  #     // heights PPS
  #     heights_sim[i, n] = normal_rng(mu_sim[i, n], sigma)[1];
  #   }
  # }

  matrix[N_pred, N_samples] mu_sim;
  real heights_sim[N_pred, N_samples];

  // POSTERIOR PREDICTIVE SAMPLING (PPS)
  for (n in 1:N_pred) {
    for (i in 1: N_samples) {
      // mu PPS
      mu_sim[n, i] = beta[i] * (weights[n] - weights_bar) + alpha[i];

    }
    // heights PPS
      heights_sim[n] = normal_rng(mu_sim[n], sigma);
  }

}
"

pred <- stan(model_code = stan_code_pred_linear, 
             iter = 1, chains = 1, 
             data = list(N_pred = length(d2$height), 
                         weights_bar = mean(d2$weight),
                         weights = d2$weight,
                         N_samples = length(post_alphas),
                         alpha = post_alphas,
                         beta = post_betas,
                         sigma = post_sigmas),
             verbose = TRUE,
             algorithm = "Fixed_param")

ext_pred <- extract(pred)


# 4H2 ---------------------------------------------------------------------

# Select out all the rows in the Howell1 data with ages below 18 years of age. 
# If you do it right, you should end up with a new data frame with 192 rows in 
# it.


# a. Fit a linear regression to these data, using quap. Present and interpret  
# the estimates. For every 10 units of increase in weight, how much taller does  
# the model predict a child gets?

# b. Plot the raw data, with height on the vertical axis and weight on the 
# horizontal axis. Superimpose the MAP regression line and 89% interval for the
# mean. Also superimpose the 89% interval for predicted heights.

# c. What aspects of the model fit concern you? Describe the kinds of 
# assumptions you would change, if any, to improve the model. You don’t have 
# to write any new code. Just explain what the model appears to be doing a bad 
# job of, and what you hypothesize would be a better model.

# 4H3 ---------------------------------------------------------------------


# 4H4 ---------------------------------------------------------------------


# 4H5 ---------------------------------------------------------------------


# 4H6 ---------------------------------------------------------------------


# 4H7 ---------------------------------------------------------------------


# 4H8 ---------------------------------------------------------------------


