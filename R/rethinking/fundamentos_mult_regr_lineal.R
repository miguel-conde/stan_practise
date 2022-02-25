library(tidyverse)

library(rethinking)

library(tidybayes)

library(dagitty)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# A - ASOCIACIONES ESPÚREAS -----------------------------------------------

data("WaffleDivorce")

d <- WaffleDivorce %>% as_tibble()
glimpse(d)

ggplot(data = d, mapping = aes(x = Marriage, y = Divorce)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = d, mapping = aes(x = MedianAgeMarriage, y = Divorce)) +
  geom_point() +
  geom_smooth(method = "lm")

d %>% select(Marriage, Divorce, MedianAgeMarriage) %>% cor()

# Marriage (rate) y MedianAgeMarriage son buenos PREDICTORES de Divorce (rate),
# hay una clara CORRELACiÓN, están ASOCIADOS  Pero ¿la relación es CAUSAL?

d <- d %>% mutate(D = standardize(Divorce),
                  M = standardize(Marriage),
                  A = standardize(MedianAgeMarriage))

sc_linear_D_A <- " // Stan code dor lineal model D ~ A
data {
  int<lower=1> n;
  vector[n] A;
  vector[n] D;
}

parameters {
  real<lower=0> sigma;
  real<lower=0> alpha;
  real beta;
  
}

transformed parameters {
  vector[n] mu;
  mu = alpha + beta * A;
}

model {

  sigma ~ exponential(1);
  alpha ~ normal(0, 0.2);
  beta ~ normal(0, 0.5); // solo con 5% de las betas mu más allá de 1 desv. 
                         // estándar de la tasa de divorcio media
  
  D ~ normal(mu, sigma);
}

generated quantities {
  vector[n] mu_sim;
  real D_sim[n];
  
  // PPS
  for (i in 1:n) {
    mu_sim[i] = alpha + beta * A[i];
    D_sim[i] = normal_rng(mu_sim[i], sigma);
  }
}
"

# PRIOR PREDICTIVE SAMPLING
prior_sigma <- rexp(1e5, 1)
prior_beta <- rnorm(1e5, 0, 0.5)
prior_alpha <- abs(rnorm(1e5, 0, 0.2))

seq_A <- seq(-2, 2, length.out = 100)

prior_mu <- sapply(seq_A, function(x) prior_alpha + prior_beta * x) %>% t()
prior_D <- prior_mu %>% apply(2, function(x) rnorm(x, prior_sigma))

prior_sample_cols <- sample(ncol(prior_mu), 100)

p <- ggplot(data = tibble(A = seq_A), aes(x = A))

for (i in 1:100) {
  print(head(prior_mu[,i]))
    p <- p + geom_line(aes(y = prior_mu[,i]))
}

p

## Fit the model

model_linear_D_A <- stan_model(model_code = sc_linear_D_A)

fit_linear_D_A <- sampling(model_linear_D_A, 
                           iter = 1000, chains = 1, 
                           data = compose_data(d %>% select(D, A)),
                           verbose = TRUE)

print(fit_linear_D_A, pars = c("alpha", "beta", "sigma"))
precis(fit_linear_D_A)

## POSTERIRO PREDICTIVE SAMPLING

post <- extract(fit_linear_D_A, c("alpha", "beta", "sigma"))
post_mu <-  sapply(seq_A, function(x) post$alpha + post$beta * x) %>% t()
prior_D <- post_mu %>% apply(2, function(x) rnorm(x, prior_sigma))

post_mu_mean <- apply(post_mu, 1, mean)
post_mu_PI <- apply(post_mu, 1, PI)

plot(D ~ A, data = d, col = rangi2)
lines(seq_A, post_mu_mean, lwd = 2)
shade(post_mu_PI, seq_A)


# 1. DAGs e independencias condicionales -----------------------------------


# THINK BEFORE REGRESS
# TESTABLE IMPLICATIONS

library(dagitty)

DAG_1 <- dagitty('
                 dag{D <- A -> M -> D}
                 ')

plot(DAG_1)

DAG_2 <- dagitty('
                 dag{D <- A -> M}
                 ')

plot(DAG_2)

impliedConditionalIndependencies(DAG_1)
impliedConditionalIndependencies(DAG_2)


# 2 - Regresión multivariable ---------------------------------------------

sc_linear_D_A_M <- " // Stan code dor lineal model D ~ A
data {
  int<lower=1> n;
  vector[n] A;
  vector[n] D;
  vector[n] M;
}

parameters {
  real<lower=0> sigma;
  real<lower=0> alpha;
  real beta_A;
  real beta_M;
  
}

transformed parameters {
  vector[n] mu;
  mu = alpha + beta_A * A + beta_M * M;
}

model {

  sigma ~ exponential(1);
  alpha ~ normal(0, 0.2);
  beta_A ~ normal(0, 0.5); 
  beta_M ~ normal(0, 0.5); 
  
  D ~ normal(mu, sigma);
}

generated quantities {
  vector[n] mu_sim;
  real D_sim;
  
  // PPS
  for (i in 1:n) {
    mu_sim[i] = alpha + beta_A * A[i]+ beta_M * A[i];
    D_sim = normal_rng(mu_sim[i], sigma);
  }
}
"

## Fit the model

model_linear_D_A_M <- stan_model(model_code = sc_linear_D_A_M)

fit_linear_D_A_M <- sampling(model_linear_D_A_M, 
                           iter = 1000, chains = 1, 
                           data = compose_data(d %>% select(D, A, M)),
                           verbose = TRUE)

print(fit_linear_D_A_M, pars = c("alpha", "beta_A", "beta_M", "sigma"))
precis(fit_linear_D_A_M)

plot(fit_linear_D_A_M, pars = c("alpha", "beta_A", "beta_M", "sigma"))

# 3 - Plotting multivariate posteriors ------------------------------------


# a. Predictor residual plots ---------------------------------------------


# b. Posterior prediction plots -------------------------------------------


# c. Counterfactual plots -------------------------------------------------


# B - RELACIONES OCULTAS --------------------------------------------------

data("milk")

d <- milk

glimpse(milk)
summary(milk)

# C - VARIABLES CATEGÓRICAS -----------------------------------------------


# EJERCICIOS --------------------------------------------------------------


# Medium ------------------------------------------------------------------


# 5M1 ---------------------------------------------------------------------

# Invent your own example of a spurious correlation. An outcome variable 
# should be correlated with both predictor variables. But when both predictors 
# are entered in the same model, the correlation between the outcome and one 
# of the predictors should mostly vanish (or at least be greatly reduced).

# SIMULACIÓN del MODELO REAL
x1 <- rnorm(1000)
# x1 -> x2
x2 <- rnorm(1000, x1) # x1 es CAUSA de x2 => estarán correladas
# x1 -> y
y  <- rnorm(1000, x1) # x1 también es CAUSA de la respuesta y => x1 e y estarán
                      # correladas; pero también estarán correladas x2 e y 
                      # debido a su origen común, aunque x2 e y no son CAUSA 
                      # una de Otra.

# Comprobémoslo:
cor(tibble(x1, x2, y))

# Si encontramos un dataset así - como en el caso de la edad (A = x1), la tasa 
# de matrimonio (M = x2) y la de Divorcio (D = y) - la primera hipótesis sería:

DAG_H1 <- dagitty(('
                 dag{x1 -> x2 -> y
                     x1 -> y}
                 '))

plot(DAG_H1)

# En este dag hay implicaciones que podemos comprobar empíricamente con los 
# datos:

# 1. x1 y x2 no son independientes (porque hay una flecha entre ellas) => 
#    => están correladas. Si lo comprobamos y resulta que no están correladas,
#    nos echaría por tierra al menos esta flecha en el modelo.
# 2. x1 e y no son independientes - idem
# 3. x2 e y no son independientes - idem

# Estas 3 ya las hemos testado arriba y el modelo hipotético ha sobrevivido.
# También podríamos haber hecho:
summary(lm(x2 ~ x1))
summary(lm(y ~ x1))
summary(lm(y ~ x2))

# En los 3 casos los coeficientes son significativos => hay correlación

# Veamos ahora independencias condicionadas
dagitty::impliedConditionalIndependencies(DAG_H1)

# No tiene - lógico porque hay flecha entre todas los posibles pares de 
# variables. 
# Comprobemos con los datos:

# ¿Siguen siendo independientes x1 y x2 si condicionamos en y?
summary(lm(x2 ~ x1 + y))

# Si, la correlación entre x2 y x1 persiste

# Y ¿x1 con y si condicionamos en x2?
summary(lm(y ~ x1 + x2))

# Si, la correlación entre x1 e y persiste.

# Por último, ¿siguen siendo independientes y y x2 si condicionamos en x1?
# En el mismo modelo anterior se ve que no (porque el coeficiente de x2 no
# es siginificativo). Por tanto, y en contra de nuestra primera hipótesis (el
# modelo DAG_H1) si hay una independencia condicional: x2 e y son independientes
# si condicionamos en x1. Y si la hipótesis es falsa, el modelo no es correcto.

# Lo modificamos:
DAG_H2 <- dagitty(('
                 dag{x1 -> x2
                     x1 -> y}
                 '))

plot(DAG_H2)

dagitty::impliedConditionalIndependencies(DAG_H2)

# Como se ve, esta segunda hipótesis lleva implícita la independencia condicional
# de que x2 e y condicionada a x1.

# 5M2 ---------------------------------------------------------------------

# Invent your own example of a masked relationship. An outcome variable should 
# be correlated with both predictor variables, but in opposite directions. And 
# the two predictor variables should be correlated with one another.

# SIMULACIÓN del MODELO REAL
x1 <- rnorm(100, 1)
# x1 -> x2
x2 <- rnorm(100, x1)
# x1 -> y <- x2
y  <- rnorm(100, x2 - x1)

d1 <- tibble(x1, x2, y)

cor(d1)

# En este caso, aunque x1 es causa de y, su correlación aparece débil, 
# enmascarada su relación por x2, que es causada también por x1 y es a su vez
# la segunda causa de y.
# De esas correlaciones podríamos pensar que el modelo correcto es:

DAG_H1 <- dagitty('
                  dag{x1 -> y <- x2
                      x1 -> x2}
                  ')
plot(DAG_H1)

# o cualquiera de sus equivalentes de Markov

equivalenceClass(DAG_H1)

# Por ejemplo:

x2 <- rnorm(100, 1)
# x2 -> x1
x1 <- rnorm(100, x2)
# x1 -> y <- x2
y  <- rnorm(100, x2 - x1)

d2 <- tibble(x1, x2, y)

cor(d2)

DAG_H2 <- dagitty('
                  dag{x1 -> y <- x2
                      x1 <- x2}
                  ')
plot(DAG_H2)

# O:

U <- rnorm(100, 1)
# U -> x2
x2 <- rnorm(100, U)
# U -> x1
x1 <- rnorm(100, U)
# x1 -> y <- x2
y  <- rnorm(100, x2 - x1)

d3 <- tibble(U, x1, x2, y)

cor(d3)

DAG_H3 <- dagitty('
                  dag{x1 -> y <- x2
                      x1 <- U -> x2}
                  ')
plot(DAG_H3)

equivalenceClass(DAG_H3)

## Modelos con solo 1 regresor:

lm1 <- lm(y ~ x1)
lm2 <- lm(y ~ x2)

summary(lm1)
summary(lm2)

old_par <- par(mfrow = c(2, 2))

plot(d3 %>% select(x1, y), main = "y ~ x1")
abline(lm1)

plot(d3 %>% select(x2, y), main = "y ~ x2")
abline(lm2)

# La asociación entre el regresor y la respuesta es más débil de lo que debería
# por el efecto del otro regresor (correlado negativamente y no incluido en el
# modelo) a través de U. Para romper esa influencia podríamos condicionar en U
# si lo conociéramos:
summary(lm(y ~ U + x1, d3))
summary(lm(y ~ U + x2, d3))

## Lo que hay que hacer es incluir los 2 regresores, ya que al condicionar en
# x1 / x2 (son chains) se rompe la influencia de x2 / x1 a través de U:

lm12 <- lm(y ~ x1 + x2)
summary(lm12)

# En este la asociación de y con los 2 regresores es mayor que en los modelos
# con solo 1 regresor. Y prácticamente la misma que en los modelos y ~ x1 + U e
# y ~ x2 + U

plot(d3 %>% select(x1, y))
title(main = "y ~ x1 + x2", 
      sub = "COUNTERFACTUAL x2 = 0")
abline(a = coef(lm12)[1], b = coef(lm12)[2])

plot(d3 %>% select(x2, y))
title(main = "y ~ x1 + x2", 
      sub = "COUNTERFACTUAL x1 = 0")
abline(a = coef(lm12)[1], b = coef(lm12)[3])

par(old_par)

# 5M3 ---------------------------------------------------------------------

# It is sometimes observed that the best predictor of fire risk is the presence 
# of firefighters—State and localities with many firefighters also have more 
# fires. Presumably firefighters do not cause fires. Nevertheless, this is not a 
# spurious correlation. Instead fires cause firefighters. Consider the same 
# reversal of causal inference in the context of the divorce and marriage data. 
# How might a high divorce rate cause a higher marriage rate? Can you think of a 
# way to evaluate this relationship, using multiple regression?

# Hard --------------------------------------------------------------------


# 5H1 ---------------------------------------------------------------------

# In the divorce example, suppose the DAG is:  M → A → D
# What are the implied conditional independencies of the graph? Are the data 
# consistent with it?

DAG_5H1 <- dagitty('
                 dag{M -> A -> D}
                 ')

plot(DAG_5H1)

impliedConditionalIndependencies(DAG_5H1)

# Al condicionar en A cerramos el pipe y M deja de influir en D


# 5H2 ---------------------------------------------------------------------

# Assuming that the DAG for the divorce example is indeed  M → A → D, fit a 
# new model and use it to estimate the counterfactual effect of halving a 
# State’s marriage rate M. Using the counterfactual example from the chapter 
# (starting on page 140) as a template.

data("WaffleDivorce")

d <- WaffleDivorce %>% as_tibble()

d <- d %>% mutate(D = standardize(Divorce),
                  M = standardize(Marriage),
                  A = standardize(MedianAgeMarriage))

# Si M -> A -> D es el DAG correcto, hemos visto en el ejercicio anterior que,
# para obtener el efecto causal de M sobre D, no debemos condicionar en A.

sc_linear_5H2 <- " // Stan code for lineal model D ~ A + M, A ~ M
data {
  int<lower=1> n;
  vector[n] M;
  vector[n] A;
  vector[n] D;
}

parameters {
  real<lower=0> sigma_DAM;
  real<lower=0> alpha_DAM;
  real beta_DAM_A;
  real beta_DAM_M;
  
  
  real<lower=0> sigma_AM;
  real<lower=0> alpha_AM;
  real beta_AM;
  
}

transformed parameters {
  vector[n] mu_DAM;
  vector[n] mu_AM;
  
  mu_DAM = alpha_DAM + beta_DAM_A * A + beta_DAM_M * M;
  
  mu_AM = alpha_AM + beta_AM * M;
}

model {

  sigma_DAM ~ exponential(1);
  alpha_DAM ~ normal(0, 0.2);
  beta_DAM_A ~ normal(0, 0.5); 
  beta_DAM_M ~ normal(0, 0.5); 
  
  D ~ normal(mu_DAM, sigma_DAM);
  
  sigma_AM ~ exponential(1);
  alpha_AM ~ normal(0, 0.2);
  beta_AM ~ normal(0, 0.5); 
  
  A ~ normal(mu_AM, sigma_AM);
}

generated quantities {
  vector[n] mu_DAM_sim;
  real D_DAM_sim[n];
  vector[n] log_lik_DAM;
  
  
  vector[n] mu_AM_sim;
  real A_AM_sim[n];
  vector[n] log_lik_AM;
  
  // PPS
  for (i in 1:n) {
    mu_DAM_sim[i] = alpha_DAM + beta_DAM_A * A[i] + beta_DAM_M * M[i];
    D_DAM_sim[i] = normal_rng(mu_DAM_sim[i], sigma_DAM);
    
    mu_AM_sim[i] = alpha_AM + beta_AM * M[i];
    A_AM_sim[i] = normal_rng(mu_AM_sim[i], sigma_AM);
    
    log_lik_DAM[i] = normal_lpdf(D[i] | mu_DAM[i], sigma_DAM);
    
    log_lik_AM[i] = normal_lpdf(A[i] | mu_AM[i], sigma_AM);
  }
}
"

model_linear_5H2 <- stan_model(model_code = sc_linear_5H2)

fit_linear_5H2 <- sampling(model_linear_5H2, 
                             iter = 1000, chains = 1, 
                             data = compose_data(d %>% select(D, A, M)),
                             verbose = TRUE)

print(fit_linear_5H2, pars = c("alpha_DAM", "beta_DAM_M", "beta_DAM_A", "sigma_DAM",
                               "alpha_AM", "beta_AM", "sigma_AM"))
precis(fit_linear_5H2)

plot(fit_linear_5H2, pars = c("alpha", "beta_M", "sigma"))

post <- extract(fit_linear_5H2)

M_seq <- seq(-2, 3, length = 30)

sim_dat <- tibble(M = M_seq)

sim_dat$A_sim <- sapply(M_seq, 
                        function(m) post$alpha_AM + post$beta_AM * m) %>% 
  apply(2, mean)

sim_dat$D_sim <- sapply(seq_along(M_seq), 
                        function(i) {
                          post$alpha_DAM + post$beta_DAM_M * M_seq[i] +
                            post$beta_DAM_A * sim_dat$A_sim[i]
                        }) %>% 
  apply(2, mean)

plot(D ~ M, data = d)
lines(D_sim ~ M, sim_dat)

lm_AM <- lm(A ~ M, d)
lm_DAM <- lm(D ~ A + M, d)

sim_D <- coef(lm_DAM)[1] + 
  (coef(lm_AM)[1] + coef(lm_AM)[2] * M_seq) * coef(lm_DAM)[2] + 
   coef(lm_DAM)[3] * M_seq
lines(M_seq, sim_D, col = "red")

effect_M_on_A <- coef(lm_AM)[1] + coef(lm_AM)[2] * M_seq
total_effect_M_on_D <- (coef(lm_AM)[2] * coef(lm_DAM)[2] + coef(lm_DAM)[3]) * M_seq

old_par <- par(mfrow = c(1, 2))
plot(M_seq, effect_M_on_A, type = "l", main = "Effect M on A")
plot(M_seq, total_effect_M_on_D, type = "l", main = "Total effect M on D")
par(old_par)

# 5H3 ---------------------------------------------------------------------

# Return to the milk energy model, m5.7. Suppose that the true causal 
# relationship among the variables is: M -> K <- N, M -> N
# Now compute the counterfactual effect on K of doubling M. You will need to 
# account for both the direct and indirect paths of causation. Use the 
# counterfactual example from the chapter (starting on page 140) as a template.

data("milk")

d <- milk %>% 
  select(M = mass, N = neocortex.perc, K = kcal.per.g) %>% 
  mutate_at(vars(M, N, K), ~ standardize(.)) %>% 
  as_tibble() %>% 
  drop_na()

dag_5H3 <- dagitty('dag{M -> K <- N
                        M -> N}')

plot(dag_5H3)

impliedConditionalIndependencies(dag_5H3)
adjustedNodes(dag_5H3)

lm_N_M <- lm(N ~ M, d)
summary(lm_N_M)

lm_K_MN <- lm(K ~ M + N, d)
summary(lm_K_MN)

seq_M <- seq(-2, 2, length.out = 30)

sims <- tibble(M = seq_M)
sims$N = predict(lm_N_M, sims)
sims$K <- predict(lm_K_MN, sims)

old_par <- par(mfrow = c(1, 3))
plot(N ~ M, data = sims, type = "l", main = "Effect of M on N")
plot(K ~ M, data = sims, type = "l", main = "Total Effect of M on K")
plot(K ~ N, data = sims, type = "l", main = "Total Effect of N on K")
par(old_par)

sims_x_2 <- sims %>% mutate(M = 2 * M)
sims_x_2$N = predict(lm_N_M, sims_x_2)
sims_x_2$K <- predict(lm_K_MN, sims_x_2)

old_par <- par(mfrow = c(1, 3))
plot(N ~ M, data = sims, type = "l", main = "Effect of M on N")
lines(N ~ M, data = sims_x_2, col = "red", main = "Effect of M on N")
plot(K ~ M, data = sims, type = "l", main = "Total Effect of M on K")
lines(K ~ M, data = sims_x_2, col = "red", main = "Total Effect of M on K")
plot(K ~ N, data = sims, type = "l", main = "Total Effect of N on K")
lines(K ~ N, data = sims_x_2, col = "red", main = "Total Effect of N on K")
par(old_par)