library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# Heights Simple Example --------------------------------------------------

# Fake data
Y <- rnorm(10, 1.5, 0.2)

fit <- stan(here::here("STAN", "height.stan"), 
            iter = 200, chains = 4, data = list(N = length(Y), y = Y))

print(fit, probs = c(0.25, 0.5, 0.75))

library(ggplot2) 

mu <- extract(fit, "mu")[[1]] 

qplot(mu)

library(shinystan) 

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)


# Heights and weights Example ---------------------------------------------

N <- 100

# Fake data
X <- rnorm(N, 60, 10)
beta <- 0.3
sigma <- 0.3

Y <- beta * log(X) + rnorm(N, 0, sigma)

fit <- stan(here::here("STAN", "heights_weights.stan"), 
            iter = 200, chains = 4, data = list(N = length(Y), y = Y, x = X))

print(fit, probs = c(0.25, 0.5, 0.75))

mu <- extract(fit, "mu")[[1]] 

qplot(mu)

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)


# Heights with transformed parameters and transformed data ----------------

# Fake data
Y <- rnorm(10, 1.6, sqrt(5.1))

fit <- stan(here::here("STAN", "height_2.stan"), 
            iter = 1000, chains = 4, data = list(N = length(Y), y = Y))

print(fit, probs = c(0.25, 0.5, 0.75))

library(ggplot2) 

mu <- extract(fit, "mu")[[1]] 

qplot(mu)

library(shinystan) 

aFit <- as.shinystan(fit) 

launch_shinystan(aFit)


# Independent samples -----------------------------------------------------

fit <- stan(here::here("STAN", "ind_samples.stan"), 
            data = list(mu = 10, kappa = 5), 
            algorithm = "Fixed_param", 
            iter = 4000, 
            chains = 1)

Y <- extract(fit, "y")[[1]] 

qplot(Y) + geom_histogram(binwidth = 2)


# Translate and compile a model for later use -----------------------------

aModel <- stan_model(here::here("STAN", "ind_samples.stan")) 

fit <- sampling(aModel, 
                data = list(mu = 10, kappa = 5), 
                algorithm = "Fixed_param", 
                iter = 4000, 
                chains = 1)

# GROUPS ------------------------------------------------------------------

# Suppose that you have individual data for three studies of individuals’ 
# heights. In particular, imagine that you have the following data in R, and you 
# want to generate a separate estimate of the mean population height for each of 
# the three cases: 
X_1 <- c( 1.53, 1.67, 1.52) 
X_2 <- c( 1.75, 1.62, 1.87, 1.95) 
X_3 <- c( 1.25, 1.75)

# What is the best way to pass this data to Stan?

### SOLUTION 1

# A nice trick is to combine all data into one long data vector. We then create 
# helper arrays in R that indicate the size of each data set and its starting 
# position in the long array: 
Y     <- c(1.53, 1.67, 1.52, 1.75, 1.62, 1.87, 1.95, 1.25, 1.75) 
S     <- c(3, 4 ,2)     # sample sizes of each study 
index <- c(1, 4, 8)     # start position of each

# Fake data
N <- 100
Y <- rnorm(N, 1.5, 0.2)
S <- c(20, 50, 30)
index = c(1, 21, 71)

fit <- stan(here::here("STAN", "heights_groups_1.stan"), 
            data = list(N = N, K = 3, Y = Y, S = S, index = index), 
            iter = 1000, 
            chains = 1, 
            control = list(adapt_delta = 0.95, stepsize = 0.01))


mu <- extract(fit, "mu")[[1]] 

print(fit, probs = c(0.25, 0.5, 0.75))

### SOLUTION 2

# An alternative way to estimate this type of model is to pass an array which 
# identifies the group to which each observation belongs:

groups = c(rep(1, 20), rep(2, 50), rep(3, 30))

fit <- stan(here::here("STAN", "heights_groups_2.stan"), 
            data = list(N = N, K = 3, Y = Y, groups = groups), 
            iter = 1000, 
            chains = 1, 
            control = list(adapt_delta = 0.95, stepsize = 0.01))


mu <- extract(fit, "mu")[[1]] 

print(fit, probs = c(0.25, 0.5, 0.75))


# WAIC, LOO-CV and other measures -----------------------------------------

# Often we want to estimate the predictive performance of a model to compare it 
# with others in order to choose between competing hypotheses. Where possible, 
# researchers should repeatedly partition their data into training and test sets 
# (that is, use explicit cross-validation). The training sets are used to fit 
# the model, and the test sets to estimate their out-of-sample predictive 
# capability. However, there are circumstances where repeated partitioning is 
# not feasible due to the computational cost of estimating a model. In these 
# cases, we can use WAIC and estimates of LOO-CV to measure the out-of-sample 
# predictive capability of a model. 
# Suppose that we generate some fake data from a Student-t distribution
N <- 10000
X <- rt(N, 5)

# We then fit two models to the data – one uses a normal sampling distribution 
# and the other assumes a Student-t sampling distribution. Here we know the 
# Student-t distribution should perform better since we used it to generate our 
# data, and hence we use this toy problem to illustrate how cross-validation, 
# WAIC and LOO-CV can be calculated.

fit_normal <- stan(here::here("STAN", "toy_normal.stan"), 
                   data = list(N = N, X = X), 
                   iter = 1000, 
                   chains = 1, 
                   control = list(adapt_delta = 0.95, stepsize = 0.01))

fit_student <- stan(here::here("STAN", "toy_student.stan"), 
                    data = list(N = N, X = X), 
                    iter = 1000, 
                    chains = 1, 
                    control = list(adapt_delta = 0.95, stepsize = 0.01))

library(loo)

logLikelihood_normal <- extract_log_lik(fit_normal, "logLikelihood")
WAIC_normal <- waic(logLikelihood_normal)

logLikelihood_student <- extract_log_lik(fit_student, "logLikelihood")
WAIC_student <- waic(logLikelihood_student)

# We see that the Student-t model has a higher estimated expected log pointwise 
# predictive density (elpd = –16,467.8, which, times –2, corresponds to a lower 
# WAIC). However, to determine whether this difference represents anything other 
# than sampling error, we compare these two models in the correct – pairwise – 
# way using:

compare(WAIC_normal, WAIC_student)
loo_compare(WAIC_normal, WAIC_student)

# we see that the difference in elpd is much greater than the standard error. 
# This suggests that there is a significant difference in performance between 
# these two models, and we prefer the Student-t model. But how do we determine 
# whether the difference is significant? One way is to calculate a z score and 
# compare with a standard normal: pr(z >=  516.2 / 73.4 )

# While we are not in general fond of Frequentist hypothesis tests, this is the 
# current state of the art here. In this case, it indicates that there is 
# basically zero probability of this difference occurring if both models are 
# equally predictive.

# The same package allows us to estimate elpd via LOO-CV without doing explicit 
# cross-validation.

LOO_normal <- loo(logLikelihood_normal)
LOO_student <- loo(logLikelihood_student)

compare(LOO_normal, LOO_student)
loo_compare(LOO_normal, LOO_student)

# This produces similar estimates of elpd as for the WAIC in this case. In both 
# cases, we find that the estimated elpd is greater for the Student-t model than 
# for the normal one, as expected. In general, we prefer the LOO-CV measure 
# since it represents a better approximation to the out-of-sample predictive 
# capability of our model. However, it is important to note that it is not 
# uncommon to get warnings about either the ‘p_waic exceeding ...’ or ‘Pareto k 
# estimates exceeding...’. While we refer the interested reader to the details 
# of the loo paper itself [40], we mention here that it is important to take 
# heed of these warnings. These warnings typically indicate that one or more of 
# the approximations used to estimate these criteria are likely violated, and 
# hence inferences about model performance cannot be trusted. In these cases, it 
# may be better to use explicit cross-validation, to which we turn our attention 
# now.


# EXPLICIT CROSS VALIDATION -----------------------------------------------

library(caret)

testIndices <- createFolds(X, k = 5, list = TRUE, returnTrain = FALSE)

kFoldCV <- function(aModel, TestIndices, X) {
  numFolds <-  length(TestIndices)
  
  # Calculate expected log pointwise predictive density
  lPointLoglikelihoodTotal <- vector()
  
  for (i in 1:numFolds) {
    print(paste0("kFold: ", i))
    XTest <-  X[ TestIndices[[i]]]
    XTrain <- X[-TestIndices[[i]]]
    
    fit <- sampling(aModel, 
                    data = list(NTest = length(XTest), 
                                XTest = XTest,
                                NTrain = length(XTrain),  
                                XTrain = XTrain))
    logLikelihood <- extract_log_lik(fit, "logLikelihood")
    
    lPointLoglikelihood <- colMeans(logLikelihood)
    lPointLoglikelihoodTotal <- c(lPointLoglikelihood,
                                  lPointLoglikelihoodTotal)
  }
  
  return(lPointLoglikelihoodTotal)
}

model_CV_normal  <- stan_model(here::here("STAN", "toy_CV_normal.stan"))
model_CV_student <- stan_model(here::here("STAN", "toy_CV_student.stan"))

lELPD_normal  <- kFoldCV(model_CV_normal,  testIndices, X)
lELPD_student <- kFoldCV(model_CV_student, testIndices, X)

sum(lELPD_normal)
sum(lELPD_student)

difference <- sum(lELPD_student) - sum(lELPD_normal)

sd <- sqrt(1000)*sd(lELPD_student-lELPD_normal)
pvalue <- 1 - pnorm(difference/sd)

pvalue

# noting that again the Student-t distribution performs better than the normal 
# distribution. (Above, we calculate the standard deviation on the difference 
# in log pointwise predictive density to calculate a z score.)


# TIDYBAYES ---------------------------------------------------------------

library(tidybayes) # https://mjskay.github.io/tidybayes/

# fit de toy_student.stan

fit %>%
  spread_draws(mu[S], sigma[S])

fit %>%
  spread_draws(mu[S], sigma[S]) %>% 
  ggplot(aes(x = mu, y = S)) +
  stat_eye()

fit %>%
  spread_draws(mu[S], sigma[S]) %>% 
  ggplot(aes(x = mu, y = S)) +
  stat_halfeyeh()

fit %>%
  spread_draws(mu[S], sigma[S]) %>% 
  ggplot(aes(x = mu, y = S)) +
  stat_dots()

fit %>%
  spread_draws(mu[S]) %>% 
  median_qi(mu)

fit %>%
  spread_draws(mu[S]) %>% 
  mean_qi(mu)
