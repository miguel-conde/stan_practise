# Ch17 - (Normal) Hierarchical Models with Predictors

# Load packages
library(bayesrules)
library(tidyverse)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(broom.mixed)


# Data --------------------------------------------------------------------


# Load data
data(cherry_blossom_sample)
running <- cherry_blossom_sample

# Remove NAs
running <- running %>% 
  select(runner, age, net) %>% 
  na.omit()
