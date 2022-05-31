library(bayesrules)
library(tidyverse)
library(lme4)

# Load data
data(cherry_blossom_sample)
running <- cherry_blossom_sample

# Remove NAs
running <- running %>% 
  select(runner, age, net) %>% 
  na.omit()


# POOLED ------------------------------------------------------------------

pooled_lm <- lm(net ~ age, running)
summary(pooled_lm)

ggplot(running, aes(age, net)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)



# NO POOLING --------------------------------------------------------------

no_pooling_lm <- lm(net ~ age:runner, running)
summary(no_pooling_lm)

running %>% filter(runner %in% 1:10) %>% 
  ggplot(aes(age, net, color = runner)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)


# MULTILEVEL --------------------------------------------------------------

multilvl_lm <- lmer(net ~ age + (age | runner), running)
coef(multilvl_lm)


coef_multilvl_lm <- coef(multilvl_lm)$runner %>% 
  as_tibble() %>% 
  rename(intercept = `(Intercept)`, slope = age) %>% 
  mutate(runner = factor(1:nrow(.)), .before = 1)


running %>% filter(runner %in% 1:10) %>% 
  ggplot(aes(age, net)) +
  ## LM POOLED
  geom_point() +
  geom_smooth(data = running, mapping = aes(age, net), color = "black", 
              size = 2, method = "lm", se = FALSE) +
  ## LM NO POOLING
  geom_point(aes(color = runner)) +
  geom_smooth(aes(color = runner), method = "lm", se = FALSE) +  
  ## LM MULTILEVEL
  # FIXED + RANDOM EFFECTS
  geom_abline(data = coef_multilvl_lm  %>% slice(1:10), 
              mapping = aes(intercept = intercept, slope = slope, color = runner), 
              linetype = 2, size = 1.2) +
  # FIXED EFFECTS (AVERAGE REGRESSION)
  geom_abline(data = coef_multilvl_lm %>% 
                summarise_if(is.numeric, mean), 
              mapping = aes(intercept = intercept, slope = slope), 
              color = "black",
              linetype = 2, size = 2) +
  xlim(40, NA)

fixef(multilvl_lm)
coef_multilvl_lm %>% 
  summarise_if(is.numeric, mean)
