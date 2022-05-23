library(bayesrules)
library(tidyverse)


## PRIOR

library(manipulate)

manipulate(
  curve(dbeta(x, alpha, beta)),
  alpha = slider(.1, 5, initial = 1, step = .1),
  beta = slider(.1, 5, initial = 1, step = .1)
        )

# As such, you've conducted 30 different polls throughout the election season. 
# Though Michelle's support has hovered around 45%, she polled at around 35% in 
# the dreariest days and around 55% in the best days on the campaign trail
# avg = alpha / (alpha + beta) = .45 => alpha / beta ~ 9 / 11
plot_beta( 9, 11, mean = TRUE, mode = TRUE)
plot_beta(27, 33, mean = TRUE, mode = TRUE)
plot_beta(45, 55, mean = TRUE, mode = TRUE)

rbeta(1000, 45, 55) %>% mean()
rbeta(1000, 45, 55) %>% var()
rbeta(1000, 45, 55) %>% sd()

## LIKELIHOOD
