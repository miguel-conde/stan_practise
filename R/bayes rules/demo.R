library(tidyverse)

N_G <- 4

# set.seed(9999)
set.seed(99999)
conf_grupos <- tibble(G = 1:N_G, 
                      beta_0 = rnorm(N_G, 5, 1),
                      beta_1 = rnorm(N_G, 3, 2),
                      sigma_y = 50)

X <- seq(0, 100, by = 10)

n_j <- sample.int(length(X)-3, N_G) + 3

complete_dataset <- lapply(seq_along(n_j), function(i) {
  X = sort(sample(X, n_j[i]))
  tibble(G = i, I = 1:n_j[i], X)
}) %>% 
  bind_rows() %>% 
  left_join(conf_grupos, by = "G") %>% 
  mutate(y = rnorm(length(X), beta_0 + beta_1 * X, sigma_y),
         G = factor(G))

dataset <- complete_dataset %>% 
  select(G, I, X, y)

ggplot(data = dataset) +
  geom_density(mapping = aes(x = y, color = G))

ggplot() +
  geom_point(data = dataset, mapping = aes(x = X, y = y, group = G, color = G)) +
  geom_smooth(data = dataset, aes(x = X, y = y), 
              color = "black", size = 2, method = "lm", se = FALSE) +
  geom_smooth(data = dataset, mapping = aes(x = X, y = y, group = G, color = G), 
              method = "lm", se = FALSE)
  
library(lem4)
