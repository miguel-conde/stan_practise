library(tidyverse)

N <- 1000
Z <- rnorm(N)

X <- rnorm(N, Z, 1.5)

Y <- rnorm(N, -3*Z, 2)

df <- tibble(X, Y, Z)

# X e Y son dependientes, pero independientes condicionalemnte en Z
# (Z es un fork)

# Por eso X e Y están correladas, aunque Y no dependa de X
cor(df)

# Modelo sin condicionar en Z. Al ser X e Y dependientes (ya que el fork queda
# abierto), el coef de X debería ser significativo
lm_X_Y <- lm(Y ~ X, df)

summary(lm_X_Y)

# Modelo condicionando en Z. Al ser X e Y independientes (ya que hemos cerrado 
# fork al condicionar en Z), el coef de X NO debería ser significativo
lm_X_Y_Z <- lm(Y ~ X + Z, df)

summary(lm_X_Y_Z)

# Sin embargo, SI ADEMÁS HUBIERA UNA RELACIÓN CAUSAL DE X sobre Y, aunque
# cerremos el fork el coeficiente de X SI sería significativo.
Y <- rnorm(N, 2*X - 3*Z, 2)

df_causal <- tibble(X, Y, Z)
lm_X_Y_Z_causal <- lm(Y ~ X + Z, df_causal)

summary(lm_X_Y_Z_causal)

## Esto nos permite comprobar a partir de  los datos si el DAG de nuestro modelo
## es correcto