library(tidyverse)
library(dagitty)

# Cerrar la puerta de atrás: bloquear los caminos confusores entre un predictor X
# y una respuesta Y

N <- 1000

# 1.1 - DAGs RELACIONES ELEMENTALES ---------------------------------------

# 1.1 - FORK --------------------------------------------------------------


# El clásico confusor. Z es la causa común a X e Y, generando una correlación
# entre ellas. 
dag_fork <- dagitty("dag {X <- Z -> Y}")

Z <- rnorm(N)
X <-  3 * Z + rnorm(N, 0, 0.1)
Y <- -2 * Z + rnorm(N, 0, 0.01)

df_fork <- tibble(X, Y, Z)

cor(df_fork)

# Si condicionamos en Z, añadir X no aporta nada. X e Y son independientes, 
# condicionalmente en Z. Hay una asociación entre X e Y inducida por Z, pero
# no una relación causal de X sobre Y (intervenir en X no cambiará Y)

coordinates(dag_fork) <- list(x = c(X = 0, Z = 1, Y = 2),
                              y = c(X = 1, Z = 0, Y = 1))
drawdag(dag_fork)

Z_fork_lm <- lm(Y ~ Z, df_fork)

summary(Z_fork_lm)

Z_y_X_fork_lm <- lm(Y ~ X + Z, df_fork)

summary(Z_y_X_fork_lm)

cor(predict(Z_fork_lm, df_fork), predict(Z_y_X_fork_lm, df_fork))
plot(predict(Z_fork_lm, df_fork), predict(Z_y_X_fork_lm, df_fork))

# Para cerrar un FORK (que la información no fluya de X a Y, que sean 
# independientes), condicionamos en Z (incluimos z) 

adjustmentSets(dag_fork, exposure = "X", outcome = "Y")
impliedConditionalIndependencies(dag_fork)

X_fork_lm <- lm(Y ~ X, df_fork)

summary(X_fork_lm)

cor(predict(X_fork_lm, df_fork), predict(Z_fork_lm, df_fork))
plot(predict(X_fork_lm, df_fork), predict(Z_fork_lm, df_fork))

cor(predict(X_fork_lm, df_fork), predict(Z_y_X_fork_lm, df_fork))
plot(predict(X_fork_lm, df_fork), predict(Z_y_X_fork_lm, df_fork))


# 1.2 - PIPE --------------------------------------------------------------

dag_pipe <- dagitty("dag {X -> Z -> Y}")

# El tratamiento X influye a Z que influye a la respuesta Y

X <- rnorm(N)
Z <-  3 * X + rnorm(N, 0, 0.1)
Y <- -2 * Z + rnorm(N, 0, 0.01)

df_pipe <- tibble(X, Y, Z)

cor(df_pipe)

# Si condicionamos en Z, también bloqueamos el path entre X e Y (o sea, si
# conocemos Z conocer X tampoco aporta nada, X e Y son independientes)

coordinates(dag_pipe) <- list(x = c(X = 0, Z = 1, Y = 2),
                              y = c(X = 0, Z = 1, Y = 2))
drawdag(dag_pipe)

Z_pipe_lm <- lm(Y ~ Z, df_pipe)

summary(Z_pipe_lm)

Z_y_X_pipe_lm <- lm(Y ~ X + Z, df_pipe)

summary(Z_y_X_pipe_lm)

cor(predict(Z_pipe_lm, df_pipe), predict(Z_y_X_pipe_lm, df_pipe))
plot(predict(Z_pipe_lm, df_pipe), predict(Z_y_X_pipe_lm, df_pipe))

# Para cerrar un PIPE, no condicionamos (no incluimos) Z

adjustmentSets(dag_pipe, exposure = "X", outcome = "Y")
impliedConditionalIndependencies(dag_pipe)

X_pipe_lm <- lm(Y ~ X, df_pipe)

summary(X_pipe_lm)

cor(predict(X_pipe_lm, df_pipe), predict(Z_pipe_lm, df_pipe))
plot(predict(X_pipe_lm, df_pipe), predict(Z_pipe_lm, df_pipe))

cor(predict(X_pipe_lm, df_pipe), predict(Z_y_X_pipe_lm, df_pipe))
plot(predict(X_pipe_lm, df_pipe), predict(Z_y_X_pipe_lm, df_pipe))


# 1.3 - COLLIDER ----------------------------------------------------------
dag_collider <- dagitty("dag {X -> Z <- Y}")

# No hay asociación entre X e Y (X _||_ Y)

X <- rnorm(N)
Y <- rnorm(N)
Z <-  3 * X -2 *Y + rnorm(N, 0, 0.1)

df_collider <- tibble(X, Y, Z)

cor(df_collider)

# A no ser que  condicionemos en Z: Si condicionamos en Z, inducimos una 
# asociación entre X e Y, abrimos el path 
# entre X e Y (o sea, si conocemos Z, X e Y NO son independientes). Pero no hay
# relaciones causales entre X e Y
# 

coordinates(dag_collider) <- list(x = c(X = 0, Z = 1, Y = 2),
                              y = c(X = 0, Z = 1, Y = 0))
drawdag(dag_collider)

Z_collider_lm <- lm(Y ~ Z, df_collider)

summary(Z_collider_lm)

Z_y_X_collider_lm <- lm(Y ~ X + Z, df_collider)

summary(Z_y_X_collider_lm)

cor(predict(Z_collider_lm, df_collider), predict(Z_y_X_collider_lm, df_collider))
plot(predict(Z_collider_lm, df_collider), predict(Z_y_X_collider_lm, df_collider))

# Para mantener cerrado un COLLIDER (X e Y independientes), NO condicionamos 
# (NO incluimos) Z

adjustmentSets(dag_collider, exposure = "X", outcome = "Y")
impliedConditionalIndependencies(dag_collider)

X_collider_lm <- lm(Y ~ X, df_collider)

summary(X_collider_lm)

cor(predict(X_collider_lm, df_collider), predict(Z_collider_lm, df_collider))
plot(predict(X_collider_lm, df_collider), predict(Z_collider_lm, df_collider))

cor(predict(X_collider_lm, df_collider), predict(Z_y_X_collider_lm, df_collider))
plot(predict(X_collider_lm, df_collider), predict(Z_y_X_collider_lm, df_collider))

# 1.4 - DESCENDANT --------------------------------------------------------




# 2 - RECETA --------------------------------------------------------------

# Fork y pipe (en Z): si condiciono en Z los cierro (X e Y independientes si 
#                     condiciono en Z)
# Collider (en Z): si condiciono en Z lo abro (X e Y dependientes si condiciono 
#                  en Z) 

# 1 - Enumerar todos los paths entre la variable causal potencial (X) y la 
# variable respuesta Y.
#
# 2 - Clasificar cada path en ABIERTO o CERRADO:
#     Un path está abierto siempre que no contenga un collider.
#
# 3 - Clasificar cada path según contenga o no una PUERTA TRASERA:
#     Si X tiene una flecha entrante, el path es una puerta trasera
#
# 4 - Para cada path con PUERTA TRASERA ABIERTA:
#     Determonar sobre qué variables condicionar para cerrar la puerta trasera


# 3 - EJEMPLOS ------------------------------------------------------------


# 3.1 - Ejemplo1  ---------------------------------------------------------

dag_2_roads <- dagitty("dag {
  U [unobserved]
  X -> Y
  X <- U <- A -> C -> Y
  U -> B <- C
  }")

coordinates(dag_2_roads) <- list(x = c(U = 0, X = 0, A = 1, B = 1, C = 2, Y = 2),
                                 y = c(U = 2, X = 0, A = 3, B = 1, C = 2, Y = 0))
drawdag(dag_2_roads)

A <- rnorm(N)
U <-  2 * A + rnorm(N, 0, 0.1)
C <- -3 * A + rnorm(N, 0.01)
B <- U + C + rnorm(N, 0, 0.1)
X <- 2 * U + rnorm(N, 0, 0.1)
Y <- -X + C + rnorm(N, 0, 0.1)

df_2_roads <- tibble(A, U, C, B, X, Y)
cor(df_2_roads)

# Paso a paso
#
# x <- U -> B <- C -> Y 
#                        Puerta trasera pq flecha entra en X
#                        Contiene collider => está cerrada
#                        No condicionar en B pq lo abriría
# X <- U <- A -> C -> Y 
#                        Puerta trasera pq flecha entra en X
#                        No contiene collider => está abierta
#                        Hay 1 fork (A) y 2 pipes (U, C)
#                        U => no observada, no podemos condicionar en 
#                             ella
#                        A, C => condicionar en 1 de ellas cierra 
#                                esta puerta
adjustmentSets(dag_2_roads, exposure = "X", outcome = "Y")
impliedConditionalIndependencies(dag_2_roads)

# Mejor condicionar en C => mejor eficiencia por más cercana a Y:
lm_2_roads <- lm(Y ~ X + C, df_2_roads)
summary(lm_2_roads)


# 3.2 - Ejemplo2 ----------------------------------------------------------

library(rethinking)

data(WaffleDivorce)
d <- WaffleDivorce %>% as_tibble()

dag_divorce <- dagitty("dag {
  MedianAgeMarriage -> Divorce
  MedianAgeMarriage -> Marriage -> Divorce
  MedianAgeMarriage <- Location -> Marriage
  Location -> WaffleHouses -> Divorce
  }")

coordinates(dag_divorce) <- list(x = c(MedianAgeMarriage = 0, 
                                       Location = 0, 
                                       Marriage = 1, 
                                       Divorce = 2, 
                                       WaffleHouses = 2),
                                 y = c(MedianAgeMarriage = 2, 
                                       Location = 0, 
                                       Marriage = 1, 
                                       Divorce = 2, 
                                       WaffleHouses = 0))
drawdag(dag_divorce)

# Backdoor Paths
#        WaffleHouses <- Location -> MedianAgeMarriage -> Divorce
#           Abierto por el fork (Location) y el pipe MedianAgeMarriage 
#           Podemos cerrarlo condicionando Location o MedianAgeMarriage
#        WaffleHouses <- Location -> Marriage -> Divorce
#           Abierto por el fork (Location) y el pipe Marriage 
#           Podemos cerrarlo condicionando Location o Marriage
#        WaffleHouses <- Location -> MedianAgeMarriage -> Marriage -> Divorce
#           Abierto por el fork (Location) y el pipe MedianAgeMarriage 
#           Podemos cerrarlo condicionando Location o MedianAgeMarriage
#        WaffleHouses <- Location -> Marriage -> MedianAgeMarriage -> Divorce
#           Cerrado por el collider (Marriage) 

impliedConditionalIndependencies(dag_divorce)
adjustmentSets(dag_divorce, exposure = "WaffleHouses", outcome = "Divorce")

divorce_lm <- lm(Divorce ~ Location, d)

summary(divorce_lm)
