# SIMULACIÓN DE CÓPULA DE CLAYTON Y SU CORRESPONDIENTE CÓPULA DE SUPERVIVENCIA
# PROYECTO FINAL ESTADÍSTICA APLICADA 3 - Otoño 2024
# Fernanda Arelle y Dara Meneses


# Paquetes necesarios
library(copula)
library(ggplot2)
library(viridis)
library(gridExtra)
library(MASS)

# Parámetros para las tasas de falla y censura
beta1 <- -0.5
beta2 <- 0.1
gamma1 <- 0.3
gamma2 <- 0.2

# Parámetros para la cópula de Clayton
alpha <- 8
clayton_cop <- claytonCopula(alpha)

# Simulación de datos
set.seed(123)
n <- 200
copula_samples <- rCopula(n, clayton_cop)
u_star <- copula_samples[, 1]
v_star <- copula_samples[, 2]

# Funciones inversas de supervivencia basándose en las marginales exponenciales
quantile_failure <- function(u, Tr, Age) {
  -log(1 - u) / (0.5 * exp(beta1 * Tr + beta2 * Age))
}

quantile_censor <- function(v, Tr, Age) {
  -log(1 - v) / (0.2 * exp(gamma1 * Tr + gamma2 * Age))
}

# Simulación de covariables
set.seed(123)
Tr <- rbinom(n, size = 1, prob = 0.5) 
Age <- runif(n, min = -10, max = 10)

# Tiempos marginales usando las funciones cuantiles
T <- quantile_failure(u_star, Tr, Age)  # Tiempos de falla
C <- quantile_censor(v_star, Tr, Age)  # Tiempos de censura dependiente
C_admin <- rexp(n, rate = 0.05)        # Censura administrativa

# Supuesto de que en promedio 45% de las observaciones fueron fallas y 35% fueron dadas de baja y el resto fueron censurados administrativamente
set.seed(123)
event_type <- sample(c("failure", "withdrawal", "admin_censor"),
                     size = n,
                     replace = TRUE,
                     prob = c(0.45, 0.35, 0.20))

# Tiempos observados e indicador de censura
T_obs <- numeric(n)
delta <- numeric(n)

for (i in 1:n) {
  if (event_type[i] == "failure") {
    T_obs[i] <- min(T[i], C[i])
    delta[i] <- as.numeric(T[i] <= C[i])  
  } else if (event_type[i] == "withdrawal") {
    T_obs[i] <- C[i]  
    delta[i] <- 0
  } else if (event_type[i] == "admin_censor") {
    T_obs[i] <- min(T[i], C_admin[i])  
    delta[i] <- as.numeric(T[i] <= C_admin[i])
  }
}

# Densidad: cópula de Clayton
density_clayton <- kde2d(u_star, v_star, n = 100)
df_clayton <- data.frame(x = rep(density_clayton$x, each = length(density_clayton$y)),
                         y = rep(density_clayton$y, length(density_clayton$x)),
                         z = as.vector(density_clayton$z))

# Densidad: cópula de supervivencia
u_surv <- 1 - u_star
v_surv <- 1 - v_star
density_surv <- kde2d(u_surv, v_surv, n = 100)
df_surv <- data.frame(x = rep(density_surv$x, each = length(density_surv$y)),
                      y = rep(density_surv$y, length(density_surv$x)),
                      z = as.vector(density_surv$z))

# Gráfica de la cópula de Clayton
plot_clayton <- ggplot() +
  geom_raster(data = df_clayton, aes(x, y, fill = z), interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(data = df_clayton, aes(x = x, y = y, z = z), linetype = "solid", color = "black") +
  geom_point(data = data.frame(u = u_star, v = v_star), aes(x = u, y = v), 
             color = "white", size = 1.5, shape = 16, alpha = 0.6) +
  labs(x = "u", y = "v", title = "Cópula de Clayton (\u03B1=8)") +
  theme_bw()

# Gráfica de la cópula de supervivencia
plot_surv <- ggplot() +
  geom_raster(data = df_surv, aes(x, y, fill = z), interpolate = TRUE) +
  scale_fill_viridis_c() +
  geom_contour(data = df_surv, aes(x = x, y = y, z = z), linetype = "solid", color = "black") +
  geom_point(data = data.frame(u = u_surv, v = v_surv), aes(x = u, y = v), 
             color = "white", size = 1.5, shape = 16, alpha = 0.6) +
  labs(x = "1-u", y = "1-v", title = "Cópula de Supervivencia (\u03B1=8)") +
  theme_bw()

grid.arrange(plot_clayton, plot_surv, ncol = 2)

# Proporciones resultantes del tipo de observación (fallas, dada de baja, censura admin.)
table(event_type) / n
