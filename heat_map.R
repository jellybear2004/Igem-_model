# === Load required libraries ===
packages <- c("ggplot2", "dplyr", "tidyr")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(packages, library, character.only = TRUE)

# === 1. Simulation function === 
simulate_CM1 <- function(X1_init, X2_init, dt = 0.001, t_end = 3) {
  # === Phase C: formation du curli ===
  dt_C <- 0.001
  t_end_C <- 5.0
  times_C <- seq(0, t_end_C, by = dt_C)
  n_steps_C <- length(times_C)

  # Constantes
  kf <- 760  # L/(mol*s)
  masse_molaire_P1 <- 14000
  masse_molaire_P2 <- 12460
  masse_molaire_C <- masse_molaire_P1 + masse_molaire_P2
  alpha1 <-0.40
  alpha2 <-0.00252

  # Initialisation
  dfC <- data.frame(
    time = times_C,
    P1 = numeric(n_steps_C),
    P2 = numeric(n_steps_C),
    C  = numeric(n_steps_C)
  )
  dfC$P1[1] <- (X1_init*alpha1) / masse_molaire_P1
  dfC$P2[1] <- (X2_init*alpha2) / masse_molaire_P2
  dfC$C[1]  <- 0

  # Simulation Phase C
  for (i in 1:(n_steps_C - 1)) {
    P1 <- dfC$P1[i]
    P2 <- dfC$P2[i]
    C  <- dfC$C[i]

    dP1 <- -kf * P1 * P2
    dP2 <- -kf * P1 * P2
    dC  <-  kf * P1 * P2

    dfC$P1[i + 1] <- max(P1 + dP1 * dt_C, 0)
    dfC$P2[i + 1] <- max(P2 + dP2 * dt_C, 0)
    dfC$C[i + 1]  <- max(C  + dC  * dt_C, 0)
  }

  # === Phase D: complexation du curli avec métal M1 ===
  dt_D <- dt
  t_end_D <- t_end
  times_D <- seq(0, t_end_D, by = dt_D)
  n_steps_D <- length(times_D)

  # Paramètres cinétiques
  kon_1  <- 1e9       # L/(mol*s) — réduit pour stabilité numérique
  koff_1 <- 0.05      # 1/s
  gamma  <- 4         # facteur correction

  # Initialisation
  df2 <- data.frame(
    time = times_D,
    C    = numeric(n_steps_D),
    CM1  = numeric(n_steps_D),
    M1   = numeric(n_steps_D)
  )

  # Récupérer valeur finale de C
  idx_end_C <- which.min(abs(dfC$time - t_end_C))
  df2$C[1] <- dfC$C[idx_end_C]


  df2$M1[1]  <- 0.001  # mol/L
  df2$CM1[1] <- 0.0

  # Simulation Phase D
  for (i in 1:(n_steps_D - 1)) {
    C   <- df2$C[i]
    M1  <- df2$M1[i]
    CM1 <- df2$CM1[i]

    dM1  <- koff_1 * CM1 - kon_1 * C * M1 * gamma 
    dCM1 <- -dM1  # conservation
    dC   <- koff_1 * CM1 / gamma - kon_1 * C * M1 

    df2$M1[i + 1]  <- max(M1  + dM1  * dt_D, 0)
    df2$CM1[i + 1] <- max(CM1 + dCM1 * dt_D, 0)
    df2$C[i + 1]   <- max(C    + dC   * dt_D, 0)
  }

  # Retourne la valeur finale de CM1
  return(df2$CM1[n_steps_D])
}

# === 2. Exécuter la grille de simulations ===
P1_vals <- seq(0.5, 10, by = 0.5)
P2_vals <- seq(0.5, 10, by = 0.5)

results <- expand.grid(P1 = P1_vals, P2 = P2_vals) %>%
  rowwise() %>%
  mutate(CM1_final = simulate_CM1(P1, P2, dt = 0.001, t_end = 1)) %>%
  ungroup()

# Nettoyer les NA (au cas où)
results <- results %>% filter(!is.na(CM1_final))

# === 3. Affichage : Heatmap ===
ggplot(results, aes(x = P1, y = P2, fill = CM1_final)) +
  geom_tile() +
  scale_fill_viridis_c(name = "final CM1 [mol/L]") +
  labs(
    title = "Heatmap: CM1 vs P1 & P2 (0.5–5 mol/L)",
    x = "X1  (g/L)",
    y = "X2  (g/L)"
  ) +
  theme_minimal()
