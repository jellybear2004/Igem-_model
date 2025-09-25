{r}
# === Load required libraries ===
packages <- c("ggplot2", "dplyr", "tidyr")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(packages, library, character.only = TRUE)

# === 1. Simulation function === 
simulate_CM1 <- function(X1_init, X2_init, dt = 0.001, t_end = 3) {
  
  mu_max1 <- 1.1
  Ks1 <- 0.7
  Y1 <- 0.8
  alpha1 <- 0.40
  
  mu_max2 <- 1.99
  Ks2 <- 0.486
  Y2 <- 0.63
  alpha2 <- 0.00252
  
  kf <- 760 # doit être en L/(mol*s), voir papier bench : en L/(mol*h) = 2736000
  masse_molaire_P1 = 14000 # Da => to check !!!
  masse_molaire_P2 =  12460 # Da => to check !!!
  masse_molaire_C = masse_molaire_P1+masse_molaire_P2
  
  kon_1 <- 10^9 
  koff_1 <- 0.05 
  aff_1 <- 10^-6
  
  gamma = 4

  #Environement A
  # Time settings
  dt_A <- 0.01
  t_end_A <- 48
  times_A <- seq(0, t_end_A, by = dt_A)
  n_steps_A <- length(times_A)

  # Create storage dataframe
  dfA <- data.frame(
    time = times_A,
    X1 = numeric(n_steps_A),    # biomass A [g/L] (nombre de bactéries)
    S1 = numeric(n_steps_A),    # substrate A [g/L]
    P1 = numeric(n_steps_A)    # protein P1 [g/L] = curli
  )

  # Initial conditions
  dfA$X1[1] <- X1_init 
  dfA$S1[1] <- 2
  dfA$P1[1] <- 0.0

  # Simulation loop
  for (i in 1:(n_steps_A - 1)) {
    # Extract current values
    X1 <- dfA$X1[i] 
    S1 <- dfA$S1[i]
    P1 <- dfA$P1[i]


   # --- Growth ---
   mu1 <- mu_max1 * (S1 / (Ks1 + S1))
    dX1 <-  mu1 * X1 
    dS1 <- (-1/Y1 * mu1) * X1

    # --- Protein production ---
    if (dfA$time[i]>2){
     if (X1-dfA$X1[i-1]==dfA$X1[i-1]-dfA$X1[i-2]){
       dP1 <- alpha1 * X1
     }
     else{
        dP1 = 0
      }
    }
   else{
     dP1 = 0
   }

    # --- Euler update ---
   dfA$X1[i+1]    <- max(X1 + dX1 * dt_A,0)
    dfA$S1[i+1]    <- max(S1 + dS1 * dt_A,0)
    dfA$P1[i+1]    <- max(P1 + dP1 * dt_A,0)
  }

  # Environnement B
  # Time settings
  dt_B <- 0.01
  t_end_B <- 24
  times_B <- seq(0, t_end_B, by = dt_B)
  n_steps_B <- length(times_B)

  # Create storage dataframe
  dfB = data.frame(
    time = times_B,
   X2 = numeric(n_steps_B),     # biomass B [g/L]
   S2 = numeric(n_steps_B),     # substrate B [g/L]
   mu = numeric(n_steps_B),     # 
   P2 = numeric(n_steps_B)     # protein P2 [g/L] = binding protein
  )

  # Initial conditions
  dfB$X2[1] <- X2_init 
  dfB$S2[1] <- 2
  dfB$P2[1] <- 0.0


  # Simulation loop
  for (i in 1:(n_steps_B - 1)) {
    # Extract current values
    X2 <- dfB$X2[i]
    S2 <- dfB$S2[i] 
    P2 <- dfB$P2[i]

    # --- Growth ---
    dfB$mu[i] <- mu_max2 * (S2 / (Ks2 + S2))
    mu2=dfB$mu[i]
    dX2 <-  mu2 * X2 
    dS2 <- (-1/Y2 * mu2) * X2

    # --- Protein production ---
    #dP2 <- if (mu2 > mu_min) alpha2 * X2 - delta2 * P2 else -delta2 * P2
    if (dfB$time[i] >2){
     dP2 <- if (X2-dfB$X2[i-1]-(dfB$X2[i-1]-dfB$X2[i-2])>0) alpha2 * X2  else 0
    }
    else{
      dP2 = alpha2 * X2 
    }

    # --- Euler update ---
    dfB$X2[i+1]    <- X2 + dX2 * dt_B
    dfB$S2[i+1]    <- S2 + dS2 * dt_B
    dfB$P2[i+1]    <- P2 + dP2 * dt_B
  }
  
  # === Phase C: formation du curli ===
  dt_C <- 0.001
  t_end_C <- 3.0
  times_C <- seq(0, t_end_C, by = dt_C)
  n_steps_C <- length(times_C)

  # Constantes
  

  # Initialisation
  dfC <- data.frame(
    time = times_C,
    P1 = numeric(n_steps_C),
    P2 = numeric(n_steps_C),
    C  = numeric(n_steps_C)
  )
  
  dfC$P1[1] <- dfA$P1[n_steps_A] / masse_molaire_P1
  dfC$P2[1] <- dfB$P2[n_steps_B] / masse_molaire_P2
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

    dM1  <- koff_1 * CM1 - kon_1 * C * M1 * gamma * (M1 / (M1 + aff_1))
    dCM1 <- -dM1  # conservation
    dC   <- koff_1 * CM1 / gamma - kon_1 * C * M1 * (M1 / (M1 + aff_1))

    df2$M1[i + 1]  <- max(M1  + dM1  * dt_D, 0)
    df2$CM1[i + 1] <- max(CM1 + dCM1 * dt_D, 0)
    df2$C[i + 1]   <- max(C    + dC   * dt_D, 0)
  }

  # Retourne la valeur finale de CM1
  return(df2$CM1[n_steps_D])
}

# === 2. Exécuter la grille de simulations ===
P1_vals <- seq(1, 10, by = 1)
P2_vals <- seq(9990, 10000, by = 1)

results <- expand.grid(P1 = P1_vals, P2 = P2_vals) %>%
  rowwise() %>%
  mutate(CM1_final = simulate_CM1(P1, P2, dt = 0.001, t_end = 0.5)) %>%
  ungroup()

# Nettoyer les NA (au cas où)
results <- results %>% filter(!is.na(CM1_final))

# === 3. Affichage : Heatmap ===
ggplot(results, aes(x = P1, y = P2, fill = CM1_final)) +
  geom_tile() +
  scale_fill_viridis_c(name = "final CM1 [mol/L]") +
  labs(
    title = "Heatmap: CM1 vs P1 & P2 (0.5–5 mol/L)",
    x = "P1 initial (mol/L)",
    y = "P2 initial (mol/L)"
  ) +
  theme_minimal()

