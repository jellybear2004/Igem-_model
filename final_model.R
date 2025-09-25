---
title: "final?"
author: "moi"
date: "`r Sys.Date()`"
output: word_document
---

```{r parameters}
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

convert = function(OD){
  biomass = OD*0.04 # based on OD calculator : https://www.calculatorsconversion.com/en/optical-density-od%e2%82%86%e2%82%80%e2%82%80-calculator/
  return(biomass)
}
```

```{r}
sequence = (seq(0,2,1))
a_initials = sequence
b_initials = sequence
for(a in 1:length(sequence)){
  var_CM1_BP = numeric(length(sequence)) 
  var_X2 = numeric(length(sequence))
  for(b in 1:length(sequence)){
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
    dfA$X1[1] <- a_initials[a] #convert(0.002)
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
     mu=numeric(n_steps_B),     # 
     P2 = numeric(n_steps_B)     # protein P2 [g/L] = binding protein
    )

    # Initial conditions
    dfB$X2[1] <- b_initials[b] #convert(0.002) 
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

    #environnement C
    # Time settings
    dt_C <- 0.001
    t_end_C <- 3.0
    times_C <- seq(0, t_end_C, by = dt_C)
    n_steps_C <- length(times_C)

    # Create storage dataframe

    dfC = data.frame(
     time = times_C,
     X1 = numeric(n_steps_C),
     P1 = numeric(n_steps_C),
     P2 = numeric(n_steps_C),
     S3 = numeric(n_steps_C),
     C = numeric(n_steps_C)
    )


    # Initial conditions
    dfC$P1[1] <- dfA$P1[dfA$time==48]/masse_molaire_P1#*(2*10^-6)/(1*10^-3)
    dfC$P2[1] <- dfB$P2[dfB$time==24]/masse_molaire_P2 # il faut encore prendre en compte les différences de  volumes, mais les résultats changent de fou!
    dfC$C[1] = 0


    # Simulation loop
    for (i in 1:(n_steps_C - 1)) {
     # Extract current values
     P1 <- dfC$P1[i]
     P2 <- dfC$P2[i]
     C = dfC$C[i]
  
     # --- reactions ---
        # Complex formation
        dP1 = - kf*P1*P2
        dP2 = - kf*P1*P2
    
        dC = kf*P1*P2
      
      # --- Euler update ---
      dfC$P1[i+1]  <- max(P1 + dP1 * dt_C, 0) # mol/L
      dfC$P2[i+1]  <- max(P2 + dP2 * dt_C, 0) # mol/L
      dfC$C[i+1]   <- max(C + dC * dt_C, 0)   # mol/L
    }

    #environnement D
    # Time settings
    dt_D <-     0.0001
    t_end_D <-  0.5
    times_D <- seq(0, t_end_D, by = dt_D)
    n_steps_D <- length(times_D)

    # Create storage dataframe

    dfD = data.frame(
      time = times_D,
      C = numeric(n_steps_D),
      CM1 = numeric(n_steps_D),    # complex-metal1 [mol/L] curli-BP-metal1
      CM2 = numeric(n_steps_D),    # complex-metal2 [mol/L] curli-BP-metal2
      M1 = numeric(n_steps_D),     # free metal 1 [mol/L]
      M2 = numeric(n_steps_D)     # free metal 2 [mol/L]
    )

    # Initial conditions
    dfD$C[1] = dfC$C[dfC$time==3]
    dfD$M1[1] <- 0.001
    #dfD$M2[1] <- 0.000
    dfD$CM1[1] <- 0.0
    #dfD$CM2[1] <- 0.0

    # Simulation loop
    for (i in 1:(n_steps_D - 1)) {
      # Extract current values
      C = dfD$C[i]
      M1 <- dfD$M1[i]
      #M2 <- dfD$M2[i]
      CM1 <- dfD$CM1[i]
      #CM2 <- dfD$CM2[i]

      # --- reactions ---
        # Complex formation

        dM1 = koff_1*CM1 - kon_1*C*M1*gamma*(M1/(M1+aff_1))
        #dM2 = koff_2*CM2 - kon_2*C*M2*gamma*(M2/(M2+aff_2))
    
        dCM1 = - koff_1*CM1 + kon_1*C*M1*gamma*(M1/(M1+aff_1))
        #dCM2 = - koff_2*CM2 + kon_2*C*M2*gamma*(M2/(M2+aff_2))
    
        dC = koff_1*CM1/gamma - kon_1*C*M1*(M1/(M1+aff_1)) #+ koff_2*CM2/gamma - kon_2*C*M2*(M2/(M2+aff_2))
      
      # --- Euler update ---
      dfD$M1[i+1]  <- max(M1 + dM1 * dt_D, 0)
      #dfD$M2[i+1]  <- max(M2 + dM2 * dt_D, 0)
      dfD$C[i+1]   <- max(C + dC * dt_D, 0)
      dfD$CM1[i+1] <- max(CM1 + dCM1 * dt_D, 0)
      #dfD$CM2[i+1] <- max(CM2 + dCM2 * dt_D, 0)
    }
    var_CM1_BP[b] = dfD$CM1[5001]
    var_X2[b] = b_initials[b]+a_initials[a]
  }
  plot(var_X2, var_CM1_BP)
  #print(dfD$CM1[5001]/(dfA$X1[1] + dfB$X2[1]))
}
```

```{r graphs environnement D}
# --- Group variables by environment ---
envD_vars <- c("C", "M1", "CM1")#, "M2", "CM2")  # Post-mixing reactions

# --- Plotting function with margin fixes ---
plot_group <- function(vars, group_title) {
  n <- length(vars)
  n_row <- ceiling(sqrt(n))
  n_col <- ceiling(n / n_row)

  # Adjust margins and text sizes for readability
#  par(mfrow = c(n_row, n_col),        # subplot layout
#      mar = c(4, 4, 2, 1),            # inner margins (bottom, left, top, right)
#      oma = c(0, 0, 4, 0))            # outer margins (for main title)
  par(mfrow = c(2,3))

  for (var in vars) {
    plot(dfD$time, dfD[[var]],
         type = "l",
         col = "blue",
         xlab = "Time [h]",
         ylab = paste0(var, " [", ifelse(grepl("^M", var), "mol/L", "mol/L"), "]"),
         main = var,
         cex.main = 1.3,    # main title size
         cex.lab = 1.1,     # axis label size
         cex.axis = 0.9,    # tick size
         las = 1)           # horizontal y-axis labels
  }

  # Add main title centered above all subplots
  mtext(group_title, outer = TRUE, line = 2, cex = 1.5, font = 2)
  par(mfrow = c(1, 1))  # reset layout
}

# --- Plot all environments ---
plot_group(envD_vars, "Environment D: Complex and Metal Binding")
```

```{r}
plot(dfA$X1[1] + dfB$X2[1], dfD$CM1[5001])
print(dfD$CM1[5001]/(dfA$X1[1] + dfB$X2[1]))
```

#temps pour bind un pourcentage des métaux
```{r}
percent = 90
temps=0
for (ion in dfD$M1){
  if (ion <= dfD$M1[1]*(1-percent/100)){
    if (temps==0){
      temps=ion
    }
  }
}
dfD$time[dfD$M1==temps]
#dfD$CM1[dfD$M1==temps]
#dfD$CM1[dfD$M1==temps] + dfD$C[dfD$M1==temps]

percent2 = 90
temps2=0
for (ion in dfD$M2){
  if (ion <= dfD$M2[1]*(1-percent2/100)){
    if (temps2==0){
      temps2=ion
    }
  }
}
dfD$time[dfD$M2==temps2]
#dfD$CM2[dfD$M2==temps2]
#dfD$CM2[dfD$M2==temps2] + dfD$C[dfD$M2==temps2]
```
