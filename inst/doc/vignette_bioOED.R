## ------------------------------------------------------------------------
library(bioOED)

## ------------------------------------------------------------------------
inactivation_model <- "Bigelow"

## ------------------------------------------------------------------------
temp_profile <- data.frame(time = c(0, 60),
                           temperature = c(30, 60))

## ------------------------------------------------------------------------
parms_fix <- c(temp_ref = 57.5, N0 = 1e6)
parms <- c(D_R = 3.9,
           z = 4.2
           )

## ------------------------------------------------------------------------
sensitivity <- sensitivity_inactivation(inactivation_model,
                                        parms, temp_profile,
                                        parms_fix)

## ------------------------------------------------------------------------
head(sensitivity)

## ---- fig.width=8, fig.height=6------------------------------------------
plot(sensitivity)

## ------------------------------------------------------------------------
parms_fix <- c(temp_ref = 57.5)
parms <- c(delta_ref = 3.9, z = 4.2, p = 1, N0 = 1e6)
temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))
pars_correlation <- calculate_pars_correlation("Mafart", parms,
                                               temp_profile, parms_fix)

## ------------------------------------------------------------------------
print(pars_correlation)

## ---- fig.width=8--------------------------------------------------------
plot(pars_correlation)

## ------------------------------------------------------------------------
temp_ref0 <- 57
lower <- 50
upper <- 70

## ------------------------------------------------------------------------
inactivation_model <- "Mafart"

## ------------------------------------------------------------------------
temp_profile <- data.frame(time = c(0, 60),
                           temperature = c(30, 60))

## ------------------------------------------------------------------------
parms <- c(delta_ref = 3.9, z = 4.2, p = 1, N0 = 1e6)
parms_fix <- c()

## ------------------------------------------------------------------------
optim_refTemp <- optimize_refTemp(temp_ref0, lower, upper,
                                  inactivation_model, parms, temp_profile,
                                  parms_fix)

## ------------------------------------------------------------------------
print(optim_refTemp$par)

## ------------------------------------------------------------------------
print(optim_refTemp$convergence)

## ------------------------------------------------------------------------
inactivation_model <- "Mafart"

## ------------------------------------------------------------------------
parms_fix <- c(temp_ref = 57.5)
parms <- c(delta_ref = 3.9,
           z = 4.2,
           p = 1,
           N0 = 1e6
           )

## ------------------------------------------------------------------------
temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))

## ------------------------------------------------------------------------
n_points <- 5

## ------------------------------------------------------------------------

opts_global <- list(maxeval=1000,  local_solver=0,
                    local_finish="DHC", local_iterprint=1)

global_OED <- inactivation_OED(inactivation_model, parms, temp_profile,
                               parms_fix, n_points, criteria = "D",
                               opts_global = opts_global)


## ------------------------------------------------------------------------
print(global_OED$optim_times)

## ---- fig.width=8--------------------------------------------------------
plot(global_OED)

## ------------------------------------------------------------------------
time_min <- 8

## ------------------------------------------------------------------------
global_OED_penalty <- inactivation_OED_penalty(inactivation_model, parms,
                                               temp_profile, parms_fix,
                                               n_points, time_min,
                                               criteria = "D", 
                                               opts_global = opts_global)

## ------------------------------------------------------------------------
print(global_OED_penalty$optim_times)

## ---- fig.width=8--------------------------------------------------------
plot(global_OED_penalty)

