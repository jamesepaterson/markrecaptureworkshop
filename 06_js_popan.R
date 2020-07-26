## ---------------------------
##
## Script name: 06_js_popan
##
## Purpose of script: Jolly-Seber POPAN models to estimate population size
##
## Author: Dr. James Paterson
##
## Date Created: 2020-07-26
##
## 
## Email: james.earle.paterson@gmail.com
##
## ---------------------------

## ----loadlibrariesdata-----------------
# install.packages("marked") # first time only
library(marked)

# Load dipper data (included in "marked")
data(dipper, package = "marked")

# Examine data structure
head(dipper)

## ----fitjsmodels-----
# Jolly-Seber models (POPAN formulation) are open population models, and 
# can be used to estimate abundance by including two more parameters than the CJS

# Additional parameters:
# Nsuper (or "superpopulation") = total number of individuals available to enter population throughout study
# pent ("probability of entry") =  the rate at which individuals enter the population from Nsuper (via births and immigration)

# WARNING: there is no adequate GOF tests for Jolly-Seber models. 
# One common method: Test equivalent structure of CJS model with R2ucare (previous tutorials).

# This tests *some* assumptions of Phi and p.
# Jolly-Seber models have an additional assumption:
# marked AND unmarked animals have same p (R2ucare doesn't test this)
# This assumption is required to estimate total abundance (sum of marked and unmarked animals in population)

# First, process data (Notice model = "JS", previous version = "CJS")
dipper.js.proc <- process.data(dipper, 
                               model = "JS", 
                               groups = "sex")

# Second, make design data (from processed data)
dipper.js.ddl <- make.design.data(dipper.js.proc)

fit.js.dipper.models <- function(){
  # Phi formulas
  Phi.dot <- list(formula=~1)
  Phi.time <- list(formula=~time)
  # p formulas
  p.dot <- list(formula=~1)
  # pent formulas. pent estimates MUST SUM to 1 (for each group).
  # This is constained using a Multinomial Logit link
  pent.time <- list(formula=~time)
  pent.sex <- list(formula=~sex)
  pent.dot <- list(formula=~1)
  # Nsuper formulas. Don't confuse "N" from model with predicted population size
  N.sex <- list(formula=~sex)
  N.dot <- list(formula=~1)
  cml <- create.model.list(c("Phi","p", "pent", "N"))
  results <- crm.wrapper(cml, data = dipper.js.proc, ddl = dipper.js.ddl,
                         external = FALSE, accumulate = FALSE, hessian = TRUE)
  
  return(results)
}

# Run function
dipper.js.models <- fit.js.dipper.models()

## ----displaymodeltable---------------------------------------------------
# Display model table
dipper.js.models

## ----predictjollyseber---------------------------------------------------

# Look at estimates of top model (row number on left of model table, or using name)
dipper.js.models[[1]]  # or dipper.js.models[["Phi.dot.p.dot.pent.dot.N.dot"]] or dipper.js.models$Phi.dot.p.dot.pent.dot.N.dot

# The estimates above are not on probability scale (or in individuals for N)
# (e.g. Phi, p on logit scale, pent on mlogit scale)
# Predict (real) values using top model
dipper.js.predicted <- predict(dipper.js.models[[1]]) # [[1]] just calls the model row according to the model table.

# Look at predictions of real parameters
dipper.js.predicted 

## ----derivedn------------------------------------------------------------
# Abundance (N) is derived from the estimated parameters
# We will estimate population size at each time by making a dataframe of estimates and calculating N
# We will use the predicted estimates from the top-performing model (in this case: "dipper.js.predicted")

# NOTE: the below method will have to be adjusted based on your final model and the number of capture events
N.derived <- data.frame(occ = c(1:7), # 7 events
                        Phi = c(rep(dipper.js.predicted$Phi$estimate, 6), NA),   # 6 survival estimates all the same
                        Nsuper = rep(dipper.js.predicted$N$estimate + nrow(dipper), 7), # Nsuper estimate + number of marked animals
                        pent = c(1-sum(dipper.js.predicted$pent$estimate), dipper.js.predicted$pent$estimate)) # Sum of all pent must be 1

# Set-up empty vector for calculating N
N.derived$N <- NA

# The inital population size (N[1]) = Nsuper * (1 - sum(all other pent estimates))
# This is because of the link function for estimating pent.
# The sum of all pent parameters MUST equal 1 (therefore, one less must be estimated)
N.derived$N[1] <- (N.derived$Nsuper[1] * N.derived$pent[1])

# Subsequent population sizes are estimated by calculating surviving individuals (N[t-1] * Phi[t]), and
# Adding new births (Nsuper * pent[t])
for(i in 2:nrow(N.derived)){
  N.derived$N[i] <- (N.derived$N[i-1]*N.derived$Phi[i-1]) + (N.derived$Nsuper[i] * N.derived$pent[i])
}

# Look at what we did
N.derived

## ----rmarkpopanderived---------------------------------
detach("package:marked", unload=TRUE) # many of the function names are the same. unload `marked`

# install.packages("RMark") # first time only
# For RMark to work, you also need mark.exe installed separately
# http://www.phidot.org/software/mark/downloads/
# this may not be easy on a Mac OS (http://www.phidot.org/software/mark/rmark/linux/)
# RMark calls "mark" to do all the work outside of R
library(RMark)

# We will use the same data but will just create the same top model (not all the other subsets)
dipper.rmark.processed <- RMark::process.data(dipper,
                                              model = "POPAN")

# Formulae for model
Phi.dot <- list(formula=~1)
p.dot <- list(formula=~1)
pent.dot <- list(formula=~1)
N.dot <- list(formula=~1)

# The argument names are similar but a little different (notice "POPAN" instead of "js")
dipper.rmark <- mark(dipper, model = "POPAN", 
                     model.parameters = list(Phi = Phi.dot, p= p.dot, 
                                                     pent = pent.dot, N = N.dot),
                     realvcv = TRUE)


# The popan.derived function of RMark estimates N 
# (plus estimates SE and 95% CI using the delta method)
dipper.derived.rmark <- popan.derived(dipper.rmark.processed,
                                      dipper.rmark)$N

## ----lookatrmarkderived--------------------------------------------------
# Look at results
dipper.derived.rmark

## ----comparepopanestimates----------------
library(ggplot2)
ggplot(N.derived, aes(x = occ, y = N)) +
  geom_path() +
  geom_path(data = dipper.derived.rmark, aes(x = Occasion+0.1, y = N), col = "blue") +
  geom_point(shape = 16, size = 3, col = "black") +
  # geom_errorbar(aes(ymin = N.lcl, ymax = N.ucl), width = 0.3) +
  geom_point(data = dipper.derived.rmark, aes(x = Occasion+0.1, y = N), col = "blue") +
  geom_errorbar(data = dipper.derived.rmark, aes(x = Occasion+0.1, ymin = LCL, ymax = UCL), width = 0.2, color = "blue") +
  scale_x_continuous(breaks = c(1:7),
                     labels = c(1:7)) +
  scale_y_continuous(limits = c(0, 126), breaks = c(0, 25, 50, 75, 100, 125)) +
  labs(x = "Capture event", y = "Estimated population size\n(and 95% confidence interval)") +
  theme_classic()

