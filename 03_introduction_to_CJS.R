# Introduction to Cormack-Jolly-Seber models
# Author: James E Paterson 
# Contact: james.earle.paterson@gmail.com
# Date created: 2020-04-26


## ----loadpackages----------------------
# install.packages("marked") # first time using the package
library(marked)  # For building CJS models

## ----cjsmodels-----------------------------------------------------------
# Load example dipper dataset
data(dipper)

# Look at data
head(dipper)

## ----basiccjs--------------------------
# Create a basic CJS model (constant survival and detection) with: 
cjs.m1 <- crm(dipper)

# Examine model and coefficient estimates
cjs.m1

# Re-fit model with precision estimates
cjs.m1 <- cjs.hessian(cjs.m1)
cjs.m1

## ----inverselogit--------------------------------------------------------
# exp(x)/(1+exp(x)) or R function "plogis"
exp(cjs.m1$results$beta$Phi)/(1+exp(cjs.m1$results$beta$Phi)) # real Phi (survival) estimate by hand
plogis(cjs.m1$results$beta$Phi) # real Phi estimate with plogis

# real estimates of Phi and p with the "predict" function
predict(cjs.m1, 
        newdata = data.frame(sex = c("Female", "Male")),
        se = TRUE)  # In this case there are no groups or covariates in the model, so "newdata" not used

## ----unequaltimes----------------------
cjs.m1.unequaltime <- crm(dipper, 
                          # time.interval vector is the interval between each capture event.
                          # for 7 capture events, there are 6 intervals
                          time.intervals = c(1,3,1,1,4,5)) 
predict(cjs.m1.unequaltime)
# Notice the predicted survival is higher now because the time length between capture events is >= 1

## ----cjsfunctions--------------------------
# It is good practice to follow these steps:
# 1. Process data
# 2. Create design data
# 3. Set-up and execute entire candidate model set
# The advantage is it is easier to add covariates (e.g. time-varying effects from weather or experiments)

# Process data (and set grouping variables)
dipper.proc <- process.data(dipper, 
                            group = "sex") # group variable has to be in the data

# Make design data (from processed data)
dipper.ddl <- make.design.data(dipper.proc)

# Outine formulas for each parameter
Phi.dot <- list(formula=~1)  # ~1 is always a constant (or single estimate)
Phi.sex <- list(formula=~sex) # This formula will have an intercept (for females) and an estimate for the difference between females and males
p.sex <- list(formula=~sex) # Be careful of case-sensitive names. Use the exact group column that was in data

# Make new model (using design data) with constant survival, but different detection probabilities between sexes
cjs.m2 <- crm(dipper.proc, 
              dipper.ddl,
              model.parameters = list(Phi = Phi.dot, 
                                      p = p.sex),
              accumulate = FALSE)
cjs.m2
# Is cjs.m2 a better fit than Phi.dot.p.dot (i.e. "cjs.m1")? (hint: look at AIC)

cjs.m3 <- crm(dipper.proc, 
              dipper.ddl,
              model.parameters = list(Phi = Phi.sex, 
                                      p = p.sex),
              accumulate = FALSE)
cjs.m3

## ----cjsfunction---------------------------
# You almost always will fit more than a few models, so manually creating each model
# would take a long time and have a high chance of making a mistake.
# It is more efficient (and tidy) to use a function to automate model creation 
# and put them in a table to rank based on AIC or QAIC
fit.dipper.cjs.models <- function(){
    
    # Apparent survival (Phi) formula
    Phi.sex.time <- list(formula=~sex*time)  # Just like in other linear models "*" includes main effects and an interaction
    Phi.time <- list(formula=~time) # differs between discrete times
    Phi.sex <- list(formula=~sex) # differs between males and females
    Phi.dot <- list(formula=~1) # constant survival
    
    # Detection probability (p) formula
    p.sex <- list(formula=~sex)  # differs between males and females
    p.time <- list(formula=~time)  # one discrete estimate of p per capture event
    p.dot <- list(formula=~1) # constant detection
    
    # Construct all combinations and put into one model table
    cml <- create.model.list(c("Phi","p")) # makes all possibile combinations of those parameter formulas
    results <- crm.wrapper(cml, 
                           data = dipper.proc, 
                           ddl = dipper.ddl,
                           external = FALSE, 
                           accumulate = FALSE, 
                           hessian = TRUE)
    return(results)
}

# Run function
dipper.cjs.models <- fit.dipper.cjs.models()

# Display model table
dipper.cjs.models

## ----floodcjs--------------------------
# Sometimes, you may want to model a variable that changes through time, 
# rather than a variable that differs between groups or individuals (e.g. sex, or some other demographic variable)

# As an example, we will model how flooding affects survival in dippers
# time varying variables can be added to the design matrix, which is part of the reason
# I suggested to get into the practice of creating processed and design data

# We will model Flood as a variable with two levels (floods during events 3, and 4)
# Note there is a difference between "time" and "Time" (very tricky!)

# Phi flood
dipper.ddl$Phi$Flood <- "noflood" # New variable set as 0
dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==3 | dipper.ddl$Phi$time==4] <- "flood" # Phi 1 is survival from time 1 to time 2

# p flood
dipper.ddl$p$Flood <- "noflood"
dipper.ddl$p$Flood[dipper.ddl$p$time==3 | dipper.ddl$p$time==4] <- "flood"

# Check it worked
head(dipper.ddl$Phi)
head(dipper.ddl$p)

# Make a model with Flood variable
Phi.Flood <- list(formula =~Flood)
p.Flood <- list(formula =~Flood)

dipper.cjs.flood <- crm(dipper.proc, 
                        dipper.ddl,
                        model.parameters = list(Phi = Phi.Flood, p = p.Flood),
                        accumulate = FALSE, hessian = TRUE)

predict(dipper.cjs.flood, se = TRUE)
# OK, our (fake) flooding value doesn't seem to strongly affect survival or detection

