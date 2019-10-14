### Mark-recapture models in R
# Author: James E Paterson james.earle.paterson@gmail.com
# Date created: 2019-09-03

# Outline
# 1. load packages
# 2. Create capture histories from "long-format" capture data
# 3. Run basic Cormack-Jolly-Seber model on dipper
# 4. More flexible and advanced CJS models
# 5. Goodness-of-fit Tests for CJS models
# 6. Jolly-Seber (POPAN) models to estimate abundance

##### 1. Load necessary packages ----------------------------------

# Load packages
library(marked)  # For creating models
library(ggplot2) # For plotting
library(reshape2) # for formatting data

##### 2. Create capture histories ------------------------------------------------

# Simulate some data (long-form)
set.seed(1111)
input.data <- data.frame(id = c(sample(1:50, 20, replace = FALSE), 
                                sample(1:50, 20, replace = FALSE), 
                                sample(1:50, 20, replace = FALSE)),
                         event = c(rep(1, 20), rep(2, 20), rep(3, 20)),
                         detect = 1)

# Convert into wide format capture histories
junk <- reshape2::melt(input.data, 
                       id.var = c("id","event"), # id = animal identifier, event = time indicator
                       measure.var="detect") # all detect = 1

y <- reshape::cast(junk, id ~ event) # Change shape to have a row for every individual, column for every event
y[is.na(y)] = 0 # Change NA's into 0's
k <- dim(y)[2] # Number of animals (rows)

# needed to deal with repeats on a capture event. May or may not be an issue in your data
y[,2:k] <- (y[,2:k]>0)*1

# Function to collapse capture histories into one vector (combine columns of events)
pasty<-function(x) 
{
  k<-ncol(x) # number of events
  n<-nrow(x) # number of individuals
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,]>0)*1
    out[i]<-paste(y, sep = "", collapse = "")
  }
  return(out)
}

capt.hist <- data.frame(ch = pasty(y[,-1]), # ch = capture history. "-1" means exclude first column, which has id
                        id = y[,1])  # id column

# Save for another analysis
save(capt.hist, file = "createdcapthist.R")

load(file = "createdcapthist.R")
# If you have a sample data set, 
# try to create capture histories in this space


##### 3. Run basic Cormack-Jolly-Seber model on dipper ------------------------

# Cormack-Jolly-Seber models estimate:
# Apparent survival (Phi), and
# Detection probability (p) for open populations

# Load example dipper dataset
data(dipper)

# Loak at data
head(dipper)

# Create a basic CJS model (constant survival and detection) with: 
cjs.m1 <- crm(dipper)

# Examine model and coefficient estimates
cjs.m1

# Re-fit model with precision estimates
cjs.m1 <- cjs.hessian(cjs.m1)
cjs.m1

# The estimates are on a log-linear scale, let's change them back.
# The beta for Phi and p was estimated with a logit link
# The inverse logit will convert beta for survival back to the real survival estimate
# exp(x)/(1+exp(x)) or R function "plogis"
exp(cjs.m1$results$beta$Phi)/(1+exp(cjs.m1$results$beta$Phi))
plogis(cjs.m1$results$beta$Phi)
plogis(cjs.m1$results$beta$p)

# Alternatively, we can use the "predict" function
predict(cjs.m1, 
        newdata = data.frame(sex = c("Female", "Male")),
        se=TRUE)  # In this case there are no groups or covariates, so "new data" not used

# IMPORTANT NOTE: we assumed all time intervals between captures = 1, 
# you can change this with a vector of time intervals
# Often, time intervals are not equal (e.g. due to weather constraints or missing surveys)
cjs.m1.unequaltime <- crm(dipper, time.intervals = c(1,3,1,1,4,5))
predict(cjs.m1.unequaltime)
# Notice the predicted survival is higher now because the time length between capture events is always >= 1

##### 4. More flexible and advanced CJS models ------------------------------

# It is good practice to follow these steps:
# 1. Process data
# 2. Create design data
# 3. Set-up and execute entire candidate model set
# The advantage is it is easier to add covariates (e.g. time-varying effects from weather or experiments)

# Process data
dipper.proc <- process.data(dipper, 
                            group = "sex")

# Make design data (from processed data)
dipper.ddl <- make.design.data(dipper.proc)

# Look at design data
dipper.ddl

# Outine formulas for each parameter
Phi.dot <- list(formula=~1)  # ~1 is always a constant (or single estimate)
Phi.sex <- list(formula=~sex)
p.sex <- list(formula=~sex) # Be careful of case-sensitive names. Use the exact group column that was in data

# Make new model (using design data) with constant survival, but different detection probabilities between sexes
cjs.m2 <- crm(dipper.proc, 
              dipper.ddl,
              model.parameters = list(Phi = Phi.dot, p = p.sex),
              accumulate = FALSE)
cjs.m2
# Is this a better fit than Phi.dot.p.dot (i.e. "model.1")? (hint: look at AIC)

cjs.m3 <- crm(dipper.proc, 
              dipper.ddl,
              model.parameters = list(Phi = Phi.sex, p = p.sex),
              accumulate = FALSE)
cjs.m3

# You almost always will fit more than a few models, so this approach would get lengthy.
# It is more efficient (and tidy) to use a function to automate model creation 
# and put them in a table to rank based on AIC or QAIC
fit.dipper.cjs.models <- function(){
    
    # Apparent survival (Phi)
    Phi.sex.time <- list(formula=~sex*time)  # Just like in other linear models "*" includes main effects and an interaction
    Phi.time <- list(formula=~time) # differs between discrete times
    Phi.sex <- list(formula=~sex) # differs between males and females
    Phi.dot <- list(formula=~1) # constant survival
    
    # Detection probability (p)
    # p.sex.time <- list(formula=~sex*time)
    p.sex <- list(formula=~sex)  # differs between males and females
    p.time <- list(formula=~time)  # one discrete estimate of p per capture event
    p.dot <- list(formula=~1) # constant detection
    
    # Construct all combinations and put into one model table
    cml <- create.model.list(c("Phi","p")) # makes all possibile combinations of those parameter formulas
    results <- crm.wrapper(cml, data = dipper.proc, ddl = dipper.ddl,
                          external = FALSE, accumulate = FALSE, hessian = TRUE)
    return(results)
}

# Run function
dipper.cjs.models <- fit.dipper.cjs.models()

# Display model table
dipper.cjs.models

# Look at estimates of top model (number on left, or using name)
dipper.cjs.models[[1]] # or dipper.cjs.models[["Phi.dot.p.dot"]] or dipper.cjs.models$Phi.dot.p.dot

# Getting "real" estimates using "predict"
# First, create new data to predict differences
newdipper <- data.frame(sex = c("Female", "Male"))

reals <- predict(dipper.cjs.models[[2]], # NOTE, i'm deliberately calling a model that isn't the "top-performing
                 # model" in order to demonstate the similar detection probabilities between males and females
                 newdata = newdipper, 
                 se = TRUE)
reals

# Plot
ggplot(reals$p, aes(x = sex, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ucl, ymax = lcl), width = 0.1) +
  labs(y = "Detection probability") +
  coord_cartesian(ylim = c(0,1))

# Dealing with identifiability
dipper.cjs.models[[12]] # Look at fully time-dependent model (model 12)
# NOTICE the huge se for the last parameters
# In this model, we actually can't disentangle the last Phi or p (only the product of both of them)
# Be careful when using fully time-dependent models because not all parameters are identifiable


# Sometimes, you may want to model a variable that changes through time, 
# rather than a variable that differs between individuals (e.g. sex, or some other demographic variable)

# As an example, we will model how flooding affects survival in dippers
# time varying variables can be added to the design matrix, which is part of the reason
# I suggested to get into the practice of creating processed and design data

# We will model Flood as a variable with two levels (floods during events 3, and 4)
# Note there is a difference between "time" and "Time" (very tricky!)

# Phi flood
dipper.ddl$Phi$Flood <- 0 # New variable set as 0
dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==3 | dipper.ddl$Phi$time==4] <- 1 # Phi 1 is survival from time 1 to time 2

# p flood
dipper.ddl$p$Flood <- 0
dipper.ddl$p$Flood[dipper.ddl$p$time==3 | dipper.ddl$p$time==4] <- 1

# Check it worked
dipper.ddl$Phi
dipper.ddl$p

# Make a model with Flood variable
Phi.Flood <- list(formula =~Flood)
p.Flood <- list(formula =~Flood)

dipper.cjs.flood <- crm(dipper.proc, 
                        dipper.ddl,
                        model.parameters = list(Phi = Phi.Flood, p = p.Flood),
                        accumulate = FALSE, hessian = TRUE)

predict(dipper.cjs.flood, se = TRUE)
# OK, our (fake) flooding value doesn't seem to strongly affect survival or detection

# If you have a variable that changes through time (or between groups), 
# try adding it to your design data here.



# But, does our general model adequately fit the data? 
# We need to do some goodness-of-fit tests

##### 5. Goodness-of-fit Tests for CJS models -------------------------------

library(R2ucare) # For Goodness of fit tests
library(dplyr)
library(magrittr)

# Cormack-Jolly-Seber model assumptions:
# 1. Every marked animal present in the population at time [t] has the same probability of recapture (p[t])
# 2. Every marked animal in the population immediately after time[t] has the same probability of surviving to time [t+1] 
# 3. Marks are not lost or missed.
# 4. All samples are instantaneous, relative to the interval between occasion [t] and [t+1], and each release is made immediately after the sample.

# Assumptions 3 and 4 aren't tested, but we can use U-CARE (via R2ucare) to test assumptions 1 and 2.
# CARE = CApture REcapture
# R2ucare is a translated version of U-CARE, and doesn't require any installation of U-CARE

# Test 2: loosely tests equal catchability. Are all animals (known to be alive) equally detected?
# Test 3: Loosely tests differences in survival. Do all animals (known to be alive) have equal survival?

# R2ucare requires data in a slightly different format (a matrix with a column for each capture event)
dipper.ch.gof <- dipper$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper))

# Females
dipper.fem.ch.gof <- dipper$ch[dipper$sex == "Female"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Female",]))

# Males
dipper.mal.ch.gof <- dipper$ch[dipper$sex == "Male"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Male",]))
  
# Overall GOF for all individuals (not separated by sex)
# In general, combining Test 2 and Test 3 can be dangerous because you would not see 
# the cause of any heterogeneity or asssumption violation
overall_CJS(dipper.ch.gof, rep(1,nrow(dipper))) # All pooled together (not considering sex differences). 
# first argument = capture history matrix, second argument = capture history frequency (vector of 1's for our example)
# These tests are always for "time-dependent" models, i.e. "Phi.time.p.time"

# Better to separate components of Test 2 and Test 3
# Below we perform tests on females, then you can try on the male subset

# Test 2 (females); we need the m-array to perform test2ct and test2cl. Does recapture depent on when they were first marked?
# Test 2CT: tests whether there is a difference in p at [t+1] between those captured and not captured at [t] (when known to be alive)
test2ct_fem <- test2ct(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))  # first argument = capture history matrix, second argument is frequency of each capture history (1 for example)
# Test2CL: is there a difference in the expected time of next recapture between individuals captured and not captured at [t] (when known to be alive)
test2cl_fem <- test2cl(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
# Notice that some tests are "none" because of low sample size; they are pooled with other capture events

# Test 3 (females)
# Test3SR: Does marking affect survival? Do individuals with previous marks have different survival than first-time captures?
test3sr_fem <- test3sr(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
# Test3sm: For animals seen again, does when they are recaptured depend on when they were marked on or before time[i]
test3sm_fem <- test3sm(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))

# The "overall" test statistic is the sum of all component tests (below for female tests)
4.985+2.041+3.250+0 # 10.276

# The same as Test1 (overall or omnibus test):
overall_CJS(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",]))) # All pooled together

# Exercise: apply the same component tests to males




# What is "chat" or the variance inflation factor (or overdispersion estimate)?
# Can help adjust model output based on extra-binomial noise
# It is also an estimate of the fit of the model to the observed data
# One estimate is overall X2/df (sum from Test2 and Test3 components)
# chat < 1: generally ignore slight deviations,
# chat = 1, perfect fit
# chat > 1, extra noise, can adjust model tables to use QAIC. chat > 3 generally represents poor fit 
# (maybe change the general model structure)
# This methods of chat estimation only applies if you're using a time-dependent model as the general model
# For other approaches, look up the median chat or bootstrap methods in Gentle Introduction to MARK

# Estimate chat for females
10.276/12  # 0.856333 for females

# Male overall
overall_CJS(dipper.mal.ch.gof, rep(1,nrow(dipper[dipper$sex == "Male",]))) # All pooled together

# Overall chat for dipper (male test stat + female test stat) / (male df + female df)
# This represents a test for the model: Phi.sex*time.p.sex*time
dipper.chat <- (10.276 + 11.062) / (12 + 9)
# 1.016095, still very close to 1 and indicates relatively good fit

detach("package:R2ucare", unload=TRUE) # Have to be careful with mutliple packages with "dipper" data set
# Here I detach R2ucare to avoid confusion

# Adjust our CJS model table to account for chat different than 1 (chat = 1.016095)
# QAIC = (-2lnLik / chat) + 2K
# This favours simpler models as chat gets larger than 1
dipper.cjs.models$model.table$chat <- dipper.chat
dipper.cjs.models$model.table$QAIC <- (dipper.cjs.models$model.table$neg2lnl/dipper.chat) + (2*dipper.cjs.models$model.table$npar)

# Note: this doesn't autmatically update the order of the table! 
# In our case, the order of supported models did not change
# To sort:
dipper.cjs.models$model.table <- dipper.cjs.models$model.table %>%
  arrange(QAIC)

dipper.cjs.models

# The order of models will favour simpler models as chat increases greater than one
# But in our case, the models are largely ordered by most simple to more complicated.

##### 6. Jolly-Seber models ---------------------------------------------

# Jolly-Seber models (POPAN formulation) are open population models, and 
# can be used to estimate abundance by including two more parameters than the CJS

# Additional parameters:
# Nhat (or "superpopulation") = total number of individuals available to enter population throughout study
# pent ("probability of entry") =  the rate at which individuals enter the population from Nhat (via births and immigration)

# Number of new animals in at time [i] = Nhat * pent[i]
# N at an time[i] = N[i-1] * phi[i] + (Nhat * pent[i])

# WARNING: there is no adequate GOF tests for Jolly-Seber models. 
# One common method: Test equivalent structure of CJS model with R2ucare (see above). 
# This tests *some* assumptions of Phi and p.
# Jolly-Seber models have an additional assumption:
# marked AND unmarked animals have same p (R2ucare doesn't test this)
# This assumption is required to estimate total abundance (sum of marked and unmarked animals in population)

# First, process data (Notice model = "js", previous version = "cjs")
dipper.js.proc <- process.data(dipper, 
                               model = "JS", 
                               groups = "sex")

# Second, make design data (from processed data)
dipper.js.ddl <- make.design.data(dipper.js.proc)

fit.js.dipper.models <- function(){
  # Phi formulas
  Phi.dot <- list(formula=~1)
  # p formulas
  p.dot <- list(formula=~1)
  # pent formulas. pent estimates MUST SUM to 1 (for each group).
  # This is constained using a Multinomial Logit link
  pent.time <- list(formula=~time)
  pent.sex <- list(formula=~sex)
  pent.dot <- list(formula=~1)
  # Nhat formulas
  N.sex <- list(formula=~sex)
  N.dot <- list(formula=~1)
  cml <- create.model.list(c("Phi","p", "pent", "N"))
  results <- crm.wrapper(cml, data = dipper.js.proc, ddl = dipper.js.ddl,
                         external = FALSE, accumulate = FALSE, hessian = TRUE)
  
  return(results)
}

# Run function
dipper.js.models <- fit.js.dipper.models()

# Display model table
dipper.js.models

# Display top model
dipper.js.models[[1]]

# Display model with pent.time
dipper.js.models[[5]]

# Predict (real) values using top model
dipper.js.predicted <- predict(dipper.js.models[[1]], newdipper, se = TRUE)

# To examine how pent sums to 1 for each group, look at model 3 where pent(~sex)
dipper.js.predicted.3 <- predict(dipper.js.models[[3]], newdipper, se = TRUE)

# Predict with pent.time (model 5)
dipper.js.predicted.5 <- predict(dipper.js.models[[5]], newdipper, se = TRUE)
dipper.js.predicted.5

# sex difference in nsuper
dipper.js.predicted.2 <-  predict(dipper.js.models[[2]], newdipper, se = TRUE)

# How do we get population size at each capture event?
# Abundance (N) is derived from the estimated parameters
# We will estimate population size at each time by making a dataframe of estimates and calculating N
# We will use the predicted estimates from the top-performing model (in this case: "dipper.js.predicted")

# NOTE: the below method will have to be adjusted based on your final model and the number of capture events
N.derived <- data.frame(occ = c(1:7),
                        Phi = c(rep(dipper.js.predicted$Phi$estimate, 6), NA),  
                        Nhat = rep(dipper.js.predicted$N$estimate + nrow(dipper), 7),
                        pent = c(1-sum(dipper.js.predicted$pent$estimate), dipper.js.predicted$pent$estimate))

# Set-up empty vectors for calculating N and confidence intervals
N.derived$N <- NA
N.derived$N.lcl <- NA
N.derived$N.ucl <- NA

# The inital population size (N[1]) = Nhat * (1 - sum(all other pent estimates))
# This is because of the link function for estimating pent.
# The sum of all pent parameters MUST equal 1 (therefore, one less must be estimated)
N.derived$N[1] <- (N.derived$Nhat[1] * N.derived$pent[1])
# Adding SE / confidence intervals around population size estimates
N.derived$N.ucl[1] <- (dipper.js.predicted$N$ucl + nrow(dipper))*(1-sum(dipper.js.predicted$pent$lcl))
N.derived$N.lcl[1] <- (dipper.js.predicted$N$lcl + nrow(dipper))*(1-sum(dipper.js.predicted$pent$ucl))

# Subsequent population sizes are estimated by calculating surviving individuals (N[t-1] * Phi[t]), and
# Adding new births (Nhat * pent[t])
for(i in 2:nrow(N.derived)){
  N.derived$N[i] <- (N.derived$N[i-1]*N.derived$Phi[i-1]) + (N.derived$Nhat[i] * N.derived$pent[i])
  N.derived$N.ucl[i] <- (N.derived$N.ucl[i-1]*dipper.js.predicted$Phi$ucl[1]) + (((dipper.js.predicted$N$ucl + nrow(dipper))) *dipper.js.predicted$pent$ucl[i-1])
  N.derived$N.lcl[i] <- (N.derived$N.lcl[i-1]*dipper.js.predicted$Phi$lcl[1]) + (((dipper.js.predicted$N$lcl + nrow(dipper))) *dipper.js.predicted$pent$lcl[i-1])
}

# Look at what we did
N.derived

# Plot population size and confidence interval
ggplot(N.derived, aes(x = occ, y = N)) +
  geom_point(shape = 16, size = 3, col = "black") +
  geom_errorbar(aes(ymin = N.lcl, ymax = N.ucl), width = 0.3) +
  geom_path() +
  labs(x = "Capture event", y = "Estimated population size\n(and 95% confidence interval)") +
  theme_classic()

# If you have your own data, try to estimate abundance using a Jolly-Seber POPAN model below
# I recommend initially starting with simple model structure



