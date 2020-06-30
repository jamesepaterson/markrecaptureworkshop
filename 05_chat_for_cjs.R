## ---------------------------
##
## Script name: 05_chat_for_cjs
##
## Purpose of script: Calculate variance inflation factors and adjust model output of CJS
##
## Author: Dr. James Paterson
##
## Date Created: 2020-06-30
##
## 
## Email: james.earle.paterson@gmail.com
##
## ---------------------------


# chat for mark-recapture models
# Author: James E Paterson 
# Contact: james.earle.paterson@gmail.com
# Date created: 2020-06-30


## ----loadlibraries---------------------
# install.packages("R2ucare") # first time only
library(R2ucare) # For Goodness of fit tests
library(dplyr) # for tidy data
library(magrittr) # for pipes

## ----reshapedata---------------------------------------------------------
# Load dipper data (in "marked", but also in R2ucare package in different format)
# I'm loading `marked` version to be consistent with other tutorials/workshops
# and for fitting models below in "marked".
data(dipper, package = "marked")

# Full data
dipper.ch.gof <- dipper$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper))

# Females only
dipper.fem.ch.gof <- dipper$ch[dipper$sex == "Female"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Female",]))

# Males only
dipper.mal.ch.gof <- dipper$ch[dipper$sex == "Male"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Male",]))

## ----overallchatestimate-------------------------------------------------
library(R2ucare)

# Male overall
overall_CJS(dipper.mal.ch.gof, rep(1,nrow(dipper[dipper$sex == "Male",])))

# Female overall
overall_CJS(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))

# Overall chat for dipper (male test stat + female test stat) / (male df + female df)
dipper.chat <- (10.276 + 11.062) / (12 + 9)

dipper.chat # 1.016095, still very close to 1 and indicates relatively good fit

## ----createmodeltable----
detach("package:R2ucare", unload=TRUE) # Have to be careful with mutliple packages with "dipper" data set
# Here I detach R2ucare to avoid confusion

library(marked)
# Process data (and set grouping variables)
dipper.proc <- process.data(dipper, 
                            group = "sex") # group variable has to be in the data

# Make design data (from processed data)
dipper.ddl <- make.design.data(dipper.proc)

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

## ----displaycjstabls-----------------------------------------------------
# Display model table
dipper.cjs.models

## ----adjustmodelselection------------------------------------------------
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

dipper.cjs.models # display models
