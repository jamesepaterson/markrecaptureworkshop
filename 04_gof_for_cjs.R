# Goodness-of-fit tests for mark-recapture models
# Author: James E Paterson 
# Contact: james.earle.paterson@gmail.com
# Date created: 2020-05-20


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

## ----overallcjsgoffem----------------------------------------------------
# first argument = capture history matrix, second argument = capture history frequency (vector of 1's for our example)
overall_CJS(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",]))) # Females only

## ----test2ct-------------------------------------------------------------
# first argument = capture history matrix, second argument is frequency of each capture history (1 for example)
test2ct_fem <- test2ct(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",]))) 
test2ct_fem

## ----test2cl-------------------------------------------------------------
test2cl_fem <- test2cl(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test2cl_fem

## ----test3sr-------------------------------------------------------------
test3sr_fem <- test3sr(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test3sr_fem

## ----test3sm-------------------------------------------------------------
test3sm_fem <- test3sm(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test3sm_fem

## ----test2test3addition--------------------------------------------------
4.985+2.041+3.250+0

