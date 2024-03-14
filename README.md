# Workshop: Mark-Recapture Analyses in R

A workshop on mark-recapture analyses in R. This workshop uses the 'marked' package to construct open-population mark-recapture models in R. All examples use the 'dipper' dataset contained within the 'marked' package.

This workshop covers:

* Formatting data into capture histories
* Fitting Cormack-Jolly-Seber models to estimate survival and detection probability
* Testing the goodness-of-fit of Cormack-Jolly-Seber models
* Fitting Jolly-Seber models to estimate population size
* Fitting known-fate mark-recapture models

Contents:
* 01_installingmarkrecappackages.R checks if required packages are installled. Installs them if they are not.
* 02_markrecaptureworkshopcode.R contains all code used in the workshop
* 03_introduction_to_CJS.R contains code for introduction to Cormack-Jolly-Seber Models
* 04_gof_for_cjs.R contains code for goodness-of-fit tests for Cormack-Jolly Seber Models
* 05_chat_for_cjs.R contains code for estimating variance inflation factor and adjusting model output for Cormack-Jolly-Seber Models
* 06_js_popan.R contains code for estimating population size using Jolly-Seber Models (POPAN formulation)
* 07_knownfate.R contains code for fitting known-fate mark-recapture models with RMark
