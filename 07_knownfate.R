## ---------------------------
##
## Script name: 07_knownfate
##
## Purpose of script: Known-fate mark-recapture models in RMark
##
## Author: James E Paterson
##
## Date Created: 2024-03-13
##
## ---------------------------

## ----loaddata-----------------------------------------------------------------------
library(dplyr) # for data organization
library(ggplot2) # for plotting
library(RMark) # for mark-recapture models

# Load known-fate mark-recapture data on American Black Ducks
data("Blackduck") # 48 observations of 5 variables

# Look at first few rows to see structure
head(Blackduck)

## ----processdesigndata-------------------------------------------------------------------------------

# Process data. model = "Known" for known-fate models.
# See function details for other arguments, such as unequal time intervals and age variables
Blackduck.process <- process.data(Blackduck,
                             model = "Known", 
                             groups = c("BirdAge"))

# Create design data
Blackduck.ddl <- make.design.data(Blackduck.process)

# Look at first part of design data for the only parameter "S"
head(Blackduck.ddl$S)

## ----fitknownfatemodels--------------------------------------------
Blackduck_models <- function() {
S.Wing_Len <- list(formula =~Wing_Len)
S.condix <- list(formula =~condix)

#create a list of all possible models using possibilities above.  
cml <- create.model.list("Known")
results <- mark.wrapper(cml, 
                        data = Blackduck.process, 
                        ddl = Blackduck.ddl, 
                        adjust = TRUE, output = FALSE)
# adjust = TRUE since some known S parameters could be at boundaries
return(results)
}

# Execute the function
Blackduck_results <- Blackduck_models()

## ----viewresults-------------------------------------------------------------------------------------
# View results
Blackduck_results

# Look at model summary for most-supported model
summary(Blackduck_results$S.condix)


## ---- predictplotsurvival----------------------------------------------------------------------------
predicted_S <- predict_real(model = Blackduck_results$S.condix,
                            parameter = "S",
                            df = Blackduck.ddl$S,
                            replicate = TRUE,
                            data = data.frame(condix = seq(from = min(Blackduck$condix),
                                                         to = max(Blackduck$condix),
                                                         by = 0.01),
                                              Wing_Len = mean(Blackduck$Wing_Len),
                                              BirdAge = "1",
                                              Weight = mean(Blackduck$Weight)),
                            se = TRUE)

ggplot(data = predicted_S, aes(x = condix, y = real)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "gray80") +
  geom_smooth(col = "black", formula = y~x, method = "loess") +
  labs(x = "Body Condition Index", y = "Predicted survival") +
  theme_classic()

# Last step: clean temporary and unused files associated with RMark and MARK in working directory
# These are written whenever fitting models in MARK via RMark.
# Includes .inp, .out, .res, and .vcv
RMark::cleanup(ask = FALSE)
