list.of.packages <- c("ggplot2", "car", "nlme", "dplyr", "reshape2", "plotly", "tidyr", "ggthemes", "lme4", "multcomp", "rcompanion", "gridExtra", "cowplot", "survival", "here", "ggpubr", "rstatix", "WRS2", "tidyverse", "readxl", "purrr", "emmeans") #add new libraries here 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load all libraries 
lapply(list.of.packages, FUN = function(X) {
  do.call("require", list(X)) 
})

sessionInfo()

