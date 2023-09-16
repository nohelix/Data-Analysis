
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# where to save graphs
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Individual graphs"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41


# Load Packages -----------------------------------------------------------

# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(broom)

# analysis & visualization
library(psycho); library(ggplot2); library(nortest); library(hrbrthemes); library(gtools)



# Load Necessary Datasets -------------------------------------------------
load(paste0(projects_folder, "run_archive.Rdata"), .GlobalEnv)
rat_decoder = data.table::fread("rat_decoder.csv")

# Remove bad data
#TODO: deal with multiple runs in a day
dataset = run_archive %>% 
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe 
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE))






