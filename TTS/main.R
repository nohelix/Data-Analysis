# NOTES -------------------------------------------------------------------
# Major re-write on 2/25/2023 for new system of trials


# Load Packages -----------------------------------------------------------
# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(glue); library(lubridate); library(broom)

# analysis & visualization
library(psycho); library(ggplot2); library(nortest); library(hrbrthemes); library(gtools); 
library(FSA); library(nlstools); library(nlraa)
# FSA provides an SE calc


# # Clear workspace ---------------------------------------------------------
# rm(list = ls())


# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# Location of the R scripts
code_folder = "Y:/GitHub/Data Analysis/TTS"

# Explict location to save files to:
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/TTS Graphs"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41

# Working directory -------------------------------------------------------
setwd(code_folder)


# Import Data -----------------------------------------
# errors, post ABRs and 'maintenance' days are automatically removed
source(glue("{code_folder}/TTS_data.R"))


# # Analysis ----------------------------------------------------------------
# # Graph Hit, FA and Trial Count from summary data
# # calculate hearing threshold (from all data) and remove any trials below hearing level.
# # ISSUE: Plots don't show.
#
# source(glue("{code_folder}/TTS_analysis.R"))
#
# # Graphing ----------------------------------------------------------------
# # ISSUE: Plots don't show.
#
# source(glue("{code_folder}/TTS_graphs.R"))

