
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# FA cutoff
FA_cutoff = .40


# Load Packages -----------------------------------------------------------

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); library(purrr); library(forcats)

# analysis & visualization
library(ggplot2); library(nortest)



# Load Necessary Datasets -------------------------------------------------
load(paste0(projects_folder, "run_archive.Rdata"), .GlobalEnv)

# # Individual Trial Data
# load(paste0(projects_folder, "Oddball-LE_archive.Rdata"), .GlobalEnv)

# Remove bad data
#TODO: deal with multiple runs in a day
dataset = run_archive %>% 
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe 
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE))


# Get core data -----------------------------------------------------------

core_columns = c("date", "rat_name", "rat_ID", 
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", "block_size", "complete_block_count", 
                 "dprime", "reaction", "FA_percent")

core_data = dataset %>% 
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% unnest_wider(stats) %>%
  select(all_of(core_columns)) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("Oddball")) %>%
  mutate(task = str_replace_all(task, "CNO \\(3mg/kg\\)", "CNO 3mg/kg"))


# Get Reaction times by rat -----------------------------------------------

Rxn_table = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unquie
  group_by(rat_ID, task, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>% 
  unique() %>%
  rename(Position = `Inten (dB)`)


# Rxn analysis ------------------------------------------------------------

# Normality testing
ad.test(Rxn_table$Rxn)

Model_data = Rxn_table %>% 
  filter(task %in% c("Rotating", "CNO 3mg/kg")) %>%
  filter(rat_ID != 100)

## Overall Model
Rxn_overall_model = aov(Rxn ~ task * Position, 
                data = Model_data)

summary(Rxn_overall_model)

kruskal.test(Rxn ~ task, 
             data = Model_data)

Rxn_overall_model_postHoc <- 
  FSA::dunnTest(Rxn ~ interaction(task), 
                data = Rxn_table,
                method = "bonf")

Rxn_overall_model_postHoc$res %>% 
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         P.adj = round(P.adj, digits = 3))



