# Get core data -----------------------------------------------------------
source("Fmr1 vs Tsc2 in LE data.R")


core_columns = c("date", "rat_name", "rat_ID", 
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", "block_size", "complete_block_count", 
                 "dprime", "reaction", "FA_percent")

core_data = dataset %>% 
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% unnest_wider(stats) %>%
  select(all_of(core_columns)) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("Fmr1-LE", "Tsc2-LE")) %>%
  # drop Oddball and Octave
  filter(! phase %in% c("Octave", "Tone-BBN"))

# Calculate Overall TH ----------------------------------------------------

# Threshold calculation calculation based on TH_cutoff intercept of fit curve
# LOESS: Local Regression is a non-parametric approach that fits multiple regressions
# see http://r-statistics.co/Loess-Regression-With-R.html
Calculate_TH <- function(df) {
  # Uncomment to see line fitting by a package which shows line
  # library(drda)
  # drda(dprime ~ dB, data = df) %>% plot
  fit = loess(dprime ~ dB, data = df)
  # plot(fit)
  TH = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")$y
  # print(TH)
  return(TH)
}

TH_table = core_data %>%
  # Omit Training & Reset days
  filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get dprimes
  unnest(dprime) %>%
  # Sort for ordered
  arrange(rat_ID, rat_name, Freq, Dur, dB) %>%
  #Prep for Calculate_TH function
  nest(data = c(dB, dprime), .by = c(rat_ID, rat_name, detail, Freq, Dur)) %>% 
  mutate(TH = map_dbl(data, Calculate_TH)) %>%
  select(-data) 


# Get Reaction times by rat -----------------------------------------------

Rxn_table = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unquie
  group_by(rat_ID, rat_name, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>% 
  unique()

# Limit to Over TH
TH_filter <- function(df) {
  ID = unique(df$ID)
  Dur = unique(df$Dur)
  kHz = unique(df$Freq)
  # kHz = if_else(kHz == "0", "BBN", paste0(kHz,"kHz")) # %>% print
  
  cuttoff = TH_table %>% # may have to use UQ to force the evaluation of the variable
    filter(Dur == Dur & rat_ID == ID & Freq == kHz) %>% 
    .$TH
  
  cuttoff = ifelse(identical(cuttoff, numeric(0)), -99, cuttoff)
  
  r = df %>%
    filter(Inten >= UQ(cuttoff))
  
  return(r)
}

Rxn_table_over_TH = Rxn_table %>% 
  # Prep for TH_filter function
  ungroup() %>%
  mutate(ID = rat_ID, Dur = `Dur (ms)`, Freq = `Freq (kHz)`, Inten = `Inten (dB)`) %>%
  nest(data = c(ID, Freq, Dur, Inten, Rxn), .by = c(rat_ID, rat_name, detail, `Freq (kHz)`, `Dur (ms)`)) %>%
  # Apply TH_filter
  mutate(data = map(data, TH_filter)) %>%
  unnest(data) %>%
  select(- any_of(c("Dur (ms)", "Freq (kHz)", "ID"))) %>%
  rename(Intensity = Inten, Duration = Dur, Frequency = Freq) %>%
  # Use rat_ID because its sure to be unquie
  group_by(rat_ID, rat_name) %>%
  # Get Averages
  mutate(Rxn = mean(Rxn, na.rm = TRUE)) %>% 
  unique()



# Time to learning --------------------------------------------------------

Learning_streak = core_data %>%
  dplyr::filter(task %in% c("Training")) %>%
  group_by(rat_ID) %>% # pregroup this so that numbering restarts at 1 for each rat
  mutate(groupid = data.table::rleid(phase, task, detail)) %>%
  group_by(rat_ID, groupid, phase, task, detail) %>%
  summarise(running_streak = n(), .groups = "drop")


# Decode for analysis -----------------------------------------------------
core_data = left_join(core_data,
                      select(rat_decoder, all_of(c("rat_ID", "line", "genotype", "sex"))),
                      by = "rat_ID")

TH_table = left_join(TH_table,
                     select(rat_decoder, all_of(c("rat_ID", "line", "genotype", "sex"))),
                     by = "rat_ID") %>%
  rename(Frequency = Freq, Duration = Dur)


Rxn_table_over_TH = left_join(Rxn_table_over_TH,
                              select(rat_decoder, all_of(c("rat_ID", "line", "genotype", "sex"))),
                              by = "rat_ID")

Rxn_table = left_join(Rxn_table,
                      select(rat_decoder, all_of(c("rat_ID", "line", "genotype", "sex"))),
                      by = "rat_ID") %>%
  rename(Frequency = `Freq (kHz)`, Duration = `Dur (ms)`, Intensity = `Inten (dB)`)

Learning_streak = left_join(Learning_streak,
                            select(rat_decoder, all_of(c("rat_ID", "line", "genotype", "sex"))),
                            by = "rat_ID")
