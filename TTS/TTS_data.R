# Load Necessary Datasets -------------------------------------------------
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)
rat_archive = fread(glue("{projects_folder}/rat_archive.csv"), select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))

# # Individual Trial Data
# load(paste0(projects_folder, "TTS_archive.Rdata"), .GlobalEnv)


# Get core data -----------------------------------------------------------

Get_Step_Size <- function(summary) {
  temp = summary$dB_step_size %>% unique

  if(length(temp) != 1) {
    step_size = Inf
    # print(temp)
  } else step_size = temp

  return(step_size)
}


core_columns = c("date", "rat_name", "rat_ID",
                 "file_name", "experiment", "phase", "task", "detail",
                 "stim_type", "analysis_type", "block_size", "step_size", "dB_min", "complete_block_count",
                 "dprime", "reaction", "FA_percent", "trial_count", "hit_percent")

#TODO: deal with multiple runs in a day, 
dataset = run_archive %>%
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE)) %>%
  # Omit Blue3 (267)
  filter(rat_ID != 267) %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("TTS")) %>%
  # Omit Training & Reset days
  filter(! phase %in% c("Octave")) %>%
  unnest_wider(stats) %>%
  rowwise %>%
  mutate(step_size = Get_Step_Size(summary),
         dB_min = summary$dB_min %>% min())

# TODO: handle multiple hearing losses
core_data = dataset %>%
  # Omit Training & Reset days
  filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # record date of hearing loss
  left_join(select(rat_archive, Rat_ID, HL_date), by = c("rat_ID" = "Rat_ID")) %>%
  # Calculate time since HL
  mutate(date = ymd(date), HL_date = ymd(HL_date),
         HL = case_when(is.na(HL_date) ~ -1,
                        date <= HL_date ~ -1,
                        date > HL_date ~ as.numeric(date - HL_date),
                        TRUE ~ -100),
         HL_state = case_when(HL == -1 ~ "baseline",
                              HL < 3 ~ "HL",
                              HL >= 3 && HL < 15  ~ "recovery",
                              HL >= 15 ~ "post-HL",
                              HL == -100 ~ "issue",
                              TRUE ~ "issue") %>% factor(levels = c("baseline", "HL", "recovery", "post-HL"), ordered = TRUE),
         BG = if_else(str_detect(file_name, pattern = "_BG_"), str_extract(file_name, pattern = "BG_.*$"), "None"),
         BG_type = case_when(str_detect(BG, pattern = "(PKN|PNK)") ~ "Pink",
                             str_detect(BG, pattern = "WN") ~ "White",
                             BG == "None" ~ "None",
                             TRUE ~ "ISSUE"),
         BG_Intensity = case_when(BG != "None" ~ str_extract(BG, pattern = "[:digit:]+"),
                                  BG == "None" ~ "None",
                                  TRUE ~ "ISSUE")) %>%
  select(-BG)


# Filter to rats who have post-HL measures
# rats_survived_to_post_HL = filter(core_data, HL_state == "post-HL") %>% .$rat_ID %>% unique %>% as.list()
# 
# core_data = filter(core_data, rat_ID %in% rats_survived_to_post_HL)

# Descriptive Stats -------------------------------------------------------
stats_table =
  core_data %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, stim_type) %>%
  # Get Averages
  summarise(trial_count = mean(trial_count, na.rm = TRUE),
            hit_percent = mean(hit_percent, na.rm = TRUE),
            FA_percent = mean(FA_percent, na.rm = TRUE),
            .groups = "drop")

stats_table_by_BG =
  core_data %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, stim_type, BG_type, BG_Intensity) %>%
  # Get Averages
  summarise(trial_count = mean(trial_count, na.rm = TRUE),
            hit_percent = mean(hit_percent, na.rm = TRUE),
            FA_percent = mean(FA_percent, na.rm = TRUE),
            .groups = "drop")


# dprime ------------------------------------------------------------------
# needed for TH and dprime table so stop at in-between
dprime_table =
  core_data %>%
  # re-nest by Frequency
  unnest(dprime)

dprime_detail_table =
  dprime_table %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, detail, Freq, Dur, dB) %>%
  # only use BBN because it was the ony one with temporal integration testing
  filter(Freq == "0") %>%
  # Get Averages
  transmute(dprime = mean(dprime, na.rm = TRUE)) %>%
  unique()

dprime_summary_table =
  dprime_table %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, Freq, Dur, dB) %>%
  # Get Averages
  transmute(dprime = mean(dprime, na.rm = TRUE)) %>%
  unique()

dprime_5step_table =
  dprime_table %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, step_size, BG_type, BG_Intensity, Freq, Dur, dB) %>%
  # Get Averages
  transmute(dprime = mean(dprime, na.rm = TRUE)) %>%
  unique()


# Calculate Overall TH ----------------------------------------------------

# Threshold calculation calculation based on TH_cutoff intercept of fit curve
# LOESS: Local Regression is a non-parametric approach that fits multiple regressions
# see http://r-statistics.co/Loess-Regression-With-R.html
Calculate_TH <- function(df) {
  # Sort for ordered - this needs to be by dB not date
  df = arrange(df, dB)
  fit = loess(dprime ~ dB, data = df)
  TH = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")$y
  return(TH)
}

TH_table =
  dprime_table %>%
  nest(dprime = c(rat_name, date, Freq, Dur, dB, dprime), .by = c(rat_ID, rat_name, Freq, HL_state, BG_type, BG_Intensity, Dur)) %>%
  # calculate TH
  rowwise() %>%
  mutate(TH = Calculate_TH(dprime)) %>%
  select(-dprime)

TH_table_detail =
  core_data %>%
  # re-nest by Frequency
  unnest(dprime) %>%
  nest(dprime = c(rat_name, date, Freq, Dur, dB, dprime), .by = c(rat_ID, rat_name, Freq, HL_state, BG_type, BG_Intensity, Dur, detail)) %>%
  # calculate TH
  rowwise() %>%
  mutate(TH = Calculate_TH(dprime)) %>%
  select(-dprime)

TH_table_steps =
  core_data %>%
  # re-nest by Frequency
  unnest(dprime) %>%
  filter(detail == "Alone") %>%
  nest(dprime = c(rat_name, date, Freq, Dur, dB, dprime), .by = c(rat_ID, rat_name, Freq, HL_state, BG_type, BG_Intensity, Dur, step_size)) %>%
  # calculate TH
  rowwise() %>%
  mutate(TH = Calculate_TH(dprime),
         step_size = as.character(step_size)) %>%
  select(-dprime) %>%
  rename(Frequency = Freq, Duration = Dur) %>%
  bind_rows(
    filter(TH_table_detail, detail == "Alone") %>%
      mutate(step_size = "Mixed") %>%
      select(-detail)
  )


# Get Reaction times by rat -----------------------------------------------

Rxn_table =
  core_data %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>%
  unique()

Rxn_table_detail =
  core_data %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>%
  unique()

CNO_data = unique(filter(core_data, task == "CNO")$rat_ID)

Rxn_table_CNO =
  core_data %>%
  filter(rat_ID %in% CNO_data) %>%
  filter(date < "2023-05-11") %>%
  mutate(detail = if_else(task == "CNO", "CNO","Baseline")) %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>%
  unique()

Rxn_5step_table =
  core_data %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, detail, step_size, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>%
  unique()


# Rename for ease of use --------------------------------------------------

dprime_summary_table = dprime_summary_table %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                                       Duration = c(contains("Dur") | contains("ms")),
                                                       Intensity = c(contains("dB")))

dprime_detail_table = dprime_detail_table %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                                      Duration = c(contains("Dur") | contains("ms")),
                                                      Intensity = c(contains("dB")))

dprime_5step_table = dprime_5step_table %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                                   Duration = c(contains("Dur") | contains("ms")),
                                                   Intensity = c(contains("dB")))

Rxn_table = Rxn_table %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                 Duration = c(contains("Dur") | contains("ms")),
                                 Intensity = c(contains("dB")))

Rxn_5step_table = Rxn_5step_table %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                 Duration = c(contains("Dur") | contains("ms")),
                                 Intensity = c(contains("dB")))

Rxn_table_detail = Rxn_table_detail %>% rename(Frequency = c(contains("Freq") | contains("kHz")),
                                 Duration = c(contains("Dur") | contains("ms")),
                                 Intensity = c(contains("dB")))

TH_table = TH_table %>% rename(Frequency = contains("Freq"),
                               Duration = c(contains("Dur") | contains("ms")))

TH_table_detail = TH_table_detail %>% rename(Frequency = contains("Freq"),
                                      Duration = c(contains("Dur") | contains("ms")))


# Rxn over TH calc --------------------------------------------------------

# Limit to Over TH
TH_filter <- function(df) {
  ID = unique(df$rat_ID)
  Dur = unique(df$Duration)
  kHz = unique(df$Frequency)
  HL = unique(df$HL_state)
  BG = unique(df$BG_type)
  BG_dB = unique(df$BG_Intensity)

  cuttoff = TH_table %>%
    filter(Duration == Dur & rat_ID == ID & Frequency == kHz & HL_state == HL & BG_type == BG & BG_Intensity == BG_dB) %>%
    .$TH

  cuttoff = ifelse(identical(cuttoff, numeric(0)), -99, cuttoff)

  r = df %>%
    filter(Intensity >= UQ(cuttoff))  # have to use UQ to force the evaluation of the variable

  return(r)
}

Rxn_table_over_TH = Rxn_table %>%
  # Prep for TH_filter function
  ungroup() %>%
  nest(data = c(rat_ID, HL_state, BG_type, BG_Intensity, Frequency, Duration, Intensity, Rxn), .by = c(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, Frequency, Duration)) %>%
  # Apply TH_filter
  mutate(data = map(data, TH_filter)) %>%
  mutate(tmp = map(data, ~ select(., Intensity, Rxn))) %>%
  unnest(tmp) %>%
  select(-data) %>%
  unique()

Rxn_table_over_TH_detail = Rxn_table_detail %>%
  # Prep for TH_filter function
  ungroup() %>%
  nest(data = c(rat_ID, HL_state, BG_type, BG_Intensity, Frequency, Duration, Intensity, Rxn), .by = c(rat_ID, rat_name, HL_state, BG_type, BG_Intensity, detail, Frequency, Duration)) %>%
  # Apply TH_filter
  mutate(data = map(data, TH_filter)) %>%
  mutate(tmp = map(data, ~ select(., Intensity, Rxn))) %>%
  unnest(tmp) %>%
  select(-data) %>%
  unique()

Rxn_table_over_TH_step_size = Rxn_5step_table %>%
  # Prep for TH_filter function
  ungroup() %>%
  nest(data = c(rat_ID, HL_state, BG_type, BG_Intensity, Frequency, Duration, Intensity, Rxn), .by = c(rat_ID, rat_name, step_size, HL_state, BG_type, BG_Intensity, Frequency, Duration)) %>%
  # Apply TH_filter
  mutate(data = map(data, TH_filter)) %>%
  mutate(tmp = map(data, ~ select(., Intensity, Rxn))) %>%
  unnest(tmp) %>%
  select(-data) %>%
  unique()





