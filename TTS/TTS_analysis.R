
# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05

  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 17)))}

}

# Descriptive Stats -------------------------------------------------------
# Trial Count
trial_count.aov.data = stats_table
trial_count.aov.data$Gaus = LambertW::Gaussianize(trial_count.aov.data$trial_count)[, 1]
trial_count.aov = aov(trial_count ~ HL_state * stim_type, data = trial_count.aov.data)

Parametric_Check(trial_count.aov)
# Normal un-transformed
summary(trial_count.aov)

trial_count.aov.stats = tidy(trial_count.aov) %>% mutate(sig = stars.pval(p.value))
# No post-Hoc testing as only stim_type is significant

# Trial Count for tones by BG
trial_count_by_BG.aov.data = stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL")),
         BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB")))
trial_count_by_BG.aov.data$Gaus = LambertW::Gaussianize(trial_count_by_BG.aov.data$trial_count)[, 1]
trial_count_by_BG.aov = aov(trial_count ~ HL_state * BG, data = trial_count_by_BG.aov.data)

Parametric_Check(trial_count_by_BG.aov)
# Normal un-transformed
summary(trial_count_by_BG.aov)

trial_count_by_BG.aov.stats = tidy(trial_count_by_BG.aov) %>% mutate(sig = stars.pval(p.value))
# No post-Hoc testing as nothing significant



# Hit% Count
hit_percent.aov.data = stats_table
hit_percent.aov.data$Gaus = LambertW::Gaussianize(hit_percent.aov.data$hit_percent)[, 1]
hit_percent.aov = aov(Gaus ~ HL_state * stim_type, data = hit_percent.aov.data)

Parametric_Check(hit_percent.aov)
# Normal with Gausian transform
car::qqPlot(hit_percent.aov.data$Gaus)

summary(hit_percent.aov)
hit_percent.aov.stats = tidy(hit_percent.aov) %>% mutate(sig = stars.pval(p.value))

# Post-Hoc
hit_percent.aov.postHoc =
  tidy(TukeyHSD(hit_percent.aov)) %>% mutate(sig = stars.pval(adj.p.value))

filter(hit_percent.aov.postHoc, sig != " ")

# Hit rate for tones by BG
hit_percent_by_BG.aov.data = stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB")))
hit_percent_by_BG.aov.data$Gaus = LambertW::Gaussianize(hit_percent_by_BG.aov.data$hit_percent)[, 1]
hit_percent_by_BG.aov = aov(Gaus ~ HL_state * BG, data = hit_percent_by_BG.aov.data)

Parametric_Check(hit_percent_by_BG.aov)
# Not normal even with transform but soooo close
car::qqPlot(hit_percent_by_BG.aov.data$hit_percent)
summary(hit_percent_by_BG.aov)

hit_percent_by_BG.aov.stats = tidy(hit_percent_by_BG.aov) %>% mutate(sig = stars.pval(p.value))
# Non-Parametric Post Hoc:
hit_percent_by_BG.aov.postHoc =
  hit_percent_by_BG.aov.data %>%
    group_by(BG) %>%
    summarise(dunns = tryCatch(FSA::dunnTest(hit_percent ~ HL_state, method = "none")$res, error = function(err){NA}),
              .groups = "drop") %>%
    unnest_wider(dunns) %>%
    filter(! is.na(Comparison)) %>%
    mutate(HL_state = NA) %>%
    relocate(HL_state, .after = BG) %>%
    bind_rows(hit_percent_by_BG.aov.data %>%
                group_by(HL_state) %>%
                summarise(dunns = tryCatch(FSA::dunnTest(hit_percent ~ BG, method = "none")$res, error = function(err){NA}),
                          .groups = "drop") %>%
                unnest_wider(dunns) %>%
                filter(! is.na(Comparison))) %>%
    mutate(P.adj = p.adjust(P.unadj, "BH"),
           sig = stars.pval(P.adj))

filter(hit_percent_by_BG.aov.postHoc, ! sig %in% c(" "))


# FA% Count
FA_percent.aov.data = stats_table
FA_percent.aov.data$Gaus = LambertW::Gaussianize(FA_percent.aov.data$FA_percent)[, 1]
FA_percent.aov = aov(FA_percent ~ HL_state * stim_type, data = FA_percent.aov.data)

Parametric_Check(FA_percent.aov)
# Normal
summary(FA_percent.aov)
FA_percent.aov.stats = tidy(FA_percent.aov) %>% mutate(sig = stars.pval(p.value))
# No post-Hoc testing as only stim_type is significant

# FA rate for tones by BG
FA_percent_by_BG.aov.data = stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB")))
FA_percent_by_BG.aov.data$Gaus = LambertW::Gaussianize(FA_percent_by_BG.aov.data$FA_percent)[, 1]
FA_percent_by_BG.aov = aov(FA_percent ~ HL_state * BG, data = FA_percent_by_BG.aov.data)

Parametric_Check(FA_percent_by_BG.aov)
# Normal
summary(FA_percent_by_BG.aov)
FA_percent_by_BG.aov.stats = tidy(FA_percent_by_BG.aov) %>% mutate(sig = stars.pval(p.value))

FA_percent_by_BG.aov.postHoc =
  tidy(TukeyHSD(FA_percent_by_BG.aov, "BG")) %>% mutate(sig = stars.pval(adj.p.value))

filter(FA_percent_by_BG.aov.postHoc, sig != " ")

# No effect of detail -----------------------------------------------------
## BBN Alone vs. Mixed duration for baseline vs. post-HL

###################
# TH Comparison
###################
detail.TH.aov.data = TH_table_detail %>% filter(Frequency == "0" & HL_state %in% c("baseline", "post-HL")) %>% filter(! is.na(TH))

detail.TH.aov.data$Gaus = LambertW::Gaussianize(detail.TH.aov.data$TH)[, 1]

detail.aov = aov(Gaus ~ detail * Duration * HL_state, data = detail.TH.aov.data)

Parametric_Check(detail.aov)
#TH very not spherical while Gaus is close

summary(detail.aov)

BBN.detail.stats = kruskal.test(TH ~ detail, data = detail.TH.aov.data)

# In both cases the detail i.e. mixed or alone is irrelevant (p > 0.11)
# so either combine or select BBN alone and 4-32 Mixed

BBN.detail.stats


###################
# Rxn Comparison
###################
detail.Rxn.aov.data = Rxn_table_over_TH_detail %>% filter(Frequency == "0" & HL_state %in% c("baseline", "post-HL")) %>%
  #limit to individuals with both states
  filter(rat_name %in% c("Orange11", "Orange12", "Green11", "Green12"))

detail.Rxn.aov.data$Gaus = LambertW::Gaussianize(detail.Rxn.aov.data$Rxn)[, 1]

detail.aov = aov(Gaus ~ detail * Duration * HL_state, data = detail.Rxn.aov.data)

Parametric_Check(detail.aov)
# Not even close to normal
summary(detail.aov)

# By intensity
detail.Rxn.aov.data %>% group_by(Intensity, Duration) %>%
  summarise(Comparison = kruskal.test(Rxn ~ detail)$data.name,
            P.unadj = kruskal.test(Rxn ~ detail)$p.value,
            .groups = "drop") %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj)) %>%
  filter(sig != " ")

# By HL_state
# No effect of Mixed vs. Alone on dprime across Hearing Loss
detail.Rxn.postHoc =
detail.Rxn.aov.data %>%
  # filter(rat_name != "Green12") %>%
  group_by(HL_state, Duration) %>%
  summarise(Comparison = kruskal.test(Rxn ~ detail)$data.name,
            P.unadj = kruskal.test(Rxn ~ detail)$p.value,
            .groups = "drop") %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj))


# By HL_state
detail.Rxn.aov.data %>% group_by(detail, Duration) %>%
  summarise(Comparison = kruskal.test(Rxn ~ HL_state)$data.name,
            P.unadj = kruskal.test(Rxn ~ HL_state)$p.value,
            .groups = "drop") %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj)) %>%
  filter(sig != " ")

detail.Rxn.postHoc

BBN.detail.Rxn.stats = kruskal.test(Rxn ~ detail, data = detail.Rxn.aov.data)

# In both cases the detail i.e. mixed or alone is irrelevant (p > 0.11)
# so either combine or select BBN alone and 4-32 Mixed

BBN.detail.Rxn.stats

###################
# d' Comparison
###################
detail.dprime.aov.data = dprime_detail_table %>% filter(Frequency == "0" & HL_state %in% c("baseline", "post-HL")) %>%
  #limit to individuals with both states
  filter(rat_name %in% c("Orange11", "Orange12", "Green11", "Green12"))

detail.dprime.aov.data$Gaus = LambertW::Gaussianize(detail.dprime.aov.data$dprime)[, 1]

detail.dprime.aov = aov(Gaus ~ detail * Duration * HL_state, data = detail.dprime.aov.data)

Parametric_Check(detail.dprime.aov)
# not close
# summary(detail.dprime.aov)

BBN.detail.dprime.stats = kruskal.test(dprime ~ detail, data = detail.dprime.aov.data)

BBN.duration.dprime.stats = kruskal.test(dprime ~ Duration, data = detail.dprime.aov.data)

BBN.HL.dprime.stats = kruskal.test(dprime ~ HL_state, data = detail.dprime.aov.data)

# Post-Hoc Dunns testing
# By HL_state
  detail.dprime.aov.data %>%
  group_by(Duration, HL_state) %>%
  summarise(dunns = tryCatch(FSA::dunnTest(dprime ~ interaction(detail), method = "none")$res, error = function(err){NA}),
            .groups = "drop") %>%
  unnest_wider(dunns) %>%
  filter(! is.na(Comparison)) %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj)) %>%
    filter(! sig %in% c(" "))

# In both cases the detail i.e. mixed or alone is irrelevant (p > 0.11)
# so either combine or select BBN alone and 4-32 Mixed

BBN.detail.dprime.stats


# Step size effect --------------------------------------------------------

step_size.dprime.aov.data = dprime_5step_table %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>% filter(Duration == 50) %>% filter(BG_type != "White") %>% filter(step_size != Inf) %>%
  mutate(step_size = as.factor(step_size))

step_size.dprime.aov.data$Gaus = LambertW::Gaussianize(step_size.dprime.aov.data$dprime)[, 1]

step_size.aov = aov(Gaus ~ step_size * Frequency * Intensity * HL_state, data = step_size.dprime.aov.data)

Parametric_Check(step_size.aov)
# not close
# summary(step_size.aov)

BBN.step_size.dprime.stats = kruskal.test(dprime ~ step_size, data = step_size.dprime.aov.data)
BBN.step_size.dprime.stats

# Do repeated Kruskal tests
step_size_factors = c("step_size", "Frequency", "HL_state", "BG_Intensity")

kruskal_results_dprime_step_size = lapply(step_size_factors, function(x) kruskal.test(reformulate(x, "dprime"), data = step_size.dprime.aov.data))

# convert to table
kruskal_results_dprime_step_size =
  as_tibble(do.call(rbind, kruskal_results_dprime_step_size)) %>%
  transmute(test = as.character(method),
            model = as.character(data.name),
            statistic = as.numeric(statistic),
            df = as.numeric(parameter),
            p.adj = p.adjust(p.value, "BH", n = 4),
            sig = stars.pval(p.adj))

kruskal_results_dprime_step_size

# Post-Hoc testing

PostHoc_Dunns_Testing <- function(df){
  df %>%
    group_by(Frequency, Intensity, BG_Intensity) %>%
    summarise(dunns = tryCatch(FSA::dunnTest(dprime ~ step_size, method = "none")$res, error = function(err){NA}),
              .groups = "drop") %>%
    unnest_wider(dunns) %>%
    filter(! is.na(Comparison)) %>%
    mutate(P.adj = p.adjust(P.unadj, "BH"),
           sig = stars.pval(P.adj)) # must escape the wild card . to be an actual period
}


step_size_postHoc_complete =
  step_size.dprime.aov.data %>%
  PostHoc_Dunns_Testing()

step_size_postHoc_complete %>% filter(! sig %in% c(" ", "."))


#######################
# Rxn Full Comparison
#######################

step_size.Rxn.aov.data = Rxn_5step_table %>% filter(Frequency == "0" & HL_state %in% c("baseline", "post-HL")) %>%
  # only rats with 10step on the 5s can be used as otherwise the ranges are different so we expect a difference
  filter(rat_name %in% c("Orange11", "Orange12", "Green11", "Green 12"))

step_size.Rxn.aov.data$Gaus = LambertW::Gaussianize(step_size.Rxn.aov.data$Rxn)[, 1]

step_size.aov = aov(Gaus ~ step_size * Frequency * Intensity * Duration * HL_state, data = step_size.Rxn.aov.data)

Parametric_Check(step_size.aov)
#TH very not spherical while Gaus is close

summary(step_size.aov)

step_size.Rxn.aov.data %>% group_by(Frequency, Duration) %>%
  summarise(P.unadj = kruskal.test(Rxn ~ step_size)$p.value, .groups = "drop") %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj))

BBN.step_size.Rxn.stats = kruskal.test(Rxn ~ step_size, data = step_size.Rxn.aov.data %>% filter(Duration == 300))

# In both cases the step_size i.e. mixed or alone is irrelevant (p > 0.11)
# so either combine or select BBN alone and 4-32 Mixed

BBN.step_size.Rxn.stats




# TH Change -------------------------------------------------------------

TH_averages = TH_table %>%
  group_by(Frequency, Duration, HL_state, BG_type, BG_Intensity) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")

BBN_TH_avg_table =
  TH_averages %>%
    filter(Frequency == 0) %>%
    # Combine Frequency and Duration to create a single key column
    mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
    # Remove redundant columns
    select(-all_of(c("Frequency", "Duration"))) %>%
    spread(key, TH)

BBN_TH_avg_table

TH_avg_50ms_table =
  TH_averages %>%
    # filter(Frequency != 0) %>%
    filter(Duration == 50) %>%
    spread(Frequency, TH)

TH_avg_50ms_table

# TH analysis -------------------------------------------------------------

TH.aov.data = TH_table %>% filter(! is.na(TH)) %>% filter(Duration == 50) %>% filter(BG_type != "White")

TH.aov.data$Gaus = LambertW::Gaussianize(TH.aov.data$TH)[, 1]

TH.aov = aov(Gaus ~ Frequency * HL_state * BG_type * BG_Intensity,
             data = TH.aov.data)

Parametric_Check(TH.aov)
#Not close
# summary(TH.aov)

# Do repeated Kruskal tests
TH_factors = c("Frequency", "HL_state", "BG_type", "BG_Intensity")

kruskal_results_TH = lapply(TH_factors, function(x) kruskal.test(reformulate(x, "TH"), data = TH.aov.data))

# convert to table
kruskal_results_TH =
  as_tibble(do.call(rbind, kruskal_results_TH)) %>%
  transmute(test = as.character(method),
            model = as.character(data.name),
            statistic = as.numeric(statistic),
            df = as.numeric(parameter),
            p.adj = p.adjust(p.value, "BH", n = 4),
            sig = stars.pval(p.adj))

kruskal_results_TH

# Post-Hoc testing

PostHoc_Dunns_Testing <- function(df){
  df %>%
  group_by(Frequency, BG_Intensity) %>%
    summarise(dunns = FSA::dunnTest(TH ~ HL_state, method = "none")$res,
              .groups = "drop") %>%
    unnest_wider(dunns) %>%
    mutate(P.adj = p.adjust(P.unadj, "BH"),
           sig = stars.pval(P.adj)) # must escape the wild card . to be an actual period
}


TH_postHoc_complete =
  TH.aov.data %>%
  PostHoc_Dunns_Testing()

TH_postHoc_timeline =
  TH.aov.data %>%
  filter(BG_type == "None") %>%
  PostHoc_Dunns_Testing()

TH_postHoc_BG =
  TH.aov.data %>%
  filter(Frequency != 0) %>%
  # filter(rat_ID %in% c(22, 15, 19, 20, 24, 17)) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  filter(BG_type %in% c("None", "Pink")) %>%
  PostHoc_Dunns_Testing() %>%
  print


filter(TH_postHoc_complete, sig != " ")


# dprime testing ----------------------------------------------------------

dprime.aov.data = dprime_summary_table %>%
  filter(str_detect(Intensity, pattern = "0$")) %>% filter(HL_state %in% c("baseline", "post-HL")) %>% filter(Duration == 50) %>% filter(BG_type != "White") %>%
  # drop Intensitites with only one rat
  filter(! (Frequency %in% c(4, 8) && Intensity == 0))

dprime.aov.data$Gaus = LambertW::Gaussianize(dprime.aov.data$dprime)[, 1]

dprime.aov = aov(Gaus ~ Frequency * HL_state * BG_type * BG_Intensity,
              data = dprime.aov.data)

Parametric_Check(dprime.aov)
#dprime very not spherical
# summary(dprime.aov)

# Do repeated Kruskal tests
dprime_factors = c("Frequency", "Intensity", "HL_state", "BG_Intensity")

kruskal_results_dprime = lapply(dprime_factors, function(x) kruskal.test(reformulate(x, "dprime"), data = dprime.aov.data))

# convert to table
kruskal_results_dprime =
  as_tibble(do.call(rbind, kruskal_results_dprime)) %>%
  transmute(test = as.character(method),
            model = as.character(data.name),
            statistic = as.numeric(statistic),
            df = as.numeric(parameter),
            p.unadj = as.numeric(p.value),
            p.adj = p.adjust(p.value, "BH"),
            sig = stars.pval(p.adj))

kruskal_results_dprime

# Post-Hoc testing
dprime_postHoc_complete =
  dprime.aov.data %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  filter(Frequency != 0) %>%
  group_by(Frequency, Intensity) %>%
  summarise(dunns = FSA::dunnTest(dprime ~ BG_Intensity, method = "none")$res,
            .groups = "drop") %>%
  unnest_wider(dunns) %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj))


filter(dprime_postHoc_complete, sig != " ")


# Reaction time -----------------------------------------------------------


# BBN Reaction time
BBN.Rxn.aov.data = Rxn_table_over_TH_detail %>% filter(Frequency == 0 & detail == "Alone") %>%
  filter(HL_state %in% c("baseline", "post-HL"))

BBN.Rxn.aov.data$Gaus = LambertW::Gaussianize(BBN.Rxn.aov.data$Rxn)[, 1]


BBN.Rxn.aov = aov(Gaus ~ Intensity * Duration * HL_state,
                  data = BBN.Rxn.aov.data)

Parametric_Check(BBN.Rxn.aov)
#Rxn very not spherical
summary(BBN.Rxn.aov)

# Do repeated Kruskal tests
Rxn_factors = c("Intensity", "Duration", "HL_state")

kruskal_results_BBN_Rxn = lapply(Rxn_factors, function(x) kruskal.test(reformulate(x, "Rxn"), data = BBN.Rxn.aov.data))

# convert to table
kruskal_results_BBN_Rxn =
  as_tibble(do.call(rbind, kruskal_results_BBN_Rxn)) %>%
  transmute(test = as.character(method),
            model = as.character(data.name),
            statistic = as.numeric(statistic),
            df = as.numeric(parameter),
            p.adj = p.adjust(p.value, "BH", n = 4),
            sig = stars.pval(p.adj))

kruskal_results_BBN_Rxn






Rxn.aov.data = Rxn_table_over_TH %>% filter(str_detect(Intensity, pattern = "0$")) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>% filter(Duration == 50) %>% filter(BG_type != "White")

Rxn.aov.data$Gaus = LambertW::Gaussianize(Rxn.aov.data$Rxn)[, 1]

Rxn.aov = aov(Gaus ~ Frequency * HL_state * BG_type * BG_Intensity,
             data = Rxn.aov.data)

Parametric_Check(Rxn.aov)
#Rxn very not spherical while Gaus is close

summary(Rxn.aov)

# Do repeated Kruskal tests
Rxn_factors = c("Frequency", "Intensity", "HL_state", "BG_type", "BG_Intensity")

kruskal_results_Rxn = lapply(Rxn_factors, function(x) kruskal.test(reformulate(x, "Rxn"), data = Rxn.aov.data))

# convert to table
kruskal_results_Rxn =
  as_tibble(do.call(rbind, kruskal_results_Rxn)) %>%
  transmute(test = as.character(method),
            model = as.character(data.name),
            statistic = as.numeric(statistic),
            df = as.numeric(parameter),
            p.unadj = as.numeric(p.value),
            p.adj = p.adjust(p.value, "BH"),
            sig = stars.pval(p.adj))

kruskal_results_Rxn

# Post-Hoc testing


Rxn_postHoc_complete =
  Rxn.aov.data %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  filter(Frequency != 0) %>%
  group_by(Frequency, Intensity) %>%
  summarise(dunns = FSA::dunnTest(Rxn ~ BG_Intensity, method = "none")$res,
            .groups = "drop") %>%
  unnest_wider(dunns) %>%
  mutate(P.adj = p.adjust(P.unadj, "BH"),
         sig = stars.pval(P.adj))


filter(Rxn_postHoc_complete, sig != " ")

