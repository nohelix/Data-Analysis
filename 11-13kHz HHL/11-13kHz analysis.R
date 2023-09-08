# Package loading ---------------------------------------------------------

# data loading/manipulation
library(readxl); library(tidyverse); library(magrittr); library(dplyr); library(tidyr); library(broom); library(glue)

# Analysis
library(FSA); library(LambertW)

# Data visualization
library(ggplot2); library(forcats); library(gtools); library(ggpmisc)

# Data Export
library(data.table)


# Variables ---------------------------------------------------------------
ProjectFolder = "C:/Users/Noelle/Box/ABR recordings/ABR Results/11-13kHz HHL"

SaveFolder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Paper Drafts/2023 11-13kHz noise expsure"

Shams = c("244", "247")


# Load Data ---------------------------------------------------------------

Graph = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                   sheet = "4-32kHz")

Graph_detailed = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                            sheet = "4-48kHz", na = "NA")


# Functions ---------------------------------------------------------------

Parametric_Test <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05

  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}

}

# Select Data ------------------------------------------------------------
# Can be summarized or not

# Combine 4-32 & 4-48 data sets
Wave1_data =
  Graph_detailed %>%
  bind_rows(Graph) %>%
  # Get average for each rat across both ears
  summarise(RMS = mean(RMS, na.rm = TRUE),
            W1.amp = mean(`W1 amp (uV)`, na.rm = TRUE),
            W1.lat = mean(`W1 latency (ms)`, na.rm = TRUE),
            .by = c(ID, Condition, Freq, Inten)) %>%
  # Label frequency and order
  mutate(Freq = case_when(Freq != "0" ~ str_glue("{Freq} kHz"),
                          Freq == "0" ~ "BBN",
                          .default = "ERROR") %>%
           ordered(c("4 kHz", "6 kHz", "8 kHz", "12 kHz", "16 kHz",
                     "24 kHz", "32 kHz", "48 kHz", "BBN"))) %>%
  # fix Timepoint (currently called Condition)
  mutate(Condition = str_replace_all(Condition, "[45] week", "4-5 week"),
         Condition = ordered(Condition, c("Baseline", "Hearing Loss", "1 day",
                                          "1 week", "2 week", "4-5 week"))) %>%
  rename(Timepoint = Condition) %>%
  # Create treatment (called Condition for ease of reading)
  mutate(Condition = if_else(ID %in% Shams, "Sham", "TTS"))

# See which rats have which data
Wave1_data %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  summarise(n_ID = length(unique(ID)), ID = paste(unique(ID), collapse = ", "),
            .by = c(Timepoint, Freq, Condition)) %>%
  spread(Freq, n_ID) %>%
  relocate(ID, .after = BBN) %>%
  print

# Anova data prep
AOV.data =
  Wave1_data

# Attempt to transform to normal
AOV.data$RMS.Gaus = LambertW::Gaussianize(AOV.data$RMS)[, 1]
AOV.data$W1.amp.Gaus = LambertW::Gaussianize(AOV.data$W1.amp)[, 1]
AOV.data$W1.lat.Gaus = LambertW::Gaussianize(AOV.data$W1.lat)[, 1]

# TH Calc -----------------------------------------------------------------

# Calculate from lowest Intensity
Wave1_TH =
  Wave1_data %>%
  group_by(ID, Timepoint, Freq) %>%
  do(filter(., Inten == min(Inten))) %>%
  ungroup() %>%
  summarise(TH = mean(Inten, na.rm = TRUE),
            .by = c(ID, Condition, Timepoint, Freq)) %>%
  mutate(Timepoint = ordered(Timepoint, c("Baseline", "Hearing Loss", "1 day",
                                          "1 week", "2 week", "4-5 week")))


Wave1_TH %>%
  group_by(Freq) %>%
  do(dplyr::filter(., TH == max(TH)))


# TH Stats ----------------------------------------------------------------

TH.aov <- aov(TH ~ Condition * Freq * Timepoint,
               data = Wave1_TH)

Parametric_Test(TH.aov)
# summary(RMS.aov)

# Kruskal Testing - very non-normal
TH.stats <-
  lapply(c("Freq", "Timepoint", "Condition"),
         function(x) kruskal.test(reformulate(x, "TH"), data = Wave1_TH)) %>%
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value))

TH.stats

TH.postHoc =
  Wave1_TH %>%
  group_by(Freq) %>%
  do(
    FSA::dunnTest(TH ~ interaction(Timepoint, Condition),
                  data = .,
                  method = "none") %>% print()
  ) %>% dplyr::select(-P.adj) %>%
  mutate(P.adj = p.adjust(P.unadj, "bonf"),
         Sig = gtools::stars.pval(P.adj),
         TimePoint1 = str_extract(Comparison, pattern = '^.+?(?=\\.)') %>%
                          ordered(c("Baseline", "Hearing Loss", "1 day",
                                    "1 week", "2 week", "4-5 week")),
         Condition1 = str_extract(Comparison, pattern = '(?<=\\.).+?(?= -)'),
         TimePoint2 = str_extract(Comparison, pattern = '(?<= - ).+?(?=\\.)') %>%
                         ordered(c("Baseline", "Hearing Loss", "1 day",
                                   "1 week", "2 week", "4-5 week")),
         Condition2 = str_extract(Comparison, pattern = '(?<= - ).+?$') %>% str_extract(., pattern = '(?<=\\.).+?$')) %>%
  dplyr::select(Freq, TimePoint1, TimePoint2, Condition1, Condition2, Z, P.adj, Sig)

TH.postHoc %>%
  filter(! Sig %in% c(" ", "."))


# TH Graphs ---------------------------------------------------------------

TH.annotation =
  TH.postHoc %>%
  ungroup %>%
  filter(! Sig %in% c(" ", ".")) %>%
  mutate(Timepoint_start = case_when(TimePoint1 < TimePoint2 ~ TimePoint1,
                                     TimePoint2 < TimePoint1 ~ TimePoint2,
                                     .default = NA),
         Timepoint_end = case_when(TimePoint1 < TimePoint2 ~ TimePoint2,
                                   TimePoint2 < TimePoint1 ~ TimePoint1,
                                   .default = NA),
         x = (as.numeric(Timepoint_start) + as.numeric(Timepoint_end))/2,
         Condition = Condition1
  ) %>%
  left_join(Wave1_TH %>%
              group_by(Freq) %>%
              filter(Timepoint == "Hearing Loss") %>%
              do(dplyr::filter(., TH == max(TH))) %>%
              dplyr::select(Freq, TH),
            relationship = "many-to-many",
            by = c("Freq" = "Freq")) %>%
  mutate(Condition = if_else(Condition == "Sham", "Sham", "11-13kHz\nNoise Exposure")) %>%
  unique


Wave1_TH %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  mutate(Condition = if_else(Condition == "Sham", "Sham", "11-13kHz\nNoise Exposure")) %>%
  ggplot(aes(x = Timepoint, y = TH, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1.2) +
  # add lines for significance
  geom_segment(data = TH.annotation %>% filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")),
               aes(x = Timepoint_start, xend = Timepoint_end, y = TH - 5, yend = TH - 5),
               position = position_jitter(height = 5, width = 0),
               show.legend = FALSE) +
  # add significance stars
  geom_text(data = TH.annotation %>% filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")),
             aes(label = Sig, x = x, y = TH), size = 7, color = "black",
             show.legend = FALSE) +
  labs(color = "Treatment",
       x = "", y = "Threshold (+/- SE)",
       caption = "Significant threshold shifts only between TTS timepoints") +
  facet_wrap( ~ Freq, scale = "fixed", ncol = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
    legend.position = c(0.8, 0.1),
    legend.background=element_blank()
  )

ggsave(filename = "TH_shifts.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 10, height = 10, units = "in", dpi = 300)


# Overall Graph -----------------------------------------------------------

# BBN
# RMS, W1 latency (ms), W1 amp (uV)
Wave1_data %>%
  filter(Inten == "80" & Freq %in% c("BBN")) %>%
  rename(`W1 amp (uV)` = W1.amp,
         `W1 latency (ms)` = W1.lat) %>%
  mutate(Condition = if_else(Condition == "Sham", "Sham", "11-13kHz\nNoise Exposure")) %>%
  gather(key = "ABR", value = "value", RMS, `W1 amp (uV)`, `W1 latency (ms)`) %>%
  ggplot(aes(x = Timepoint, y = value, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 2) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1) +
  labs(title = "80dB Intensity Only for Broadband",
       color = "Treatment",
       y = "") +
  facet_wrap( ~ ABR, scales = "free", nrow = 3) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

ggsave(filename = "ABR_BBN_80dB.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# 16 kHz
Wave1_data %>%
  filter(Inten == "80" & Freq %in% c("16 kHz")) %>%
  rename(`W1 amp (uV)` = W1.amp,
         `W1 latency (ms)` = W1.lat) %>%
  mutate(Condition = if_else(Condition == "Sham", "Sham", "11-13kHz\nNoise Exposure")) %>%
  gather(key = "ABR", value = "value", RMS, `W1 amp (uV)`, `W1 latency (ms)`) %>%
  ggplot(aes(x = Timepoint, y = value, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 2) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1) +
  labs(title = "80dB Intensity Only for 16 kHz tone",
       color = "Treatment",
       y = "") +
  facet_wrap( ~ ABR, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

ggsave(filename = "ABR_16kHz_80dB.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# RMS ANOVA -------------------------------------------------------------------

RMS.aov <- aov(RMS.Gaus ~ Condition * Freq * Inten * Timepoint,
               data = AOV.data)

Parametric_Test(RMS.aov)
# summary(RMS.aov)

# Kruskal Testing - very non-normal
RMS.stats <-
  lapply(c("Freq", "Inten", "Timepoint", "Condition"),
         function(x) kruskal.test(reformulate(x, "RMS"), data = AOV.data)) %>%
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "bonf"),
           sig = gtools::stars.pval(adj.p.value))

RMS.stats

RMS.postHoc =
  AOV.data %>%
  group_by(Freq) %>%
  do(
    FSA::dunnTest(RMS ~ interaction(Timepoint, Condition),
                  data = .,
                  method = "none") %>% print()
  ) %>% dplyr::select(-P.adj) %>%
  mutate(P.adj = p.adjust(P.unadj, "bonf"),
         Sig = gtools::stars.pval(P.adj),
         TimePoint1 = str_extract(Comparison, pattern = '^.+?(?=\\.)'),
         Condition1 = str_extract(Comparison, pattern = '(?<=\\.).+?(?= -)'),
         TimePoint2 = str_extract(Comparison, pattern = '(?<= - ).+?(?=\\.)'),
         Condition2 = str_extract(Comparison, pattern = '(?<= - ).+?$') %>% str_extract(., pattern = '(?<=\\.).+?$')) %>%
  dplyr::select(Freq, TimePoint1, TimePoint2, Condition1, Condition2, Z, P.adj, Sig)

RMS.postHoc %>%
    filter(! Sig %in% c(" ", ".")) %>%
    View


# RMS Graph ---------------------------------------------------------------

# All Frequencies
Wave1_data  %>%
  filter(Condition == "TTS") %>%
  ggplot(aes(x = Timepoint, y = RMS, color = Inten, group = Inten)) +
  geom_smooth(aes(color = 100, group = 100)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line") +
  # scale_x_continuous(limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) +
  # labs(x = "Sound Intensity (dB)",
  # y = "Signal-to-Noise Ratio (RMS)") +
  facet_wrap( ~ Freq, nrow = 2) +
  labs(title = "") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

# Input-Output functions
Wave1_data  %>%
  filter(Condition == "TTS") %>%
  filter(Timepoint %in% c("Baseline", "Hearing Loss", "1 day", "2 week", "4-5 week")) %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  ggplot(aes(x = Inten, y = RMS, color = Timepoint, group = Timepoint)) +
  geom_smooth(se = FALSE) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  scale_x_continuous(n.breaks = 10) +
  labs(Title = "Input-Output for Wave 1 Amplitude",
       x = "Sound Intensity (dB)",
       y = "Signal to Noise (RMS)") +
  facet_wrap( ~ Freq, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    # axis.title.x = element_text(hjust = 0.45),
    legend.position = c(0.8, 0.1),
    legend.background=element_blank()
  )

ggsave(filename = "ABR_RMS_IO.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# W1 Change data ----------------------------------------------------------

# Percent (%) change from baseline
Wave1_percent_change =
  Wave1_data %>%
  group_by(ID, Condition, Freq, Inten) %>%
  # Timepoints = c("Baseline", "Hearing Loss", "1 day", "1 week", "2 week", "4-5 week")
  do(
    # Hearing Loss
    RMS_HL = filter(., Timepoint == "Hearing Loss")$RMS / filter(., Timepoint == "Baseline")$RMS,
    W1.amp_HL = filter(., Timepoint == "Hearing Loss")$W1.amp / filter(., Timepoint == "Baseline")$W1.amp,
    W1.lat_HL = filter(., Timepoint == "Hearing Loss")$W1.lat / filter(., Timepoint == "Baseline")$W1.lat,
    # 1 day
    RMS_1d = filter(., Timepoint == "1 day")$RMS / filter(., Timepoint == "Baseline")$RMS,
    W1.amp_1d = filter(., Timepoint == "1 day")$W1.amp / filter(., Timepoint == "Baseline")$W1.amp,
    W1.lat_1d = filter(., Timepoint == "1 day")$W1.lat / filter(., Timepoint == "Baseline")$W1.lat,
    # 1 week
    RMS_1w = filter(., Timepoint == "1 week")$RMS / filter(., Timepoint == "Baseline")$RMS,
    W1.amp_1w = filter(., Timepoint == "1 week")$W1.amp / filter(., Timepoint == "Baseline")$W1.amp,
    W1.lat_1w = filter(., Timepoint == "1 week")$W1.lat / filter(., Timepoint == "Baseline")$W1.lat,
    # 2 weeks
    RMS_2w = filter(., Timepoint == "2 week")$RMS / filter(., Timepoint == "Baseline")$RMS,
    W1.amp_2w = filter(., Timepoint == "2 week")$W1.amp / filter(., Timepoint == "Baseline")$W1.amp,
    W1.lat_2w = filter(., Timepoint == "2 week")$W1.lat / filter(., Timepoint == "Baseline")$W1.lat,
    # 4-5 weeks
    RMS_5w = filter(., Timepoint == "4-5 week")$RMS / filter(., Timepoint == "Baseline")$RMS,
    W1.amp_5w = filter(., Timepoint == "4-5 week")$W1.amp / filter(., Timepoint == "Baseline")$W1.amp,
    W1.lat_5w = filter(., Timepoint == "4-5 week")$W1.lat / filter(., Timepoint == "Baseline")$W1.lat,
  ) %>%
  ungroup() %>%
  mutate(Freq = ordered(Freq, c("4 kHz", "6 kHz", "8 kHz", "12 kHz", "16 kHz",
                                "24 kHz", "32 kHz", "48 kHz", "BBN")),
         # Hearing Loss
         RMS_HL = as.numeric(RMS_HL) * 100,
         W1.amp_HL = as.numeric(W1.amp_HL) * 100,
         W1.lat_HL = as.numeric(W1.lat_HL) * 100,
         # 1 day
         RMS_1d = as.numeric(RMS_1d) * 100,
         W1.amp_1d = as.numeric(W1.amp_1d) * 100,
         W1.lat_1d = as.numeric(W1.lat_1d) * 100,
         # 1 week
         RMS_1w = as.numeric(RMS_1w) * 100,
         W1.amp_1w = as.numeric(W1.amp_1w) * 100,
         W1.lat_1w = as.numeric(W1.lat_1w) * 100,
         # 2 week
         RMS_2w = as.numeric(RMS_2w) * 100,
         W1.amp_2w = as.numeric(W1.amp_2w) * 100,
         W1.lat_2w = as.numeric(W1.lat_2w) * 100,
         # 4-5 weeks
         RMS_5w = as.numeric(RMS_5w) * 100,
         W1.amp_5w = as.numeric(W1.amp_5w) * 100,
         W1.lat_5w = as.numeric(W1.lat_5w) * 100,
         )

# Change in ABR W1 attributes from either:
#     avr = average of last 2 timepionts or final = just the last timepoint
Wave1_final_change =
  Wave1_data %>%
    group_by(ID, Condition, Freq, Inten) %>%
    do(
      # Average 2 & 5 week time points
      RMS_avg_change = mean(filter(., Timepoint %in% c("2 week", "4-5 week"))$RMS) -
        filter(., Timepoint == "Baseline")$RMS,
      W1.amp.avg_change = mean(filter(., Timepoint %in% c("2 week", "4-5 week"))$W1.amp) -
        filter(., Timepoint == "Baseline")$W1.amp,
      W1.lat.avg_change = mean(filter(., Timepoint %in% c("2 week", "4-5 week"))$W1.lat) -
        filter(., Timepoint == "Baseline")$W1.lat,
      # only use the final time point
      RMS_final_change = filter(., Timepoint == "4-5 week")$RMS - filter(., Timepoint == "Baseline")$RMS,
      W1.amp.final_change = filter(., Timepoint == "4-5 week")$W1.amp - filter(., Timepoint == "Baseline")$W1.amp,
      W1.lat.final_change = filter(., Timepoint == "4-5 week")$W1.lat - filter(., Timepoint == "Baseline")$W1.lat
    ) %>%
    ungroup() %>%
    mutate(Freq = ordered(Freq, c("4 kHz", "6 kHz", "8 kHz", "12 kHz", "16 kHz",
                       "24 kHz", "32 kHz", "48 kHz", "BBN")),
           RMS_avg_change = as.numeric(RMS_avg_change),
           W1.amp.avg_change = as.numeric(W1.amp.avg_change),
           W1.lat.avg_change = as.numeric(W1.lat.avg_change),
           RMS_final_change = as.numeric(RMS_final_change),
           W1.amp.final_change = as.numeric(W1.amp.final_change),
           W1.lat.final_change = as.numeric(W1.lat.final_change))


# W1 amp ANOVA ----------------------------------------------------------------

W1.amp.aov <- aov(W1.amp.Gaus ~ Condition * Freq * Inten * Timepoint,
                  data = AOV.data)

Parametric_Test(W1.amp.aov)
# summary(RMS.aov)

# Kruskal Testing - very non-normal
W1.amp.stats <-
  lapply(c("Freq", "Inten", "Timepoint", "Condition"),
         function(x) kruskal.test(reformulate(x, "W1.amp"), data = AOV.data)) %>%
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value))

W1.amp.stats

W1.amp.postHoc =
  AOV.data %>%
  group_by(Freq) %>%
  do(
    FSA::dunnTest(W1.amp ~ interaction(Timepoint, Condition),
                  data = .,
                  method = "none") %>% print()
  ) %>% dplyr::select(-P.adj) %>%
  mutate(P.adj = p.adjust(P.unadj, "bonf"),
         Sig = gtools::stars.pval(P.adj),
         TimePoint1 = str_extract(Comparison, pattern = '^.+?(?=\\.)'),
         Condition1 = str_extract(Comparison, pattern = '(?<=\\.).+?(?= -)'),
         TimePoint2 = str_extract(Comparison, pattern = '(?<= - ).+?(?=\\.)'),
         Condition2 = str_extract(Comparison, pattern = '(?<= - ).+?$') %>% str_extract(., pattern = '(?<=\\.).+?$')) %>%
  dplyr::select(Freq, TimePoint1, TimePoint2, Condition1, Condition2, Z, P.adj, Sig)

W1.amp.postHoc %>%
  filter(! Sig %in% c(" ", ".")) %>%
  View


# W1 Amp Graph ------------------------------------------------------------

# Input-Output functions
Wave1_data  %>%
  filter(Condition == "TTS") %>%
  filter(Timepoint %in% c("Baseline", "Hearing Loss", "1 day", "2 week", "4-5 week")) %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  ggplot(aes(x = Inten, y = W1.amp, color = Timepoint, group = Timepoint)) +
    geom_smooth(se = FALSE) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                 fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                 geom = "errorbar", width = 0.1) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
    scale_x_continuous(n.breaks = 10) +
    labs(Title = "Input-Output for Wave 1 Amplitude",
         x = "Sound Intensity (dB)",
         y = "Amplitude (\u03BCV)") +
    facet_wrap( ~ Freq, scales = "free_x", ncol = 2) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      # axis.title.x = element_text(hjust = 0.45),
      legend.position = c(0.8, 0.1),
      legend.background=element_blank()
    )

ggsave(filename = "ABR_W1 amp_IO.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

Wave1_final_change  %>%
  filter(Condition == "TTS") %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  ggplot(aes(x = Inten, y = W1.amp.avg_change, color = Freq, group = Freq)) +
  geom_smooth(se = FALSE) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1,
               position = position_dodge(0.3)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3,
               position = position_dodge(0.3)) +
  # stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line",
  #              position = position_dodge(0.3)) +
  scale_x_continuous(n.breaks = 10) +
  labs(title = "Permanent change in Wave 1 Amplitute",
       x = "Sound Intensity (dB)",
       y = "Amplitude (\u03BCV)",
       color = "Frequency") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    # axis.title.x = element_text(hjust = 0.45),
    legend.position = c(0.15, 0.25),
    # legend.background=element_blank()
  )

ggsave(filename = "ABR_W1 amp_change.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 9, height = 6, units = "in", dpi = 300)

# % change
Wave1_percent_change %>%
  filter(Inten %in% c(40, 60, 80)) %>%
  filter(! Freq %in% c("BBN")) %>%
  # reformat for graphing:
  dplyr::select(1:4, starts_with("W1.amp")) %>%
  rename_with( ~ str_remove(., pattern = "W1.amp_")) %>%
  gather(key = "Timepoint", value = "W1.amp", 5:ncol(.)) %>%
  mutate(Timepoint = case_match(Timepoint, "HL" ~ "Hearing Loss", "1d" ~ "1 day", "1w" ~ "1 week", "2w" ~ "2 week", "5w" ~ "4-5 week"),
         Freq = str_remove(Freq, " kHz") %>% as.numeric(),
         Timepoint = if_else(ID %in% Shams & Timepoint == "4-5 week", "Sham", Timepoint) %>%
            ordered(c("Sham", "Hearing Loss", "1 day", "1 week", "2 week", "4-5 week")),
         # Condition = if_else(ID %in% Shams & Timepoint == "Sham", "TTS", Condition),
         ) %>%
  filter(Condition == "TTS") %>%
  filter(Timepoint %in% c("Baseline", "Hearing Loss", "1 day", "2 week", "4-5 week", "Sham")) %>%
  ggplot(aes(x = Freq, y = W1.amp, color = Timepoint, group = Timepoint)) +
    # shaded region of noise exposure
    annotate('rect', xmin = 11, xmax = 13, ymin = -Inf, ymax = Inf, alpha = .3, fill = "pink") +
    # 100% i.e. baseline levels
    geom_hline(yintercept = 100, linewidth = 1.5, linetype = "longdash", color = "darkred") +
    # annotate('text', label = "Baseline", x = 48, y = 110, color = "darkred") +
    geom_text(data = data.frame(Freq = 44, W1.amp = 110, Inten = 40, Timepoint = "Sham"),
              label = "Baseline", color = "darkred") +
    # geom_smooth(se = FALSE) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                 fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                 geom = "errorbar", width = 0.1,
                 position = position_dodge(0.1)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3,
                 position = position_dodge(0.1)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1,
                 position = position_dodge(0.1)) +
    scale_x_continuous(breaks = c(4, 6, 8, 12, 16, 24, 32, 48), trans = 'log2') +
    labs(title = "Change in Wave 1 Amplitute from Baseline following noise exposure",
         x = "Frequency (kHz)",
         y = "% change from Baseline",
         color = "Timepoint") +
    facet_wrap( ~ Inten, scales = "free", nrow = 3) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      # axis.title.x = element_text(hjust = 0.45),
      # legend.position = c(0.15, 0.25),
      # legend.background=element_blank()
    )

ggsave(filename = "ABR_W1 amp_percent change.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 9, units = "in", dpi = 300)

# W1 lat ANOVA ----------------------------------------------------------------

W1.lat.aov <- aov(W1.lat.Gaus ~ Condition * Freq * Inten * Timepoint,
                  data = AOV.data)

Parametric_Test(W1.lat.aov)
# summary(RMS.aov)

# Kruskal Testing - very non-normal
W1.lat.stats <-
  lapply(c("Freq", "Inten", "Timepoint", "Condition"),
         function(x) kruskal.test(reformulate(x, "W1.lat"), data = AOV.data)) %>%
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value))

W1.lat.stats

W1.lat.postHoc =
  AOV.data %>%
  group_by(Freq) %>%
  do(
    FSA::dunnTest(W1.lat ~ interaction(Timepoint, Condition),
                  data = .,
                  method = "none") %>% print()
  ) %>% dplyr::select(-P.adj) %>%
  mutate(P.adj = p.adjust(P.unadj, "bonf"),
         Sig = gtools::stars.pval(P.adj),
         TimePoint1 = str_extract(Comparison, pattern = '^.+?(?=\\.)'),
         Condition1 = str_extract(Comparison, pattern = '(?<=\\.).+?(?= -)'),
         TimePoint2 = str_extract(Comparison, pattern = '(?<= - ).+?(?=\\.)'),
         Condition2 = str_extract(Comparison, pattern = '(?<= - ).+?$') %>% str_extract(., pattern = '(?<=\\.).+?$')) %>%
  dplyr::select(Freq, TimePoint1, TimePoint2, Condition1, Condition2, Z, P.adj, Sig)

W1.lat.postHoc %>%
  filter(! Sig %in% c(" ", ".")) %>%
  View

# W1 lat Graph --------------------------------------------------------------

# Input-Output functions
Wave1_data  %>%
  filter(Condition == "TTS") %>%
  filter(Timepoint %in% c("Baseline", "Hearing Loss", "1 day", "2 week", "4-5 week")) %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  ggplot(aes(x = Inten, y = W1.lat, color = Timepoint, group = Timepoint)) +
  geom_smooth(se = FALSE) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  scale_x_continuous(n.breaks = 10) +
  labs(Title = "Input-Output for Wave 1 Latency",
       x = "Sound Intensity (dB)",
       y = "Latency (ms)") +
  facet_wrap( ~ Freq, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    # axis.title.x = element_text(hjust = 0.45),
    legend.position = c(0.8, 0.1),
    legend.background=element_blank()
  )

ggsave(filename = "ABR_W1 lat_IO.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

Wave1_final_change  %>%
  filter(Condition == "TTS") %>%
  filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
  ggplot(aes(x = Inten, y = W1.lat.avg_change, color = Freq, group = Freq)) +
  geom_smooth(se = FALSE) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1,
               position = position_dodge(0.3)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3,
               position = position_dodge(0.3)) +
  # stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line",
  #              position = position_dodge(0.3)) +
  scale_x_continuous(n.breaks = 10) +
  labs(title = "No permanent change in Wave 1 Latency",
       x = "Sound Intensity (dB)",
       y = "Latency (ms)",
       color = "Frequency") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    # axis.title.x = element_text(hjust = 0.45),
    legend.position = c(0.9, 0.8),
    # legend.background=element_blank()
  )

ggsave(filename = "ABR_W1 lat_change.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 9, height = 6, units = "in", dpi = 300)

# % change
Wave1_percent_change %>%
  filter(Inten %in% c(40, 60, 80)) %>%
  filter(! Freq %in% c("BBN")) %>%
  # reformat for graphing:
  dplyr::select(1:4, starts_with("W1.lat")) %>%
  rename_with( ~ str_remove(., pattern = "W1.lat_")) %>%
  gather(key = "Timepoint", value = "W1.lat", 5:ncol(.)) %>%
  mutate(Timepoint = case_match(Timepoint, "HL" ~ "Hearing Loss", "1d" ~ "1 day", "1w" ~ "1 week", "2w" ~ "2 week", "5w" ~ "4-5 week"),
         Freq = str_remove(Freq, " kHz") %>% as.numeric(),
         Timepoint = if_else(ID %in% Shams & Timepoint == "4-5 week", "Sham", Timepoint) %>%
           ordered(c("Sham", "Hearing Loss", "1 day", "1 week", "2 week", "4-5 week")),
         # Condition = if_else(ID %in% Shams & Timepoint == "Sham", "TTS", Condition),
  ) %>%
  filter(Condition == "TTS") %>%
  filter(Timepoint %in% c("Baseline", "Hearing Loss", "1 day", "2 week", "4-5 week", "Sham")) %>%
  ggplot(aes(x = Freq, y = W1.lat, color = Timepoint, group = Timepoint)) +
    # shaded region of noise exposure
    annotate('rect', xmin = 11, xmax = 13, ymin = -Inf, ymax = Inf, alpha = .3, fill = "pink") +
    # 100% i.e. baseline levels
    geom_hline(yintercept = 100, linewidth = 1.5, linetype = "longdash", color = "darkred") +
    # annotate('text', label = "Baseline", x = 48, y = 110, color = "darkred") +
    geom_text(data = data.frame(Freq = 44, W1.lat = 110, Inten = 40, Timepoint = "Sham"),
              label = "Baseline", color = "darkred") +
    # geom_smooth(se = FALSE) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                 fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                 geom = "errorbar", width = 0.1,
                 position = position_dodge(0.3)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3,
                 position = position_dodge(0.3)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1,
                 position = position_dodge(0.3)) +
    scale_x_continuous(breaks = c(4, 6, 8, 12, 16, 24, 32, 48), trans = 'log2') +
    labs(title = "Change in Wave 1 Latency from Baseline following noise exposure",
         x = "Frequency (kHz)",
         y = "% change from Baseline",
         color = "Timepoint") +
    facet_wrap( ~ Inten, scales = "fixed", nrow = 3) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      # axis.title.x = element_text(hjust = 0.45),
      # legend.position = c(0.15, 0.25),
      # legend.background=element_blank()
    )

ggsave(filename = "ABR_W1 latency_percent change.jpg",
       path = SaveFolder,
       plot = last_plot(),
       width = 8, height = 9, units = "in", dpi = 300)


# Export Data -------------------------------------------------------------

fwrite(Wave1_data, file = glue("{SaveFolder}/Wave1_data.csv"), row.names = FALSE)
fwrite(Wave1_TH, file = glue("{SaveFolder}/Wave1_threshold_data.csv"), row.names = FALSE)
fwrite(Wave1_final_change, file = glue("{SaveFolder}/Wave1_permenant_change_data.csv"), row.names = FALSE)
fwrite(Wave1_percent_change, file = glue("{SaveFolder}/Wave1_percent_change_data.csv"), row.names = FALSE)
