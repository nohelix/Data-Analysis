# Package loading ---------------------------------------------------------

# data loading/manipulation
library(readxl); library(tidyverse); library(magrittr); library(dplyr); library(tidyr); library(broom); library(glue)

# Analysis
library(FSA); library(LambertW)

# Data visualization
library(ggplot2); library(forcats); library(gtools); library(ggpmisc); library(hrbrthemes)


# Variables ---------------------------------------------------------------
ProjectFolder = "C:/Users/Noelle/Box/ABR recordings/ABR Results/11-13kHz HHL"

Shams = c("244", "247")


# Load Data ---------------------------------------------------------------

Graph = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                   sheet = "4-32kHz")

Graph_detailed = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                            sheet = "4-48kHz", na = "NA")

ABR.TH_short = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                           sheet = "4-32_Thresholds", na = "NA")

ABR.TH_Detailed = read_excel(glue("{ProjectFolder}/Summary.xlsx"),
                              sheet = "4-48_Thresholds", na = "NA")


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
  Wave1_data %>%
  # filter(Condition == "TTS") %>%
  # Normalize to Baseline -
  # TODO: check if baseline exists and discard if it doesn't?
  group_by(ID, Freq, Inten, Condition) %>%
  do(mutate(., RMS.normal = RMS/filter(., Timepoint == "Baseline")$RMS,
            W1.amp.normal = W1.amp/filter(., Timepoint == "Baseline")$W1.amp,
            W1.lat.normal = W1.lat/filter(., Timepoint == "Baseline")$W1.lat) %>% print)

# Attempt to transform to normal
AOV.data$RMS.Gaus = LambertW::Gaussianize(AOV.data$RMS)[, 1]
AOV.data$W1.amp.Gaus = LambertW::Gaussianize(AOV.data$W1.amp)[, 1]
AOV.data$W1.lat.Gaus = LambertW::Gaussianize(AOV.data$W1.lat)[, 1]


# TH ----------------------------------------------------------------------

# Calculate from lowest Intensity
Wave1_TH =
  Wave1_data %>%
    group_by(ID, Timepoint, Freq) %>%
    do(filter(., Inten == min(Inten))) %>%
    ungroup() %>%
    summarise(TH = mean(Inten, na.rm = TRUE),
              .by = c(ID, Timepoint, Freq))

Wave1_TH %>%
  group_by(Freq) %>%
  do(dplyr::filter(., TH == max(TH)))


# Graph -------------------------------------------------------------------

# RMS, W1 latency (ms), W1 amp (uV)
Wave1_data  %>%
  ggplot(aes(x = Timepoint, y = W1.lat, color = Inten, group = Inten)) +
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
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )


Wave1_data %>%
  filter(Inten == "80" & Freq %in% c("BBN")) %>%
  gather(key = "ABR", value = "value", RMS, W1.lat, W1.amp) %>%
  ggplot(aes(x = Timepoint, y = value, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1.2) +
  labs(title = "80dB Intensity Only for Broadband",
       color = "Treatment") +
  facet_wrap( ~ ABR, scale = "free", nrow = 3) +
  # theme_classic() +
  theme_ipsum_es() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

ggsave(filename = "ABR_80dB.jpg",
       path = ProjectFolder,
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


# Graph -------------------------------------------------------------------

# RMS, W1 latency (ms), W1 amp (uV)
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
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

Wave1_data %>%
  filter(Inten == "80") %>%
  gather(key = "ABR", value = "value", RMS, `W1 latency (ms)`, `W1 amp (uV)`)%>%
  ggplot(aes(x = Timepoint, y = value, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1.2) +
  labs(title = "80dB Intensity Only for all Frequencies") +
  facet_wrap( ~ ABR, scale = "free", nrow = 3) +
  # theme_classic() +
  theme_ipsum_es() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

ggsave(filename = "ABR_80dB.jpg",
       path = ProjectFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


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
    mutate(adj.p.value = p.adjust(p.value, "BH"),
           sig = gtools::stars.pval(adj.p.value))

W1.amp.postHoc <-
  FSA::dunnTest(W1.amp ~ interaction(Freq, Timepoint, Condition),
                data = AOV.data,
                method = "bonf")

# print(postHoc, dunn.test.results = TRUE)

W1.amp.postHoc <-
  W1.amp.postHoc$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         Freq1 = str_split_fixed(Comp1, '\\.', 3)[,1],
         TimePoint1 = str_split_fixed(Comp1, '\\.', 3)[,2],
         Condition1 = str_split_fixed(Comp1, '\\.', 3)[,3],
         Freq2 = str_split_fixed(Comp2, '\\.', 3)[,1],
         TimePoint2 = str_split_fixed(Comp2, '\\.', 3)[,2],
         Condition2 = str_split_fixed(Comp2, '\\.', 3)[,3]) %>%
  select(-Comparison, -Comp1, -Comp2)


W1.amp.postHoc %>%
  filter(! Sig %in% c(" ", ".")) %>%
  View

# W1 Graph -------------------------------------------------------------------

# RMS, W1 latency (ms), W1 amp (uV)
Wave1_data  %>%
  ggplot(aes(x = Timepoint, y = `W1 latency (ms)`, color = Inten, group = Inten)) +
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
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

Wave1_data %>%
  filter(Inten == "80") %>%
  gather(key = "ABR", value = "value", RMS, `W1 latency (ms)`, `W1 amp (uV)`)%>%
  ggplot(aes(x = Timepoint, y = value, group = Condition, color = Condition)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.1) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point", size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "line", linewidth = 1.2) +
  labs(title = "80dB Intensity Only for all Frequencies") +
  facet_wrap( ~ ABR, scale = "free", nrow = 3) +
  # theme_classic() +
  theme_ipsum_es() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
  )

ggsave(filename = "ABR_BBN_80dB.jpg",
       path = ProjectFolder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)
