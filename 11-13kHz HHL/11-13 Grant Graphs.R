# Package loading ---------------------------------------------------------

# data loading/manipulation
library(readxl); library(tidyverse); library(magrittr); library(dplyr);
library(tidyr); library(broom); library(glue)

# Data visualization
library(ggplot2); library(forcats); library(gtools); library(ggpmisc);
library(FSA) # for SE error bars

# Graph saving
library(svglite)

# Working directory -------------------------------------------------------
ProjectFolder = "C:/Users/Noelle/Box/ABR recordings/ABR Results/11-13kHz HHL/"


# Select Graphing Data ----------------------------------------------------

# Get data from excel spreadsheet
ABR_data = read_excel(glue("{ProjectFolder}Summary.xlsx"),
                      sheet = "4-32kHz", na = "NA")
ABR_data_detailed = read_excel(glue("{ProjectFolder}Summary.xlsx"),
                               sheet = "4-48kHz", na = "NA")

HHL_ABR_data =
  ABR_data_detailed %>%
  bind_rows(ABR_data) %>%
  # Get average for each rat across both ears
  summarise(RMS = mean(RMS, na.rm = TRUE),
            W1_amp = mean(`W1 amp (uV)`, na.rm = TRUE),
            W1_lat = mean(`W1 latency (ms)`, na.rm = TRUE),
            .by = c(ID, Condition, Freq, Inten)) %>%
  # Label frequency and order
  mutate(Freq = case_when(Freq != "0" ~ str_glue("{Freq} kHz"),
                          Freq == "0" ~ "BBN",
                          .default = "ERROR") %>%
           ordered(c("4 kHz", "6 kHz", "8 kHz", "12 kHz", "16 kHz",
                      "24 kHz", "32 kHz", "48 kHz", "BBN"))) %>%
  # fix Condition
  mutate(Condition = str_replace_all(Condition, "[45] week", "4-5 week"),
         Condition = ordered(Condition, c("Baseline", "Hearing Loss", "1 day",
                                          "1 week", "2 week", "4-5 week")))


HHL_ABR_data %>%
  summarise(n_ID = length(unique(ID)), ID = paste(sort(unique(ID)), collapse = ", "),
            .by = c(Condition, Freq)) %>%
  spread(Freq, n_ID) %>%
  relocate(ID, .after = BBN) %>%
  print


# RMS Grant Graph ---------------------------------------------------------------
# Signal-to-Noise ratio for BBN, significant for 70-90dB

RMS_graph =
  HHL_ABR_data  %>%
    filter(Freq %in% c("4 kHz", "8 kHz", "16 kHz", "32 kHz", "BBN")) %>%
    ggplot(aes(x = Inten, y = RMS, color = Condition, group = Condition)) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - se(x),
                 fun.max = function(x) mean(x) + se(x),
                 geom = "errorbar", width = 2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_x_continuous(limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) +
    labs(x = "Sound Intensity (dB)",
         y = "Signal-to-Noise Ratio (RMS)") +
    facet_wrap( ~ Freq, scales = "fixed", nrow = 3) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "white"),
      legend.position = c(0.8, 0.15)
    )
print(RMS_graph)

ggsave("11-13kHz_HHL_RMS.jpg",
       dpi = 300, width = 9, height = 8, units = "in",
       plot = RMS_graph, # or an explicit ggplot object name
       path = ProjectFolder)

ggsave("11-13kHz_HHL_RMS.svg",
       dpi = 300, width = 9, height = 8, units = "in",
       plot = RMS_graph, # or an explicit ggplot object name
       path = ProjectFolder)

RMS_graph_detailed =
  HHL_ABR_data  %>%
    filter(Condition %in% c("Baseline", "Hearing Loss", "4-5 week")) %>%
    ggplot(aes(x = Inten, y = RMS, color = Condition, group = Condition)) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - se(x),
                 fun.max = function(x) mean(x) + se(x),
                 geom = "errorbar", width = 3) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_x_continuous(limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) +
    labs(x = "Sound Intensity (dB)",
         y = "Signal-to-Noise Ratio (RMS)") +
    facet_wrap( ~ Freq, nrow = 2) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      panel.grid.major.x = element_line(color = "white"),
      legend.position = c(0.9, 0.2)
    )
print(RMS_graph_detailed)

ggsave("11-13kHz_HHL_RMS_4-48kHz.jpg",
       dpi = 300, width = 9, height = 5, units = "in",
       plot = RMS_graph_detailed, # or an explicit ggplot object name
       path = ProjectFolder)

# RMS ANOVA -------------------------------------------------------------------
AOV.data <- To_Graph %>%
            group_by(ID, Freq, Inten, Condition) %>%
            summarise(RMS = mean(RMS),
                      W1 = mean(`W1 amp (uV)`),
                      .groups = "drop")


AOV.data.detailed <- Graph_detailed %>%
              group_by(ID, Freq, Inten, Condition) %>%
              summarise(RMS = mean(RMS),
                        W1 = mean(`W1 amp (uV)`),
                        .groups = "drop")

RMS.aov <- aov(LambertW::Gaussianize(AOV.data$RMS) ~ Condition * Inten * Freq,
               data = AOV.data)

shapiro.test(RMS.aov$residuals)$p.value

summary(RMS.aov)

TukeyHSD(RMS.aov, "Condition", ordered = TRUE) %>% tidy() %>%
  select(-null.value, -conf.low, -conf.high) %>%
  mutate(adj.p.value = round(adj.p.value, digits = 4),
         sig = stars.pval(adj.p.value))

# W1 Grant Graph ---------------------------------------------------------------
# Wave 1 Amplitude for BBN to match Fmr1 KO result.

To_Graph  %>%
  ggplot(aes(x = Inten, y = `W1 amp (uV)`, color = Condition, group = Condition)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 3, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) +
  labs(x = "Sound Intensity (dB)",
       y = expression("Wave 1 Amplitued (\u00b5V)")) +
  facet_wrap( ~ Freq, nrow = 2) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
    legend.position = c(0.9, 0.2)
  )

ggsave("11-13kHz_HHL_W1lat.jpg",
       plot = last_plot(), # or an explicit ggplot object name
       path = ProjectFolder)

Graph_detailed %>%
  ggplot(aes(x = Inten, y = `W1 amp (uV)`, color = Condition, group = Condition)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 3, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) +
  labs(x = "Sound Intensity (dB)",
       y = expression("Wave 1 Amplitued (\u00b5V)")) +
  facet_wrap( ~ Freq, nrow = 2) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white"),
    legend.position = c(0.9, 0.2)
    )


ggsave("11-13kHz_HHL_W1lat_4-48kHz.jpg",
       plot = last_plot(), # or an explicit ggplot object name
       path = ProjectFolder)

# W1 Amp ANOVA -------------------------------------------------------------------
W1amp.aov <- aov(LambertW::Gaussianize(AOV.data$W1) ~ Condition * Inten * Freq,
                 data = AOV.data)

is_parametric = shapiro.test(W1amp.aov$residuals)$p.value > 0.05

if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
{writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(W1amp.aov$residuals)$p.value %>% round(digits = 3)))}

summary(W1amp.aov)

TukeyHSD(W1amp.aov, "Freq", ordered = TRUE) %>%
  tidy() %>%
  mutate(sig = stars.pval(adj.p.value))

TukeyHSD(W1amp.aov, "Condition", ordered = TRUE) %>% tidy() %>%
  select(-null.value, -conf.low, -conf.high) %>%
  mutate(adj.p.value = round(adj.p.value, digits = 4),
         sig = stars.pval(adj.p.value))


# 80dB Graph ---------------------------------------------------------------

Data = To_Graph
Freq = "80"

To_Graph_Table.1week <-
  Data %>%
  mutate(Freq = fct_relevel(Freq, c("4", "8", "16", "32", "BBN"))) %>%
  filter(Inten == !!Freq) %>%
  filter(Condition %in% c("Baseline", "1 week")) %>%
  group_by(Freq) %>%
  summarise(p = wilcox.test(`W1 amp (uV)` ~ Condition,
                            exact = FALSE,
                            alternative = "greater")$p.value %>% round(digits = 3),
            .groups = "drop") %>%
  mutate(BH = p.adjust(p, method = "BH") %>% round(digits = 3),
         Sig = stars.pval(`BH`)) %>%
  select(-p)

To_Graph_Table.2week <-
  Data %>%
  mutate(Freq = fct_relevel(Freq, c("4", "8", "16", "32", "BBN"))) %>%
  filter(Inten == !!Freq) %>%
  filter(Condition %in% c("Baseline", "2 week")) %>%
  group_by(Freq) %>%
  summarise(p = wilcox.test(`W1 amp (uV)` ~ Condition,
                            exact = FALSE,
                            alternative = "greater")$p.value %>% round(digits = 3),
            .groups = "drop") %>%
  mutate(BH = p.adjust(p, method = "BH") %>% round(digits = 3),
         Sig = stars.pval(`BH`)) %>%
  select(-p)

To_Graph_Table.5week <-
  Data %>%
  mutate(Freq = fct_relevel(Freq, c("4", "8", "16", "32", "BBN"))) %>%
  filter(Inten == !!Freq) %>%
  filter(Condition %in% c("Baseline", "4-5 week")) %>%
  group_by(Freq) %>%
  summarise(p = wilcox.test(`W1 amp (uV)` ~ Condition,
                            exact = FALSE,
                            alternative = "greater")$p.value %>% round(digits = 3),
            .groups = "drop") %>%
  mutate(BH = p.adjust(p, method = "BH") %>% round(digits = 3),
         Sig = stars.pval(`BH`)) %>%
  select(-p)

# wilcox.test(`W1 amp (uV)` ~ Condition,
#             exact = FALSE, correct = FALSE,
#             subset = Condition %in% c("Baseline", "1 week"),
#             data = filter(To_Graph, Inten == "80")
# )


Data  %>%
  filter(Inten == !!Freq) %>%
  mutate(Freq = fct_relevel(Freq, c("4", "8", "16", "32", "BBN"))) %>%
  ggplot(aes(x = Condition, y = `W1 amp (uV)`, color = Freq, group = Freq)) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - se(x),
                 fun.max = function(x) mean(x) + se(x),
                 geom = "errorbar", width = 0.1, position = position_dodge(0.03)) +
    stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.03)) +
    stat_summary(fun = mean, geom = "line", position = position_dodge(0.03)) +
    labs(x = "",
         y = expression("Wave 1 Amplitued (\u00b5V)"),
         color = "Frequency",
         title = "Wave 1 Amplitude at 80dB") +
    annotate(geom = 'table',
             x = 4.5, y = 0.0,
             label = list(To_Graph_Table.1week)) +
    annotate(geom = 'table',
             x = 5.5, y = 0.0,
             label = list(To_Graph_Table.2week)) +
    annotate(geom = 'table',
             x = 6.5, y = 0.0,
             label = list(To_Graph_Table.5week)) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      panel.grid.major.x = element_line(color = "white"),
      legend.position = c(0.8, 0.9)
    )

ggsave("11-13kHz_HHL_80dB.jpg",
       plot = last_plot(), # or an explicit ggplot object name
       path = ProjectFolder)

Graph_detailed  %>%
  # bind_rows(Graph) %>%
  filter(Inten == "80") %>%
  ggplot(aes(x = Condition, y = `W1 amp (uV)`, color = Freq, group = Freq)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 0.1, position = position_dodge(0.03)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.03)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.03)) +
  labs(x = "",
       y = expression("Wave 1 Amplitued (\u00b5V)"),
       color = "Frequency",
       title = "Wave 1 Amplitude at 80dB") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "white")
  )


# Thresholds --------------------------------------------------------------

ABR.TH_short <- read_excel(glue("{ProjectFolder}Summary.xlsx"),
                      sheet = "4-32_Thresholds", na = "NA")

ABR.TH_Detailed <- read_excel(glue("{ProjectFolder}Summary.xlsx"),
                              sheet = "4-48_Thresholds", na = "NA")

ABR.TH = ABR.TH_Detailed %>%
          bind_rows(ABR.TH_short)

ABR.TH %>%
  group_by(Condition, Ear, Freq) %>%
  summarise(n = length(unique(ID)), ID = paste(unique(ID), collapse = ", "),
            .groups = "drop") %>%
  print

ABR.TH %>%
  group_by(Condition, Freq) %>%
  summarise(n = length(unique(ID)),
            Run_TH = mean(Run_TH, na.rm = TRUE),
            RMS_TH = mean(RMS_TH, na.rm = TRUE),
            Wave_TH = mean(Wave_TH, na.rm = TRUE),
            .groups = "drop")

ABR.TH.aov = ABR.TH  %>%
  filter(Freq %in% c("4", "8", "16", "32", "BBN")) %>%
  # column 6 = RMS cuttoff, column 7 = Wave pick cutoff
  gather(key = "Measure", value = "TH", 7) %>%
  filter(!(is.na(TH)))

TH.aov <- aov(LambertW::Gaussianize(ABR.TH.aov$TH) ~ Condition * Freq,
                 data = ABR.TH.aov)

# TH.aov <- aov(TH ~ Condition * Freq,
#               data = ABR.TH.aov)

shapiro.test(TH.aov$residuals)$p.value

summary(TH.aov)

TukeyHSD(TH.aov, "Condition", ordered = TRUE) %>% tidy() %>%
  select(-null.value, -conf.low, -conf.high) %>%
  mutate(adj.p.value = round(adj.p.value, digits = 4),
         sig = stars.pval(adj.p.value))



# Save data ---------------------------------------------------------------

save(Graph, Graph_detailed, To_Graph, file = "11-13 ABR data.Rdata")
