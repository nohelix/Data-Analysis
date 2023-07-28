
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# Location of the R scripts
code_folder = "Y:/GitHub/Data Analysis/Fmr1 vs Tsc2 in LE"

# Explict location to save files to:
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Liz"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41

# Drop TP3?
drop_TP3 = TRUE

# Working directory -------------------------------------------------------
setwd(code_folder)


# Load Packages -----------------------------------------------------------

# data loading
library(data.table); library(glue)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(broom);

# analysis & visualization
library(psycho); library(ggplot2); library(FSA)

# Graphing Functions _-------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}

n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}


# Load Necessary Datasets -------------------------------------------------
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)
rat_archive = fread(glue("{projects_folder}/rat_archive.csv"), 
                    select = c("Rat_ID", "DOB", "Sex", "Genotype"))

# Remove bad data
dataset = run_archive %>% 
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe 
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE))

core_columns = c("date", "rat_name", "rat_ID", 
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", "block_size", "complete_block_count", 
                 "dprime", "reaction", "FA_percent")

core_data = dataset %>% 
  # make columns with assignment data accessible; expands the dataframe
  unnest_wider(assignment) %>% 
  # Only keep relevant Experiments
  filter(experiment %in% c("Fmr1-LE", "Tsc2-LE")) %>%
  # drop Oddball and Octave
  filter(! phase %in% c("Octave", "Tone-BBN")) %>%
  # Omit Training & Reset days
  filter(! task %in% c("Training", "Reset")) %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(stats) %>%
  # limit columns to those we plan to use
  select(all_of(core_columns)) %>%
  # record date of hearing loss
  left_join(select(rat_archive, Rat_ID, Genotype), by = c("rat_ID" = "Rat_ID"))

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
  # Omit days with a higher FA% than the cutoff, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get dprimes
  unnest(dprime) %>%
  # Sort for ordered
  arrange(rat_ID, rat_name, Freq, Dur, dB, Genotype) %>%
  #Prep for Calculate_TH function
  nest(data = c(dB, dprime), .by = c(rat_ID, rat_name, detail, Freq, Dur, Genotype)) %>% 
  mutate(TH = map_dbl(data, Calculate_TH)) %>%
  select(-data) %>%
  rename(Duration = Dur, Frequency = Freq) %>%
  # Create neat genotype and line from Genotype
  mutate(line = str_remove(Genotype, pattern = "_WT|_KO|_Het") %>% str_replace(pattern = "_", replacement = "-"),
         genotype = str_extract(Genotype, pattern = "WT|KO|Het") )

TH_table %>%
  # use only comparable data
  filter(Duration == 50 & detail == "Alone") %>%
  summarise(TH = mean(TH, na.rm = TRUE),
            .by = c(Frequency, Duration, genotype, line)) %>%
  spread(genotype, TH, drop = TRUE)
  
# Note that the highest average TH is 50.0

TH_stats_Tsc2 =
  TH_table %>%
    # use only comparable data
    filter(Duration == 50 & detail == "Alone" & Frequency != 0) %>%
    filter(line == "Tsc2-LE") %>%
    mutate(Frequency = as.factor(Frequency)) %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
    aov(TH ~ genotype * Frequency, data = .) %>%
    # # Normal
    # .$residuals %>% shapiro.test
    # # summary
    # summary()
    # Post-hoc
    TukeyHSD() %>%
    tidy() %>%
    filter(term == "genotype:Frequency") %>%
    mutate(Sig = gtools::stars.pval(adj.p.value),
           F1 = str_extract(.$contrast, '\\:[:digit:]+?\\-') %>% str_remove("\\:")%>% str_remove("\\-"),
           F2 = str_extract(.$contrast, '\\:[:digit:]+?$') %>% str_remove("\\:"),
           G1 = str_extract(.$contrast, '^.+?\\:') %>% str_remove("\\:"),
           G2 = str_extract(.$contrast, '-.+?\\:') %>% str_remove("-") %>% str_remove("\\:")
           ) %>%
    filter(F1 == F2 & G1 != G2)

TH_stats_Fmr1 =
  TH_table %>%
  # use only comparable data
  filter(Duration == 50 & detail == "Alone" & Frequency != 0) %>%
  filter(line == "Fmr1-LE") %>%
  mutate(Frequency = as.factor(Frequency)) %>%
  aov(TH ~ genotype * Frequency, data = .) %>%
  # # Normal
  # .$residuals %>% shapiro.test
  # # summary
  # summary()
  # Post-hoc
  TukeyHSD() %>%
  tidy() %>%
  filter(term == "genotype:Frequency") %>%
  mutate(Sig = gtools::stars.pval(adj.p.value),
         F1 = str_extract(.$contrast, '\\:[:digit:]+?\\-') %>% str_remove("\\:")%>% str_remove("\\-"),
         F2 = str_extract(.$contrast, '\\:[:digit:]+?$') %>% str_remove("\\:"),
         G1 = str_extract(.$contrast, '^.+?\\:') %>% str_remove("\\:"),
         G2 = str_extract(.$contrast, '-.+?\\:') %>% str_remove("-") %>% str_remove("\\:")
  ) %>%
  filter(F1 == F2 & G1 != G2)
  



# TH Graphs ---------------------------------------------------------------

graph_line = "Tsc2-LE"

TH_table %>%
  # use only comparable data
  filter(Duration == 50 & detail == "Alone") %>%
  filter(Frequency != 0) %>%
  # limit to one line
  filter(line ==  graph_line) %>%
  mutate(Frequency = glue("{Frequency} kHz") %>% 
           factor(levels = c("4 kHz", "8 kHz", "16 kHz", "32 kHz"))) %>%
  # {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO", "Het"))) %>%
  ggplot(aes(x = genotype, 
             y = TH, fill = genotype, group = genotype)) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8, color = "grey40") +
  # stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
  scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       # caption = if_else(drop_TP3 & graph_line == "Tsc2-LE", "Without Het F TP3", ""),
       fill = "Genotype") +
  facet_wrap( ~ Frequency, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = "none"
  )

# Save last graph
ggsave(filename = glue("TH Plot {graph_line}.jpg"), # name of file
       path = save_folder, # location where file will save
       plot = last_plot(),
       width = 6, height = 3, units = "in", dpi = 300) # dimensions of saved file

# Get Reaction times by rat -----------------------------------------------
# NOTE: this is slow

Rxn_table = core_data %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # use only comparable data
  filter(`Dur (ms)` == 50 & detail == "Alone") %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`, Genotype) %>%
  # Get Averages
  summarise(Rxn = mean(Rxn, na.rm = TRUE) * 1000, .groups = "drop") %>% 
  # Rename the columns for easier calling
  rename(Intensity = `Inten (dB)`, Duration = `Dur (ms)`, Frequency = `Freq (kHz)`) %>%
  # Create neat genotype and line from Genotype
  mutate(line = str_remove(Genotype, pattern = "_WT|_KO|_Het") %>% str_replace(pattern = "_", replacement = "-"),
         genotype = str_extract(Genotype, pattern = "WT|KO|Het") )

Rxn_table %>%
  filter(Intensity > 50) %>%
  filter(Frequency != 0) %>%
  # Remove the 5's as they tend to make the line ugly
  # filter(! str_detect(Intensity, pattern = "5$")) %>%
  mutate(Frequency = glue("{Frequency} kHz") %>% 
           factor(levels = c("4 kHz", "8 kHz", "16 kHz", "32 kHz")),
         genotype = factor(genotype, levels = c("KO", "WT", "Het"))) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  summarise(Rxn = mean(Rxn, na.rm = TRUE),
            .by = c(Frequency, Duration, genotype, line)) %>%
  spread(genotype, Rxn, drop = TRUE)

Rxn.stats =
  Rxn_table %>%
    filter(Intensity >= 45) %>%
    filter(Frequency != 0) %>%
    mutate(Frequency = glue("{Frequency} kHz") %>% 
             factor(levels = c("4 kHz", "8 kHz", "16 kHz", "32 kHz")),
           genotype = factor(genotype, levels = c("KO", "WT", "Het"))) %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
    # aov(Rxn ~ genotype * Frequency, data = .) %>%
    # Non-normal
    # .$residuals %>% shapiro.test
    FSA::dunnTest(Rxn ~ interaction(Frequency, line, genotype), 
                  data = .,
                  method = "bonf") %>%
    .$res %>% 
    as_tibble() %>%
    select(-P.unadj) %>%
    mutate(Sig = gtools::stars.pval(P.adj),
           F1 = str_extract(.$Comparison, '[:digit:]+? kHz\\.') %>% str_remove("\\."),
           F2 = str_extract(.$Comparison, ' - [:digit:]+? kHz\\.') %>% str_remove(" - ") %>% str_remove("\\."),
           L1 = str_extract(.$Comparison, '\\..+?\\.') %>% str_remove_all("\\."),
           L2 = str_extract(.$Comparison, ' - [:digit:]+? kHz\\..+?\\.') %>% str_remove(" - [:digit:]+? kHz\\.") %>% str_remove("\\."),
           G1 = str_extract(.$Comparison, '-LE.+?(\\W)') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
           G2 = str_extract(.$Comparison, '-LE\\.(Het|KO|WT)$') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
           P.adj = round(P.adj, digits = 3),
           x = if_else(L1 == "Fmr1-LE", 1.5, 3.5),
           Rxn = 545) %>%
    filter(F1 == F2, L1 == L2) %>%
    select(all_of(c("F1", "G1", "G2", "L1", "x", "Rxn", "Z", "P.adj", "Sig")))  %>%
    rename(genotype = "G1", line = "L1", Frequency = "F1") %>%
    mutate(Frequency = factor(Frequency, levels = c("4 kHz", "8 kHz", "16 kHz", "32 kHz"))) %>%
    arrange(Frequency)
  
print(Rxn.stats)

# Rxn Graphs --------------------------------------------------------------

graph_line = "Tsc2-LE"

# All Frequencies Graph
Rxn_table %>%
  filter(Intensity >= 30) %>%
  filter(Frequency != 0) %>%
  filter(line ==  graph_line) %>%
  # Remove the 5's as they tend to make the line ugly
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  mutate(Frequency = as.factor(Frequency)) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  ggplot(aes(x = Intensity, y = Rxn)) +
  stat_summary(aes(color = genotype, group = interaction(genotype)), 
               fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(aes(color = genotype, group = interaction(genotype)),
               fun = mean,
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(aes(color = genotype, group = interaction(genotype)), 
               fun = mean, geom = "line", position = position_dodge(1)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       caption = if_else(drop_TP3 & graph_line == "Tsc2-LE", "Without Het F TP3", ""),
       color = "Genotype") +
  # scale_linetype_manual(values = c("Fmr1-LE" = "solid", "Tsc2-LE" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  # facet_wrap(~ line, scales = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) +
  theme(legend.background = element_blank(),
        legend.position = c(0.9, 0.8))


# Save last graph
ggsave(filename = glue("Rxn graph {graph_line}.jpg"), # name of file
       path = save_folder, # location where file will save
       plot = last_plot(),
       width = 6, height = 6, units = "in", dpi = 300) # dimensions of saved file

# Average Rxn box plots
Rxn_table %>%
  filter(Intensity >= 50) %>%
  filter(Frequency != 0) %>%
  filter(line ==  graph_line) %>%
  mutate(Frequency = glue("{Frequency} kHz") %>% 
           factor(levels = c("4 kHz", "8 kHz", "16 kHz", "32 kHz")),
         genotype = factor(genotype, levels = c("KO", "WT", "Het"))) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  ggplot(aes(x = genotype, y = Rxn, 
             fill = genotype, group = genotype)) +
    geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8, color = "grey40") +
    # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
    scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
    # You can add annotations:
    geom_text(data = filter(Rxn.stats, line ==  graph_line), 
              aes(x = 1.5, label = Sig), size = 12, show.legend = FALSE) +
    ylim(150, 600) +
    labs(x = "",
         y = "Reaction time (ms, mean +/- SE)",
         # caption = if_else(drop_TP3 & graph_line == "Tsc2-LE", "Without Het F TP3", ""),
         fill = "Genotype") +
    facet_wrap( ~ Frequency, ncol = 5) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
      legend.position = "none"
    )

# Save last graph
ggsave(filename = glue("Rxn plot {graph_line}.jpg"), # name of file
       path = save_folder, # location where file will save
       plot = last_plot(),
       width = 6, height = 3, units = "in", dpi = 300) # dimensions of saved file
