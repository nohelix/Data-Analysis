source("Fmr1 vs Tsc2 in LE data.R")
library(ggpmisc) # for easy adding of tables to graph
library(glue)

# Get core data -----------------------------------------------------------
octave_core_columns = c("date", "rat_name", "rat_ID",
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", "complete_block_count", 
                 "FA_detailed", "reaction", "FA_percent", "hit_percent")

octave_core_data = dataset %>% 
  # remove invalid data
  filter(invalid != "TRUE") %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% unnest_wider(stats) %>%
  select(all_of(octave_core_columns)) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("Fmr1-LE", "Tsc2-LE")) %>%
  # drop Oddball and Octave
  filter(phase %in% c("Octave"))


# Decode rat data ---------------------------------------------------------
rat_archive = fread(glue::glue("{projects_folder}/rat_archive.csv"), select = c("Rat_ID", "Sex", "Genotype"))

octave_core_data = left_join(octave_core_data, rat_archive, by = c("rat_ID" = "Rat_ID")) %>%
  # make columns by line and genotype
  mutate(line = str_extract(Genotype, pattern = "^.*(?=(_|-)LE)"),
         genotype = str_extract(Genotype, pattern = "(?<=LE_).*$"))

# Discrimination Data -----------------------------------------------------

discrimination_data = 
  octave_core_data %>%
  filter(analysis_type == "Octave") %>%
  # get reaction time
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  # get d'
  unnest(FA_detailed) %>%
  # rename Frequency
  rename(Frequency = `Freq (kHz)`) %>%
  # calculate octave steps - string extract w/ regex is to get the go frequency
  mutate(octave_fraction = log(as.numeric(str_extract(file_name, pattern = "[:digit:]+?(?=-.+?kHz)"))/Frequency)/log(2),
         octave_steps = abs(round(octave_fraction * 12))) %>%
  # determine if 1/8 scale or zoomed in 1/16 scale
  mutate(Type = case_when(all((octave_steps %% 2) == 0) ~ "Broad",
                          any((octave_steps %% 2) != 0) ~ "Zoom",
                          .default = "Error"),
         Range = R.utils::seqToHumanReadable(octave_steps) %>% str_extract(pattern = "[:digit:]+-[:digit:]+"),
         .by = c(date, rat_ID))
  
  
# TH ----------------------------------------------------------------------

Calculate_TH_Octave <- function(df) {
  fit = loess(dprime ~ octave_steps, data = df)
  # plot(fit)
  TH = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")$y
  # print(TH)
  return(TH)
}

octave_TH_table =
  discrimination_data %>%
  # Sort for ordered (start w/ ID which is unique)
  arrange(rat_ID, rat_name, detail, octave_steps) %>%
  # Prep for Calculate_TH function
  nest(data = c(octave_steps, Frequency, dprime), .by = c(rat_ID, rat_name, Sex, Genotype, line, genotype)) %>%
  # mutate not summarise to keep the other columns
  mutate(TH = map_dbl(data, Calculate_TH_Octave)) %>%
  select(-data) 


# FA % --------------------------------------------------------------------

discrimination_FA_table =
  discrimination_data %>%
  group_by(rat_ID, rat_name, Sex, Genotype, line, genotype, detail, octave_steps) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

discrimination_FA_table_by_type =
  discrimination_data %>%
  group_by(rat_ID, rat_name, Sex, Genotype, line, genotype, detail, Type, octave_steps, Range) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

# Training stuff ----------------------------------------------------------

octave_training_data = 
  octave_core_data %>%
    filter(task == "Training") %>%
    summarise(days = n(), Genotype = unique(Genotype), 
              genotype = unique(genotype), line = unique(line), 
              .by = c(rat_ID, detail))

# octave_training_stats =
  aov(days ~ Genotype * detail, 
      data = octave_training_data) %>%
    # # Normal
    # .$residuals %>% shapiro.test()
    summary() %>%
    print
  
stats_Fmr1_learning = 
  t.test(days ~ Genotype, paired = FALSE, alternative = "two.sided",
         data = filter(octave_training_data, line == "Fmr1" & detail == "Normal"))
# stats_Fmr1_reverse = 
  # t.test(days ~ Genotype, paired = FALSE, alternative = "two.sided",
  #        data = filter(octave_training_data, line == "Fmr1" & detail == "Reversed"))
  
stats_Tsc2_learning = 
  t.test(days ~ Genotype, paired = FALSE, alternative = "two.sided",
         data = filter(octave_training_data, line == "Tsc2" & detail == "Normal"))
# stats_Tsc2_reverse = 
  # t.test(days ~ Genotype, paired = FALSE, alternative = "two.sided",
  #        data = filter(octave_training_data, line == "Tsc2" & detail == "Reversed"))

stats_learning =
  tidy(stats_Fmr1_learning) %>%
  mutate(line = "Fmr1", detail = "Normal", x = 1.5) %>%
  bind_rows(tidy(stats_Tsc2_learning) %>%
              mutate(line = "Tsc2", detail = "Normal", x = 3.5)) %>%
  rename(t = statistic, df = parameter) %>%
  mutate(Genotype = NA, genotype = NA,
         adj.p.value = p.adjust(p.value, n = 4, method = "bonferroni"),
         Sig = stars.pval(adj.p.value))

# N tables ----------------------------------------------------------------

discrimination_FA_n_table = 
  group_by(discrimination_FA_table, line, detail, genotype) %>% 
    summarise(n = length(unique(rat_ID)), .groups = "drop") %>% 
    transmute(Genotype = paste(line, genotype, sep = "."), detail = detail, n = n) %>%
  arrange(Genotype, detail)

discrimination_FA_n_table_by_Type = 
  group_by(discrimination_FA_table_by_type, line, detail, Type, genotype) %>% 
  summarise(n = length(unique(rat_ID)), .groups = "drop") %>% 
  transmute(Genotype = paste(line, genotype, sep = "."), 
            detail = detail, n = n, Type = Type) %>%
  arrange(Genotype, detail, Type)

octave_training_table =
  octave_core_data %>%
  filter(task == "Training") %>%
  group_by(line, detail, genotype) %>%
  summarise(n = length(unique(rat_ID)), .groups = "drop") %>% 
  transmute(Genotype = paste(line, genotype, sep = "."), 
            detail = detail, n = n) %>%
  arrange(Genotype, detail)
  

# Reversal ----------------------------------------------------------------

octave_summary_data =
  octave_core_data %>%
  summarise(date = "Average", day = 0, rat_name = unique(rat_name), 
            experiment = unique(experiment), genotype = unique(genotype), line = unique(line),
            phase = unique(phase), task = unique(task), detail = unique(detail),
            complete_block_count = mean(complete_block_count, na.rm = TRUE),
            FA_percent = mean(FA_percent, na.rm = TRUE), hit_percent = mean(hit_percent, na.rm = TRUE),
            .by = c(rat_ID, task, detail))

octave_holding_normal_short_data =
  octave_core_data %>%
  filter(task == "Holding" & detail == "Normal") %>%
  group_by(rat_ID) %>%
    #select only the last 3 days of holding prior to reversal
  do(arrange(., desc(date)) %>% head(n = 3) %>% mutate(day = row_number() - 4))

octave_reversal_data =
  octave_core_data %>%
  filter(task == "Training" & detail == "Reversed") %>%
    # set date as relative
  group_by(rat_ID) %>%
    do(arrange(., date) %>% mutate(day = row_number())) %>%
    # add last 3 days of Normal Holding
  rows_append(octave_holding_normal_short_data) %>%
    # add averages for Normal Holding
  mutate(date = as.character(date)) %>%
  rows_append(filter(octave_summary_data, task == "Holding" & detail == "Normal"))

# Graph -------------------------------------------------------------------

Octave_graph_all =
  ggplot(filter(discrimination_FA_table,
                ! detail %in% c("Reversed", "Reversed, Punished")), 
         aes(x = octave_steps, y = FA_percent_detailed * 100,
             color = genotype, fill = line, linetype = detail,
             group = interaction(detail, line, genotype))) +
    geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
    # mean for genotypes across all frequencies
    stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
    # mean for each frequency by genotype
    stat_summary(aes(shape = line), fun = mean, geom = "point",
                 position = position_dodge(.2), size = 2, stroke = 2) +
    # add labels by x-axis
    geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
                            tone = "No Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
              aes(label = tone), size = 5, show.legend = FALSE, vjust = 2, hjust = -0.2,
              family = "EconSansCndReg") +
    geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
                            tone = "Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
              aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
              family = "EconSansCndReg") +
    coord_cartesian(clip = "off") +
    # table on graph
    annotate(geom = "table", x = 12, y = 100,
             label = list(discrimination_FA_n_table),
             table.theme = ttheme_gtplain(
               padding = unit(c(1, 0.75), "char")
             )) +
    scale_x_continuous(breaks = seq(0, 12, by = 2)) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
    scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    labs(x = "Octave Step",
         y = "False Alarm %",
         title = "Discrimination across an octave",
         fill = "Line", shape = "Line",
         color = "Genotype") +
    theme_ipsum_es() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  theme(legend.key.width = unit(1.25,"cm"))

print(Octave_graph_all)

Octave_graph =
  ggplot(discrimination_FA_table %>%
           filter(detail %in% c("Normal", "Reversed", "Reversed, Punished")), 
         aes(x = octave_steps, y = FA_percent_detailed * 100,
             color = genotype, fill = line, linetype = detail,
             group = interaction(detail, line, genotype))) +
  geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
  # mean for genotypes across all frequencies
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
  # mean for each frequency by genotype
  stat_summary(aes(shape = line), fun = mean, geom = "point",
               position = position_dodge(.2), size = 2, stroke = 2) +
  # add labels by x-axis
  geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
                          tone = "No Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
            aes(label = tone), size = 5, show.legend = FALSE, vjust = 2, hjust = -0.2,
            family = "EconSansCndReg") +
  geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
                          tone = "Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
            aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
            family = "EconSansCndReg") +
  coord_cartesian(clip = "off") +
  # table on graph
  annotate(geom = "table", x = 12, y = 100,
           label = list(discrimination_FA_n_table %>%
                          filter(detail %in% c("Normal", "Reversed"))),
           table.theme = ttheme_gtplain(
             padding = unit(c(1, 0.75), "char")
           )) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
  scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_ipsum_es() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  theme(legend.key.width = unit(1.25,"cm"))

print(Octave_graph)

Octave_graph_by_Type =
  ggplot(discrimination_FA_table_by_type %>%
           filter(detail == "Normal"), 
         aes(x = octave_steps, y = FA_percent_detailed * 100,
             color = genotype, fill = line, linetype = detail,
             group = interaction(detail, line, genotype))) +
    geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
    # mean for genotypes across all frequencies
    stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
    # mean for each frequency by genotype
    stat_summary(aes(shape = line), fun = mean, geom = "point",
                 position = position_dodge(.2), size = 2, stroke = 2) +
    # add labels by x-axis
    # geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
    #                         tone = "No Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
    #           aes(label = tone), size = 5, show.legend = FALSE, vjust = 2, hjust = -0.2,
    #           family = "EconSansCndReg") +
    geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
                            tone = "Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
              aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
              family = "EconSansCndReg") +
    # # Graph lines for each individual
    # Note, alpha <1 not working (lines disappearing) for some reason
    # geom_line(aes(group = interaction(detail, line, genotype, as.factor(rat_ID)))) + 
    coord_cartesian(clip = "off") +
    scale_x_continuous(breaks = seq(0, 12, by = 2)) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
    scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    facet_wrap(~ Type, scales = "free_x", ncol = 2) +
    labs(x = "Octave Step",
         y = "False Alarm %",
         title = "Discrimination across an octave",
         fill = "Line", shape = "Line",
         color = "Genotype", linetype = "Detail") +
    theme_ipsum_es()

print(Octave_graph_by_Type)

discrimination_FA_n_table_by_Type %>%
  filter(detail == "Normal") %>%
  mutate(n = paste("n =", n)) %>%
  arrange(Type) %>%
  print

Octave_graph_by_Range =
  ggplot(discrimination_FA_table_by_type %>%
           filter(detail == "Normal" & Range != "5-7" & !is.na(Range)), 
         aes(x = octave_steps, y = FA_percent_detailed * 100,
             color = genotype, fill = Range, linetype = detail,
             group = interaction(Range, detail, line, genotype))) +
    geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
    # mean for genotypes across all frequencies
    stat_summary(aes(group = interaction(detail, line, genotype)),
                 fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
    # mean for each frequency by genotype
    stat_summary(aes(shape = line), fun = mean, geom = "point",
                 position = position_dodge(.2), size = 2, stroke = 2) +
    # add labels by x-axis
    # geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
    #                         tone = "No Go", genotype = "WT", line = "Fmr1", detail = "Normal"),
    #           aes(label = tone), size = 5, show.legend = FALSE, vjust = 2, hjust = -0.2,
    #           family = "EconSansCndReg") +
    geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
                            tone = "Go", genotype = "WT", line = "Fmr1", detail = "Normal", Range = NA),
              aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
              family = "EconSansCndReg") +
    # # Graph lines for each individual
    # Note, alpha <1 not working (lines disappearing) for some reason
    # geom_line(aes(group = interaction(detail, line, genotype, as.factor(rat_ID)))) + 
    coord_cartesian(clip = "off") +
    scale_x_continuous(breaks = seq(0, 12, by = 2)) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
    scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    labs(x = "Octave Step",
         y = "False Alarm %",
         title = "Discrimination across an octave by zoom grouping",
         fill = "Range", shape = "Line",
         color = "Genotype", linetype = "Detail") +
    theme_ipsum_es()

print(Octave_graph_by_Range)


# Learning graphs ---------------------------------------------------------

Octave_learning_plot =
  ggplot(data = octave_training_data,
         aes(x = fct_relevel(Genotype, c("Fmr1-LE_KO", "Fmr1-LE_WT", "Tsc2_LE_WT", "Tsc2_LE_Het")), 
             y = days, 
             fill = genotype, color = line, 
             group = Genotype)) +
    geom_boxplot(linewidth = 1) +
    # add individual points (layer 2)
    geom_point(aes(group = interaction(rat_ID, Genotype)), 
               shape = 1, show.legend = FALSE) +
    # Add p values labels (layer 3)
    geom_text(data = stats_learning, 
              aes(x = x, y = 78, label = glue("p = {round(adj.p.value, digits = 2)}")), 
              size = 5, show.legend = FALSE) +
    scale_fill_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    facet_wrap(~ detail) +
    labs(x = "Genotype",
         y = "# of days of trials to criterion",
         title = "Learning",
         fill = "Genotype", color = "Line") +
    theme_ipsum_es() +
    theme(panel.grid.major.x = element_line(colour="white", size=0.5))

print(Octave_learning_plot)

Octave_graph_Reversal_learning =
  ggplot(filter(octave_reversal_data),
         aes(x = day, 
             y = FA_percent * 100,
             # y = hit_percent * 100,
             color = genotype, fill = line,
             group = interaction(line, genotype))) +
    # Add criterion line
    geom_hline(aes(yintercept = 20), linewidth = 1.5, linetype = "dashed", color = "goldenrod") +
    # Individual lines
    geom_line(aes(group = interaction(line, genotype, rat_name)))+
    # mean for genotypes across all frequencies
    stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
    # mean for each frequency by genotype
    stat_summary(aes(shape = line), fun = mean, geom = "point",
                 position = position_dodge(.2), size = 2, stroke = 2) +
    # Mark the Day 0 is average
    geom_text(aes(x = 1.8, y = 15, label = "Average")) +
    # Add n table
    annotate(geom = "table", x = 60, y = 100,
             label = list(octave_training_table %>%
                            filter(detail %in% c("Reversed"))),
             table.theme = ttheme_gtplain(
               padding = unit(c(1, 0.75), "char")
             )) +
    scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
    scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    labs(x = "Days on Reversal",
         y = "False Alarm %",
         title = "Reversal",
         fill = "Line", shape = "Line",
         color = "Genotype") +
    theme_ipsum_es()

print(Octave_graph_Reversal_learning)

Octave_graph_Reversal_learning_hit =
  ggplot(filter(octave_reversal_data, day < 10),
         aes(x = day, 
             y = hit_percent * 100,
             color = genotype, fill = line,
             group = interaction(line, genotype))) +
  # Add criterion line
  geom_hline(aes(yintercept = 90), linewidth = 1.5, linetype = "dashed", color = "goldenrod") +
  # Individual lines
  geom_line(aes(group = interaction(line, genotype, rat_name)))+
  # mean for genotypes across all frequencies
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
  # mean for each frequency by genotype
  stat_summary(aes(shape = line), fun = mean, geom = "point",
               position = position_dodge(.2), size = 2, stroke = 2) +
  # Mark the Day 0 is average
  geom_text(aes(x = 0, y = 101, label = "Average")) +
  # Add n table
  annotate(geom = "table", x = 8, y = 40,
           label = list(octave_training_table %>%
                          filter(detail %in% c("Reversed"))),
           table.theme = ttheme_gtplain(
             padding = unit(c(1, 0.75), "char")
           )) +
  scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
  scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
  labs(x = "Days on Reversal",
       y = "False Alarm %",
       title = "Reversal",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_ipsum_es()

print(Octave_graph_Reversal_learning_hit)
