source("Fmr1 vs Tsc2 in LE data.R")
library(ggpmisc) # for themes
library(gridExtra)

save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Oddball individual graphs/"

n_at_mean <- function(x){
  #there are 3 positions so every rat is there in triplicate
  return(data.frame(y = mean(tail(x, n = 1)), label = paste0("n = ", length(x)/3)))
}

# Get core data -----------------------------------------------------------
# load rat_archive
rat_archive = fread(glue::glue("{projects_folder}/rat_archive.csv"), select = c("Rat_ID", "Sex", "Genotype"))

oddball_core_columns = c("date", "rat_name", "rat_ID", "Sex", "Genotype", "line", "genotype",
                        "file_name", "experiment", "phase", "task", "detail", 
                        "stim_type", "analysis_type", "complete_block_count", 
                        "FA_detailed", "reaction", "FA_percent")

oddball_core_data = dataset %>% 
  # decode rats
  left_join(rat_archive, by = c("rat_ID" = "Rat_ID")) %>%
  # Only keep relevant rats
  filter(rat_ID %in% filter(rat_archive, str_detect(Genotype, pattern = "Tsc2_LE|Fmr1-LE"))$Rat_ID) %>%
  # make columns by line and genotype
  mutate(line = str_extract(Genotype, pattern = "^.*(?=(-|_)LE)"),
         genotype = str_extract(Genotype, pattern = "(?<=LE_).*$")) %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% unnest_wider(stats) %>%
  select(all_of(oddball_core_columns)) %>%
  # drop Oddball and Octave
  filter(experiment %in% c("Oddball")) %>%
  mutate(challenge = case_when(task == "BG test" ~ "Background (30dB)",
                               task == "Probe trials" ~ "Probe trials",
                               task == "Uneven odds, 5 most frequent" ~ "Middle Odds",
                               task == "Uneven odds, 6 most frequent" ~ "End Odds",
                               task == "Uneven odds" ~ "End Odds",
                               task == "Reset" ~ "Reset",
                               task == "Catch trials" ~ "Catch",
                               task == "Base case" & detail == "Round 2" ~ "Base Case 2",
                               TRUE ~ str_extract(analysis_type, pattern = "(?<=\\().*(?=\\))")))
  # modify to find day after catch trials
  # group_by(rat_ID) %>% 
  # do(arrange(., date) %>% mutate(challenge = case_when(
  #   dplyr::lag(task %in% c("Catch trials", "Probe trials")) & task != "Catch trials" ~ "Day After Catch",
  #   TRUE ~ challenge)))



# Reaction ----------------------------------------------------------------

oddball_reaction_by_frequency = 
  oddball_core_data %>%
    # Omit Training & Reset days
    dplyr::filter(! task %in% c("Training")) %>%
    unnest(reaction) %>%
    rename(position = `Inten (dB)`) %>%
    mutate(frequency = str_extract(file_name, pattern = "^[:digit:]+?(?=kHz)") %>% 
             factor(levels = c("4", "8", "16", "32")),) %>%
    group_by(rat_ID, rat_name, position, frequency, genotype, line, challenge) %>%
    summarise(reaction = mean(Rxn, na.rm = TRUE) * 1000,
              .groups = "drop") %>%
    group_by(rat_ID, rat_name, frequency, genotype, line, challenge) %>%
    do(mutate(., reaction_norm = reaction/filter(., position == min(position))$reaction)) %>%
    ungroup

oddball_reaction = 
  oddball_core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "BG test")) %>%
  unnest(reaction) %>%
  rename(position = `Inten (dB)`) %>%
  group_by(rat_ID, rat_name, position, genotype, line, challenge) %>%
  summarise(reaction = mean(Rxn, na.rm = TRUE) * 1000,
            .groups = "drop") %>%
  group_by(rat_ID, rat_name, genotype, line, challenge) %>%
  do(mutate(., reaction_norm = reaction/filter(., position == min(position))$reaction)) %>%
  ungroup


oddball_FA = 
  oddball_core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "BG test")) %>%
  unnest(FA_detailed) %>%
  mutate(frequency = str_extract(file_name, "[:digit:]+(?=kHz)")) %>%
  group_by(rat_ID, rat_name, position, genotype, line, challenge, frequency) %>%
  summarise(FA = mean(FA_percent_detailed, na.rm = TRUE)*100,
            .groups = "drop") %>%
  group_by(rat_ID, rat_name, genotype, line, challenge) %>%
  # do(mutate(., reaction_norm = reaction/filter(., position == min(position))$reaction)) %>%
  ungroup


oddball_reaction_n_table = oddball_reaction %>%
  filter(! challenge %in% c("Probe trials", "Day After Catch", "Background (30dB)")) %>%
  group_by(line, genotype, challenge) %>% 
  summarise(n = length(unique(rat_ID)), .groups = "drop") %>%
  mutate(text_color = case_when(
    genotype == "WT" ~ "black", 
    genotype == "Het" ~ "blue", 
    genotype == "KO" ~ "red",
    TRUE ~ "green")) %>%
  arrange(genotype, desc(line), desc(challenge)) %>%
  rename(Line = line, Genotype = genotype, Condition = challenge)


# Graph -------------------------------------------------------------------

# All rats individual graphs

# 139, 145, 153
individual_graphs =
  oddball_reaction_by_frequency %>% 
    filter(challenge %in% c("Standard", "Base Case 2", "Middle Odds", "End Odds", "Background", "Reset", "Catch")) %>%
    # filter(rat_name == "RP3") %>%
    group_by(rat_ID, rat_name) %>%
    do(oddball_single_rat_graph = 
      ggplot(data = .,
             aes(x = position, y = reaction_norm,
                 color = challenge, fill = frequency,
                 group = challenge)) +
      # geom_smooth(se = FALSE, linewidth = 2) +
      # mean for genotypes across all frequencies
      stat_summary(fun = mean, geom = "line", linewidth = 2) +
      # mean for each frequency by genotype
      # stat_summary(aes(linetype = frequency, group = interaction(frequency, challenge)),
      #              geom = "line", fun = mean) +
      # geom_point(aes(shape = frequency), size = 5) +
      stat_summary(aes(shape = frequency, group = interaction(frequency, challenge)),
                   fun = mean, geom = "point", size = 5, stroke = 2) +
      scale_shape_manual(values = c("4" = 21, "8" = 22, "16" = 23, "32" = 24)) +
      scale_color_manual(values = c("Standard" = "black", "Base Case 2" = "grey30",
                                    "Background" = "tan4", "Catch" = "darkgreen",
                                    "End Odds" = "violetred", "Middle Odds" = "royalblue",
                                    "Reset" = "goldenrod")) +
      scale_x_continuous(breaks = seq(2, 6, by = 1)) +
      labs(title = paste0(unique(.$rat_name), " (#", unique(.$rat_ID), ") ", 
                          unique(.$line), " ", unique(.$genotype)),
           x = "Position of different 'go tone'",
           y = "Reaction time",
           fill = "Frequency", shape = "Frequency", linetype = "Frequency",
           color = "Condition") +
      theme_ipsum_es()
    )

# print(filter(individual_graphs, rat_name == "")$oddball_single_rat_graph)
# # Save individual graphs
# apply(individual_graphs, 1,
#       function(df) ggsave(filename = glue("Oddball {df$rat_name}.jpg"), # name of file
#                           path = save_folder, # location where file will save
#                           plot = df$oddball_single_rat_graph,
#                           width = 6, height = 4, units = "in", dpi = 300))
# individual_graphs$oddball_single_rat_graph

condition_to_graph = c("Standard")

oddball_frequency_graph =
  ggplot(oddball_reaction_by_frequency %>%
           filter(challenge %in% condition_to_graph), 
         aes(x = position, y = reaction,
             color = interaction(line, genotype), fill = frequency,
             group = interaction(line, genotype))) +
    # geom_smooth(se = FALSE, linewidth = 2) +
  # mean for genotypes across all frequencies
    stat_summary(fun = mean, geom = "line", linewidth = 2) +
  # mean for each frequency by genotype
    stat_summary(aes(linetype = frequency, group = interaction(frequency, line, genotype)), 
                 geom = "line", fun = mean) +
    # geom_point(aes(shape = frequency), size = 5) +
    stat_summary(aes(shape = frequency, group = interaction(frequency, line, genotype)), 
                 fun = mean, geom = "point", size = 5) +
    scale_color_manual(values = c("Tsc2.WT" = "grey40", "Tsc2.Het" = "blue", 
                                  "Fmr1.WT" = "black", "Fmr1.KO" = "red")) +
    scale_shape_manual(values = c("4" = 21, "8" = 22, "16" = 23, "32" = 24)) +
    scale_x_continuous(breaks = seq(2, 6, by = 1)) +
    labs(x = "Position of different 'go tone'",
         y = "Reaction time",
         fill = "Frequency", shape = "Frequency", linetype = "Frequency",
         color = "Genotype") +
    # table on graph
    annotate(geom = "table", x = 5.9, y = 290, 
             label = list(filter(oddball_reaction_by_frequency, challenge %in% condition_to_graph) %>%
                            group_by(line, genotype, challenge) %>% 
                            summarise(n = length(unique(rat_ID)), .groups = "drop") %>%
                            select(line, genotype, challenge, n)
                          ),
             table.theme = ttheme_gtplain(
               padding = unit(c(1, 0.9), "char")
             ),
             vjust = -2, hjust = -0.25) +
    theme_ipsum_es()

# print(oddball_frequency_graph)

oddball_frequency_basecase =
  ggplot(oddball_reaction_by_frequency %>%
           filter(challenge %in% c("Standard", "Base Case 2")), 
         aes(x = position, y = reaction_norm,
             color = interaction(line, genotype), 
             shape = challenge, linetype = challenge,
             group = interaction(line, genotype, challenge))) +
  # geom_smooth(se = FALSE, linewidth = 2) +
  # mean for genotypes across all frequencies
  stat_summary(fun = mean, geom = "line", linewidth = 2) +
  # mean for each frequency by genotype
  stat_summary(geom = "line", fun = mean) +
  # geom_point(aes(shape = frequency), size = 5) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  scale_color_manual(values = c("Tsc2.WT" = "grey40", "Tsc2.Het" = "blue", 
                                "Fmr1.WT" = "black", "Fmr1.KO" = "red")) +
  # scale_shape_manual(values = c("4" = 21, "8" = 22, "16" = 23, "32" = 24)) +
  scale_x_continuous(breaks = seq(2, 6, by = 1)) +
  labs(x = "Position of different 'go tone'",
       y = "Reaction time", shape = "Challenge",
       color = "Genotype") +
  # table on graph
  annotate(geom = "table", x = 5.9, y = 1, 
           label = list(filter(oddball_reaction_by_frequency, challenge %in% condition_to_graph) %>%
                          group_by(line, genotype, challenge) %>% 
                          summarise(n = length(unique(rat_ID)), .groups = "drop") %>%
                          select(line, genotype, challenge, n)
           ),
           table.theme = ttheme_gtplain(
             padding = unit(c(0.9, 0.9), "char")
           ),
           vjust = -2, hjust = -0.25) +
  theme_ipsum_es()

print(oddball_frequency_basecase)


oddball_BG_graph_table = 
  oddball_reaction_n_table %>%
  filter(Condition %in% c("Standard", "Background")) %>%
  arrange(Line, Genotype, Condition)

oddball_graph =
  oddball_reaction %>%
  filter(challenge %in% c("Standard", "Background")) %>%
  ggplot(aes(x = position, y = reaction_norm, 
             color = genotype, fill = line, linetype = challenge,
             group = interaction(line, genotype, challenge))) +
    # geom_smooth(se = FALSE, linewidth = 2) +
    stat_summary(geom = "line", fun = mean, linewidth = 2) +
    stat_summary(aes(shape = line), geom = "point", fun = mean, size = 3, stroke = 3) +
    # geom_point(aes(shape = line), size = 3, stroke = 3) +
    scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
    scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
    scale_linetype_manual(values = c("Standard" = "solid", "Background" = "dotdash", "End Odds" = "dotted")) +
    scale_x_continuous(breaks = seq(2, 6, by = 1)) +
    labs(x = "Position of different 'go tone'",
         y = "Reaction time\nNormalized by rat",
         fill = "Line", shape = "Line",
         color = "Genotype",
         linetype = "Condition") +
    # table on graph
    annotation_custom(
      tableGrob(select(oddball_BG_graph_table, -text_color),
                theme = ttheme_gtplain(
                            core = list(fg_params = list(col = matrix(oddball_BG_graph_table$text_color, 
                                                                      nrow(oddball_BG_graph_table),
                                                                      ncol(oddball_BG_graph_table)-1))),
                            ), 
                rows = NULL), 
      xmax = 4.6, ymax = 0.95
      ) +
    theme_ipsum_es() +
    guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
    theme(legend.key.width = unit(1.5,"cm"))

print(oddball_graph)

Oddball_odds_graph_table = 
  oddball_reaction_n_table %>%
  filter(Condition %in% c("Standard", "End Odds", "Middle Odds")) %>%
  arrange(Line, Genotype, desc(Condition))

oddball_odds_graph =
  oddball_reaction %>%
  filter(challenge %in% c("Standard", "End Odds", "Middle Odds")) %>%
  ggplot(aes(x = position, y = reaction_norm, 
             color = genotype, fill = line, linetype = challenge,
             group = interaction(line, genotype, challenge))) +
  # geom_smooth(se = FALSE, linewidth = 2) +
  stat_summary(geom = "line", fun = mean, linewidth = 2) +
  stat_summary(aes(shape = line), geom = "point", fun = mean, size = 3, stroke = 3) +
  # geom_point(aes(shape = line), size = 3, stroke = 3) +
  scale_shape_manual(values = c("Tsc2" = 21, "Fmr1" = 24)) +
  scale_fill_manual(values = c("Tsc2" = "slategrey", "Fmr1" = "tan4")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
  scale_linetype_manual(values = c("Standard" = "solid", "Middle Odds" = "dotdash", "End Odds" = "dotted")) +
  scale_x_continuous(breaks = seq(2, 6, by = 1)) +
  labs(x = "Position of different 'go tone'",
       y = "Reaction time\nNormalized by rat",
       fill = "Line", shape = "Line",
       color = "Genotype",
       linetype = "Condition") +
  # table on graph
  annotation_custom(
    tableGrob(select(Oddball_odds_graph_table, -text_color),
              theme = ttheme_gtplain(
                core = list(fg_params = list(col = matrix(Oddball_odds_graph_table$text_color, 
                                                          nrow(Oddball_odds_graph_table),
                                                          ncol(Oddball_odds_graph_table)-1))),
              ), 
              rows = NULL), 
    xmax = 4.6, ymax = 0.95
  ) +
  theme_ipsum_es() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  theme(legend.key.width = unit(1.5,"cm"))

print(oddball_odds_graph)
