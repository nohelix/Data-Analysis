
# Groupings ---------------------------------------------------------------

# TTS Pilot Group - 50ms data only b/c mostly mixed
# 2 rats "Green2" (16) & "Orange1" (21) have incomplete data set
Group_TTS_pilot = c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# 1 rat "Green12" (191) had threshold shift
Group_1 = c(190, 191, 192, 193)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# 1 rat "Blue2" (266) died during baseline
# 1 rat "Blue3" (267) never got to FA criterion and was dropped
Group_2 = c(265, 266, 267, 268)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# in progress
Group_3 = c(338, 339, 340, 341, 342, 343)

# Get Data ----------------------------------------------------------------

# Get pilot BBN threshold data
BBN_alone_TH_pilot = 
  TH_table_detail %>%
  # limit to individuals with baseline & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # remove group 3 (Blues) that have no temporal integration & in-progress group
  filter(rat_ID < 200) %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None")

# Get pilot BBN reaction time data
BBN_alone_Rxn_pilot =
  Rxn_table_detail %>%
  # limit to individuals with baseine & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # remove group 3 (Blues) that have no temporal integration & in-progress group
  filter(rat_ID < 200) %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None") 

# All data
BBN_TempInt_TH = 
  TH_table_detail %>%
  # limit to individuals with baseline & post-HL
  filter(! rat_name %in% c("Green2", "Orange1", "Blue2")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # Remove struggling rat
  filter(rat_name != "Green24") %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None") %>%
  # Add group numbers to look for effects between groups
  mutate(group = case_when(rat_ID %in% Group_TTS_pilot ~ "Pilot",
                           rat_ID %in% Group_1 ~ "Group 1",
                           rat_ID %in% Group_2 ~ "Group 2",
                           rat_ID %in% Group_3 ~ "Group 3",
                           .default = "Unknown") %>% ordered(levels = c("Pilot", "Group 1", "Group 2", "Group 3")))

BBN_TempInt_Rxn =
  Rxn_table_detail %>%
  # limit to individuals with baseine & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None")%>%
  # Add group numbers to look for effects between groups
  mutate(group = case_when(rat_ID %in% Group_TTS_pilot ~ "Pilot",
                           rat_ID %in% Group_1 ~ "Group 1",
                           rat_ID %in% Group_2 ~ "Group 2",
                           rat_ID %in% Group_3 ~ "Group 3",
                           .default = "Unknown") %>% ordered(levels = c("Pilot", "Group 1", "Group 2", "Group 3")))

# Calculate Change in Rxn -------------------------------------------------

BBN_TempInt_Rxn_change =
  BBN_TempInt_Rxn %>%
  group_by(rat_ID, rat_name, group, Frequency, Duration, Intensity) %>%
  do(Rxn = filter(., HL_state == "baseline")$Rxn - filter(., HL_state == "post-HL")$Rxn %>%
       ifelse(identical(., numeric(0)), NA_integer_, .)) %>% 
  mutate(Rxn = ifelse(identical(Rxn, numeric(0)), NA_integer_, Rxn))


# Export for Ben ----------------------------------------------------------
# 
# BBN_alone_TH_pilot %>% 
#     fwrite(paste0(save_folder, "BBN_alone_TH_pilot_", Sys.Date(),".csv"), row.names = FALSE)
# 
# BBN_alone_Rxn_pilot %>% 
#   fwrite(paste0(save_folder, "BBN_alone_Rxn_pilot_", Sys.Date(),".csv"), row.names = FALSE)


# Threshold Graph ---------------------------------------------------------

# Boxplot - BBN Alone vs. Mixed
BBN_alone_TH_pilot %>%
# BBN_TempInt_TH %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  ggplot(aes(x = as.factor(Duration), y = TH, fill = HL_state, group = interaction(Duration, HL_state))) +
    geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
    stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
    scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "legend") +
    labs(title = "BBN thresholds by durations",
         x = "Stimulus Duration (ms)",
         y = "Threshold (dB)",
         fill = "Condition",
         caption = "Dropped individual with permanent threshold shifts. We only have data for the 50ms alone for all rats."
         ) +
    # facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.85, 0.93),
      legend.background=element_blank(),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    )

# ggsave(filename = "TH_BBN_pilot.jpg",
#        path = save_folder,
#        plot = last_plot(),
#        width = 8, height = 6, units = "in", dpi = 300)

BBN_TempInt_TH %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  ggplot(aes(x = as.factor(Duration), y = TH, 
             color = HL_state, fill = group, 
             group = interaction(Duration, group, HL_state))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_color_manual(values = c("coral", "grey30"), guide = "legend") +
  labs(title = "BBN thresholds by durations",
       x = "Stimulus Duration (ms)",
       y = "Threshold (dB)"
  ) +
  # facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.position = c(0.85, 0.93),
    legend.background=element_blank(),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  )


# Reaction Graph ----------------------------------------------------------

BBN_alone_Rxn_pilot %>%
# BBN_TempInt_Rxn %>%
  #limit to individuals with all 3 durations
  # filter(rat_name %in% c("Orange11", "Orange12", "Green11")) %>%
  filter(Intensity < 90) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = HL_state, group = interaction(HL_state))) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
  geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
              # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha)
  ) +
  # geom_text(data = Rxn_annotations, aes(label = sig), size = 10, show.legend = FALSE) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "",
       title = "Reaction curves for BBN before and after Hearing Loss",
       caption = "Following hearing loss, there is a loss of temproal integration. We only have data for the 50ms alone for all rats."
       ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  ) +
  facet_wrap( ~ Duration, ncol = 3, scales = "fixed") +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

# ggsave(filename = glue("BBN_alone_Rxn_pilot.jpg"),
#        path = save_folder,
#        plot = last_plot(),
#        width = 8, height = 6, units = "in", dpi = 300)


BBN_TempInt_Rxn %>%
  filter(Intensity < 95 & Intensity > 19) %>%
  filter(HL_state == "baseline") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  ggplot(aes(x = Intensity, y = Rxn, 
             color = as.factor(Duration), shape = as.factor(Duration),
             group = interaction(HL_state, as.factor(Duration)))) +
    geom_line(aes(group = interaction(HL_state, as.factor(Duration), rat_ID))) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - FSA::se(x),
                 fun.max = function(x) mean(x) + FSA::se(x),
                 geom = "errorbar", width = 1, position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
    geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
    labs(x = "Intensity (dB)",
         y = "Reaction time (ms, mean +/- SE)",
         color = "Duration", shape = "Duration",
         title = "Reaction curves for BBN showing temporal integration at baseline"
         ) +
    scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
      panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    ) +
    facet_wrap( ~ group, ncol = 2, scales = "fixed") +
    theme(legend.position = c(0.9, 0.8),
          legend.background=element_blank())

BBN_TempInt_Rxn %>%
  filter(Intensity < 90 & Intensity > 19) %>%
  filter(HL_state == "post-HL") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  ggplot(aes(x = Intensity, y = Rxn, 
             color = as.factor(Duration), shape = as.factor(Duration),
             group = interaction(HL_state, as.factor(Duration)))) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
  geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Duration", shape = "Duration",
       title = "Reaction curves for BBN for temporal integration after Hearing Loss"
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  ) +
  facet_wrap( ~ group, ncol = 2, scales = "fixed") +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

# BBN_TempInt_Rxn_change %>%
#   filter(Intensity < 90) %>%
#   mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
#            factor(levels = c("4", "8", "16", "32", "BBN"))) %>%
#   ggplot(aes(x = Intensity, y = Rxn, 
#              color = as.factor(Duration), shape = as.factor(Duration),
#              group = interaction(as.factor(Duration)))) +
#   stat_summary(fun = mean,
#                fun.min = function(x) mean(x) - FSA::se(x),
#                fun.max = function(x) mean(x) + FSA::se(x),
#                geom = "errorbar", width = 1, position = position_dodge(1)) +
#   stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
#   geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
#   labs(x = "Intensity (dB)",
#        y = "Change in Reaction time (ms, mean +/- SE)",
#        color = "Duration", shape = "Duration",
#        title = "Change in reaction curves for BBN for temporal integration after Hearing Loss"
#   ) +
#   scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#     panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#   ) +
#   facet_wrap( ~ group, ncol = 2, scales = "fixed") +
#   theme(legend.position = c(0.9, 0.8),
#         legend.background=element_blank())

