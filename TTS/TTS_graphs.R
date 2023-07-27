
# Variables ---------------------------------------------------------------

# Graphing ----------------------------------------------------------------
n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}

# Descriptive Stats -------------------------------------------------------
# Trial Count
stats_table %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  ggplot(aes(x = HL_state, y = trial_count, fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BA38", "#00BFC4")) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  labs(title = "Trial Count",
       x = "Hearing Loss Phase (through time)",
       y = "Trial Count",
       caption = parse(text = glue(
         "'No significant effect of hearing loss states ('*f[{filter(trial_count.aov.stats, term == 'HL_state')$df}]*' = {
         round(filter(trial_count.aov.stats, term == 'HL_state')$statistic, digits = 2)}, p = {
         round(filter(trial_count.aov.stats, term == 'HL_state')$p.value, digits = 2)}). Significantly fewer trials of tones ('*f[{filter(trial_count.aov.stats, term == 'stim_type')$df}]*' = {
         round(filter(trial_count.aov.stats, term == 'stim_type')$statistic, digits = 2)}, p = {
         round(filter(trial_count.aov.stats, term == 'stim_type')$p.value, digits = 2)})'"
       ))) +
  facet_wrap( ~ stim_type, ncol = 5, scales = "fixed", strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Trial Count.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# Trial Count by BG
stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL")),
         BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB"))) %>%
  ggplot(aes(x = HL_state, y = trial_count, fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  labs(title = "Trial Count",
       x = "Hearing Loss Phase",
       y = "Trial Count",
       caption = parse(text = glue(
         "'No significant effect of hearing loss states ('*f[{filter(trial_count_by_BG.aov.stats, term == 'HL_state')$df}]*' = {
         round(filter(trial_count_by_BG.aov.stats, term == 'HL_state')$statistic, digits = 2)}, p = {
         round(filter(trial_count_by_BG.aov.stats, term == 'HL_state')$p.value, digits = 2)}) or background ('*f[{filter(trial_count_by_BG.aov.stats, term == 'BG')$df}]*' = {
         round(filter(trial_count_by_BG.aov.stats, term == 'BG')$statistic, digits = 2)}, p = {
         round(filter(trial_count_by_BG.aov.stats, term == 'BG')$p.value, digits = 2)})'"
       ))) +
  facet_wrap( ~ BG, ncol = 5, scales = "fixed", strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Trial Count by BG.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# Hit % Count
# TODO:redo with new stats
stats_table %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  ggplot(aes(x = HL_state, y = (hit_percent * 100), fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BA38", "#00BFC4")) +
  labs(title = "Hit%",
       x = "Hearing Loss Phase (through time)",
       y = "Correct responses on go trials (Hits) %",
       caption = parse(text = glue(
         "'Significantly worse hit rate immediatedly after hearing loss ('*p[max]*' (BH adjusted) = {
         round(max(filter(hit_percent.aov.postHoc, term == 'HL_state' & str_detect(contrast, pattern = '-HL$'))$adj.p.value), digits = 2)}). Significantly fewer trials of tones ('*f[{filter(hit_percent.aov.stats, term == 'stim_type')$df}]*' = {
         round(filter(hit_percent.aov.stats, term == 'stim_type')$statistic, digits = 2)}, p = {
         round(filter(hit_percent.aov.stats, term == 'stim_type')$p.value, digits = 2)})'"
         ))) +
  facet_wrap( ~ stim_type, ncol = 5, scales = "fixed", strip.position = "bottom") +
  # # label only the tone facet
  # geom_text(data = data.frame(x = 2, y = 101, stim_type = "tone", HL_state = "HL",
  #                             label = unique(filter(hit_percent.aov.postHoc, term == 'HL_state' & str_detect(contrast, pattern = '-HL$'))$sig)),
  #           aes(x = x, y = y, label = label), size = 8) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Hit percent.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

label = tibble(
  HL_state = "baseline", FA_percent = FA_cutoff, stim_type = "BBN", label = glue("Guessing cutoff {FA_cutoff * 100}%")
)


# Hit rate by BG
stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL")),
         BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB"))) %>%
  ggplot(aes(x = HL_state, y = (hit_percent * 100), fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE,
               position = position_dodge(1), vjust = 2, size = 2) +
  # within BG
  geom_text(data = filter(hit_percent_by_BG.aov.postHoc, is.na(HL_state) & sig != " "),
            aes(x = 1.5, y = 95, label = sig), size = 8) +
  geom_segment(data = filter(hit_percent_by_BG.aov.postHoc, is.na(HL_state) & sig != " "),
               aes(x = 1, xend = 2, y = 94.8, yend = 94.8)) +
  # Between BG
  geom_text(data = filter(hit_percent_by_BG.aov.postHoc, !is.na(HL_state) &
                            sig != " " & str_detect(Comparison, pattern = "BBN", negate = TRUE)) %>%
              mutate(BG = str_remove(Comparison, pattern = " - Tones with no BG")),
            aes(y = 99, label = sig), size = 8) +
  geom_segment(data = tibble(BG = "Tones with no BG", HL_state = "baseline"),
               aes(x = 1, xend = Inf, y = 100, yend = 100)) +
  geom_segment(data = tibble(BG = "Pink noise at 30dB", HL_state = "baseline"),
               aes(x = -Inf, xend = Inf, y = 100, yend = 100)) +
  geom_segment(data = tibble(BG = "Pink noise at 50dB", HL_state = "baseline"),
               aes(x = -Inf, xend = 1, y = 100, yend = 100)) +
  labs(title = "Hit rate",
       x = "Hearing Loss Phase",
       y = "Hit rate",
      ) +
  facet_wrap( ~ factor(BG, levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                      "Pink noise at 50dB", "White noise at 50dB")),
              ncol = 5, scales = "fixed", strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Hit rate by BG.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# FA % Count
stats_table %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  ggplot(aes(x = HL_state, y = (FA_percent * 100), fill = HL_state, group = HL_state)) +
  geom_hline(yintercept = (FA_cutoff * 100), color = "blue", linetype = "longdash") +
  geom_text(data = label, aes(label = label), color = "blue", nudge_y = 2) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BA38", "#00BFC4")) +
  labs(title = "False Alarm %",
       x = "Hearing Loss Phase (through time)",
       y = "Incorrect responses on no-go trials (False Alarms) %",
       caption = parse(text = glue(
         "'Significantly more false alarms on tones ('*f[{filter(FA_percent.aov.stats, term == 'stim_type')$df}]*' = {
         round(filter(FA_percent.aov.stats, term == 'stim_type')$statistic, digits = 2)}, p = {
         round(filter(FA_percent.aov.stats, term == 'stim_type')$p.value, digits = 2)}).'"
       ))) +
  facet_wrap( ~ stim_type, ncol = 5, scales = "fixed", strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "FA percent.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# FA rate by BG
stats_table_by_BG %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL")),
         BG = case_when(
           BG_type == "None" & stim_type == "BBN" ~ "BBN with no BG",
           BG_type == "None" & stim_type == "tone" ~ "Tones with no BG",
           TRUE ~ paste0(BG_type, " noise at ", BG_Intensity, "dB")
         ) %>% factor(levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                 "Pink noise at 50dB", "White noise at 50dB"))) %>%
  ggplot(aes(x = HL_state, y = (FA_percent * 100), fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE,
               position = position_dodge(1), vjust = 2, size = 2) +
  # Between BG
  geom_text(data = filter(FA_percent_by_BG.aov.postHoc, sig != " " &
                            str_detect(contrast, pattern = "BBN|White", negate = TRUE)) %>%
              mutate(BG = str_remove(contrast, pattern = "Pink noise at 50dB-"),
                     HL_state = "baseline"),
            aes(y = 40, label = sig), size = 8, nudge_x = 0.5) +
  geom_segment(data = tibble(BG = "Tones with no BG", HL_state = "baseline"),
               aes(x = 1.5, xend = Inf, y = 40, yend = 40)) +
  geom_segment(data = tibble(BG = "Pink noise at 30dB", HL_state = "baseline"),
               aes(x = -Inf, xend = Inf, y = 40, yend = 40)) +
  geom_segment(data = tibble(BG = "Pink noise at 50dB", HL_state = "baseline"),
               aes(x = -Inf, xend = 1.5, y = 40, yend = 40)) +
  labs(title = "FA rate",
       x = "Hearing Loss Phase",
       y = "FA rate",
       caption = parse(text = glue(
         "'Significantly more false alarms with 50dB ('*f[{filter(FA_percent.aov.stats, term == 'stim_type')$df}]*' = {
         round(filter(FA_percent.aov.stats, term == 'stim_type')$statistic, digits = 2)}, p = {
         round(filter(FA_percent.aov.stats, term == 'stim_type')$p.value, digits = 2)}).'"
       ))) +
  facet_wrap( ~ factor(BG, levels = c("BBN with no BG", "Tones with no BG", "Pink noise at 30dB",
                                      "Pink noise at 50dB", "White noise at 50dB")),
              ncol = 5, scales = "fixed", strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "FA rate by BG.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# Thresholds --------------------------------------------------------------
# Boxplot - BBN Alone vs. Mixed
TH_table_detail %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
                      factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(BG == "None" & Frequency == "BBN") %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  ggplot(aes(x = as.factor(Duration), y = TH, color = detail, fill = HL_state, group = interaction(Duration, detail))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_text(aes(label = rat_name), alpha = 0.3, position = position_dodge(1)) +
  # stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "none") +
  scale_color_manual(values = c("Alone" = "black", "Mixed" = "azure3"), guide = "legend") +
  labs(title = "BBN thresholds for alone or mixed durations",
       x = "Stimulus Duration (ms)",
       y = "Threshold (dB)",
       color = "Stimulus\npresentation",
       caption = parse(text = glue("'No primary effect presentation style on threshold ('*chi^2*' = {
                                   round(BBN.detail.stats$statistic, digits = 2)}, p = {
                                   round(BBN.detail.stats$p.value, digits = 2)})'"))) +
  facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  )

ggsave(filename = "TH_BBN alone vs. mixed.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# Drop detail
TH_graph_data =
  TH_table %>%
  filter(Duration == 50) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(BG != "White\n50dB")


# Boxplot - Timeline of recovery
TH_graph_data %>%
  filter(BG == "None") %>%
  ggplot(aes(x = HL_state, y = TH, fill = HL_state, group = HL_state)) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_fill_manual(values = c("#F8766D", "#C77CFF", "#00BA38", "#00BFC4")) +
  labs(title = "Changes in threshold due to hearing loss",
       x = "Hearing Loss Phase",
       y = "Threshold (dB)",
       caption = "No background noise") +
  facet_wrap( ~ Frequency, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "TH_timeline.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# Boxplot
complete_HL_PKN = TH_graph_data %>% filter(BG == "Pink\n30dB" & HL_state == "post-HL") %>% .$rat_ID %>% unique()

TH_annotations =
  TH_postHoc_BG %>%
  mutate(BG = if_else(BG_Intensity == "None", "None", paste0("Pink\n", BG_Intensity, "dB"))) %>%
  filter(sig != " ") %>%
  left_join(TH_averages %>% filter(Duration == 50 & BG_type %in% c("None", "Pink")), by = c("Frequency", "BG_Intensity")) %>%
  group_by(Frequency, BG) %>%
  summarise(TH = max(TH), P.adj = round(unique(P.adj), digits = 2), .groups = "drop") %>%
  mutate(Frequency = as.character(Frequency),
         sig = stars.pval(P.adj))

TH_graph_data %>%
  filter(Frequency != "BBN") %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # filter(rat_ID %in% c(22, 15, 19, 20, 24, 17)) %>%
  ggplot(aes(x = Frequency, y = TH, fill = HL_state, group = interaction(Frequency, HL_state))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  geom_text(data = TH_annotations, aes(label = sig, fill = NULL, group = NULL), size = 10, nudge_y = 8) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  labs(title = "Changes in threshold in noise due to hearing loss",
       x = "Frequency (kHz)",
       y = "Threshold (dB)",
       fill = "Hearing Loss",
       caption = glue("Only 8kHz has significant difference
                      30dB: (Z = {round(filter(TH_postHoc_BG, sig != ' ' & BG_Intensity == 30)$Z, digits = 2)}, p = {
                      round(filter(TH_postHoc_BG, sig != ' ' & BG_Intensity == 30)$P.adj, digits = 2)}); 50dB: (Z = {
                      round(filter(TH_postHoc_BG, sig != ' ' & BG_Intensity == 50)$Z, digits = 2)}, p = {
                      round(filter(TH_postHoc_BG, sig != ' ' & BG_Intensity == 50)$P.adj, digits = 2)})")) +
  facet_wrap( ~ BG, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "TH_BG with HL.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# Boxplot - BBN Alone vs. step size
TH_table_steps %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(BG == "None" & Frequency == "BBN") %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  ggplot(aes(x = as.factor(Duration), y = TH, 
             color = as.factor(step_size), fill = HL_state, 
             group = interaction(Duration, step_size))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "none") +
  scale_color_manual(values = c("5" = "grey70", "10" = "grey30")) +
  labs(title = "BBN thresholds (single duration) by step size",
       x = "Stimulus Duration (ms)",
       y = "Threshold (dB)",
       color = "Stimulus\nStep Size") +
  facet_wrap( ~ HL_state, ncol = 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  )

# Boxplot - Baseline BBN Alone vs. step size
TH_table_steps %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(BG == "None" & Frequency == "BBN") %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(HL_state == "Baseline") %>%
  ggplot(aes(x = as.factor(step_size), y = TH, 
             color = as.factor(step_size), fill = HL_state, 
             group = interaction(step_size))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "none") +
  scale_color_manual(values = c("5" = "grey70", "10" = "grey30")) +
  labs(title = "BBN thresholds (single duration) by step size",
       x = "Stimulus Duration (ms)",
       y = "Threshold (dB)",
       color = "Stimulus\nStep Size") +
  facet_wrap( ~ Duration, ncol = 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  )


# Reaction times ----------------------------------------------------------

Rxn_Plotter <- function(plot) {
    plot +
    # geom_boxplot(aes(group = interaction(Intensity, detail))) +
    geom_smooth(se = FALSE, na.rm = TRUE) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - FSA::se(x),
                 fun.max = function(x) mean(x) + FSA::se(x),
                 geom = "errorbar", width = 1, position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
    # stat_summary(fun = mean, geom = "line", position = position_dodge(1)) +
    scale_color_manual(values = c("#F8766D", "#00BFC4")) +
    labs(x = "Intensity (dB)",
         y = "Reaction time (ms, mean +/- SE)",
         color = "") +
    scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    )
}

Rxn_Filter <- function(df){
  df %>%
  # filter(step_size = 10 & str_detect(Intensity, pattern = "0$")) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% factor(levels = c("4", "8", "BBN", "16", "32")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "White noise at 50dB", "Pink noise at 50dB"))) %>%
  filter(BG != "White noise at 50dB")
}

# BBN by Duration and HL_state
Rxn_annotations =
  detail.Rxn.postHoc %>%
  filter(sig != " ") %>%
  mutate(Intensity = 83,
         Rxn = 130,
         detail = "Mixed")

Rxn_table_over_TH_detail %>%
  Rxn_Filter() %>%
  filter(Frequency == "BBN") %>%
  #limit to individuals with both states
  filter(rat_name %in% c("Orange11", "Orange12", "Green11", "Green12")) %>%
  filter(rat_name != "Green12") %>%
  filter(Intensity < 90) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = HL_state, linetype = detail, shape = detail, group = interaction(detail, HL_state))) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
  geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
              method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha)
              ) +
  geom_text(data = Rxn_annotations, aes(label = sig), size = 10, show.legend = FALSE) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "") +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) +
  facet_wrap( ~ Duration, ncol = 3, scales = "fixed") +
  labs(subtitle = "Reaction curves for BBN before and after Hearing Loss",
       title = "Following hearing loss, there is a loss of temproal integration. We only have data for the 50ms alone for all rats.",
       caption = parse(text = glue("'Primary effect of presentation on reaction time (n = 4, '*chi^2*' = {round(BBN.detail.Rxn.stats$statistic, digits = 2)}, p = {round(BBN.detail.Rxn.stats$p.value, digits = 2)}). Significant only at baseline for 300ms. (Graph n = 3)'"))) +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

ggsave(filename = glue("Rxn_BBN detail and duration.jpg"),
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# BBN
Rxn_table_detail %>%
  Rxn_Filter() %>%
  filter(Frequency == "BBN") %>%
  # filter(rat_ID != 191) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = detail, group = detail)) %>%
  Rxn_Plotter() +
  facet_wrap( ~ HL_state, ncol = 3, scales = "fixed") +
  scale_color_manual(values = c("black", "grey70")) +
  labs(title = "Reaction curves for BBN before and after Hearing Loss",
       caption = parse(text = glue("'Non-significant difference at any Intesnity or Duration (300ms '*chi^2*' = {round(BBN.detail.Rxn.stats$statistic, digits = 2)}, p = {round(BBN.detail.Rxn.stats$p.value, digits = 2)})'"))) +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

ggsave(filename = glue("Rxn_BBN.jpg"),
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# BBN
Rxn_5step_table %>%
  Rxn_Filter() %>%
  filter(Frequency == "BBN") %>%
  # filter(rat_ID != 191) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = as.factor(step_size), group = as.factor(step_size))) %>%
  Rxn_Plotter() +
  facet_wrap( ~ HL_state, ncol = 3, scales = "fixed") +
  scale_color_manual(values = c("10" = "salmon", "5" = "green")) +
  labs(title = "Reaction curves for BBN before and after Hearing Loss",
       color = "Step Size") +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

# Single BG Graph
Single_BG_Grapher <- function(single_BG) {
  annotation = TH_averages %>% Rxn_Filter %>% filter(BG == single_BG) %>% filter(Duration == 50)

  graph =
    Rxn_table %>%
    Rxn_Filter() %>%
    filter(BG == single_BG) %>%
    ggplot(aes(x = Intensity, y = Rxn, color = HL_state, group = HL_state)) %>%
    Rxn_Plotter() +
    geom_vline(data = annotation, aes(xintercept = TH, color = HL_state, group = HL_state), linetype = "longdash", show.legend = FALSE) +
    facet_wrap( ~ Frequency, nrow = 2, scales = "free_x") +
    labs(title = glue("Reaction curves for {if_else(single_BG == 'None', 'no', single_BG)} background"),
         caption = parse(text = glue("'Non-significant difference at any background ('*chi^2*' = {round(filter(kruskal_results_Rxn, str_detect(model, pattern = 'HL_state'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_Rxn, str_detect(model, pattern = 'HL_state'))$p.adj, digits = 2)})'"))) +
    theme(legend.position = c(0.88, 0.22),
          legend.background=element_blank())

  print(graph)

ggsave(filename = glue("Rxn_BG_{single_BG}.jpg"),
       path = save_folder,
       plot = graph,
       width = 8, height = 6, units = "in", dpi = 300)

return(tibble_row(single_BG = single_BG, graph = graph, status = "Done"))
}

lapply(c("None", "Pink noise at 30dB", "Pink noise at 50dB"), Single_BG_Grapher) %>% bind_rows()

# All BGs for a single frequency Graph
Single_Frequency_Grapher <- function(single_frequency) {
  annotation = TH_averages %>% Rxn_Filter %>% filter(Frequency == single_frequency) %>% filter(Duration == 50)

  graph =
  Rxn_table %>%
    Rxn_Filter() %>%
    filter(Frequency == single_frequency) %>%
    ggplot(aes(x = Intensity, y = Rxn, color = HL_state, group = HL_state)) %>%
    Rxn_Plotter() +
    geom_vline(data = annotation, aes(xintercept = TH, color = HL_state, group = HL_state), linetype = "longdash", show.legend = FALSE) +
    facet_wrap( ~ BG, ncol = 3, scales = "free") +
    labs(title = glue("Reaction curves for {single_frequency}kHz"),
         caption = parse(text = glue("'Non-significant difference at any background ('*chi^2*' = {round(filter(kruskal_results_Rxn, str_detect(model, pattern = 'HL_state'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_Rxn, str_detect(model, pattern = 'HL_state'))$p.adj, digits = 2)})'"))) +
    theme(legend.position = "bottom",
          legend.background=element_blank())

  print(graph)

  ggsave(filename = glue("Rxn_{single_frequency}kHz.jpg"),
         path = save_folder,
         plot = graph,
         width = 11, height = 4.8, units = "in", dpi = 300)

return(tibble_row(single_frequency = single_frequency, graph = graph, status = "Done"))
}

lapply(c(4, 8, 16, 32), Single_Frequency_Grapher) %>% bind_rows()

# dprime curve ------------------------------------------------------------

# BBN Duration by Detail
dprime_detail_table %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
                              factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(Intensity <= 60) %>%
  ggplot(aes(x = Intensity, y = dprime, color = as.factor(Duration), linetype = detail, shape = detail, group = interaction(detail, Duration))) +
  geom_hline(yintercept = 1.5, linetype = "longdash", alpha = 0.5) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
              method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  ) +
  facet_wrap( ~ HL_state, ncol = 3, scales = "free_x") +
  scale_color_manual(values = c("magenta", "purple", "black")) +
  labs(x = "Intensity (dB)",
       y = "dprime (ms, mean +/- SE)",
       color = "Duration",
       linetype = "Simulus:\n\nPresentation",
       shape = "Simulus:\n\nPresentation",
       subtitle = "No significant effect of temporal integration (p > 0.19).",
       caption = parse(text = glue("'Primary effect of presentation style on ('*chi^2*' = {round(BBN.detail.dprime.stats$statistic, digits = 2)}, p = {round(BBN.detail.dprime.stats$p.value, digits = 15)}) and duration ('*chi^2*' = {round(BBN.duration.dprime.stats$statistic, digits = 2)}, p = {round(BBN.duration.dprime.stats$p.value, digits = 2)}), but no effect of hearing loss ('*chi^2*' = {round(BBN.HL.dprime.stats$statistic, digits = 2)}, p = {round(BBN.HL.dprime.stats$p.value, digits = 2)}).'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "dprime for BBN by Duration and Detail.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)


# drop white noise
# one plot for each Frequency, line color for each BG, line type for HL_state
dprime_graph_data =
  dprime_summary_table %>%
  # filter(str_detect(Intensity, pattern = "0$")) %>%
  filter(! (Frequency %in% c(4, 8) && Intensity == 0)) %>%
  filter(! (Frequency %in% c(32) && Intensity == -15)) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  filter(Duration == 50) %>% filter(BG_type != "White") %>%
  mutate(BG = if_else(BG_type == "None", "None", paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "White noise at 50dB", "Pink noise at 50dB"))) %>%
  filter(Frequency != 0)

dprime_graph_data %>%
  filter(Intensity <= 60) %>%
  filter(rat_name != "Green12") %>%
  ggplot(aes(x = Intensity, y = dprime, color = BG, linetype = HL_state, shape = HL_state, group = interaction(BG, HL_state))) +
    geom_hline(yintercept = 1.5, linetype = "longdash", alpha = 0.5) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - FSA::se(x),
                 fun.max = function(x) mean(x) + FSA::se(x),
                 geom = "errorbar", width = 1, position = position_dodge(0)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
    # stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
    geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
                # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
    ) +
    scale_color_manual(values = c("black", "salmon", "violetred")) +
    facet_wrap( ~ Frequency, ncol = 2, scales = "free_x") +
    labs(x = "Intensity (dB)",
         y = "dprime (ms, mean +/- SE)",
         color = "Background",
         linetype = "Hearing Loss",
         shape = "Hearing Loss",
         # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
         ) +
    scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    )

ggsave(filename = glue("dprime in BG all.jpg"),
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# more viewable
dprime_graph_data %>%
  filter(rat_name != "Green12") %>%
  filter(Intensity <= 80) %>%
  ggplot(aes(x = Intensity, y = dprime, color = BG, linetype = HL_state, shape = HL_state, group = interaction(BG, HL_state))) +
  geom_hline(yintercept = 1.5, linetype = "longdash", alpha = 0.5) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  # stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  geom_smooth(se = FALSE, na.rm = TRUE,
              # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
              # method = "nls", formula = y ~ SSquadp3xs(x, a, b, jp),
  ) +
  scale_color_manual(values = c("black", "salmon", "violetred")) +
  labs(x = "Intensity (dB)",
       y = "dprime (ms, mean +/- SE)",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)}), Non-significant effect of hearing loss ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'HL'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'HL'))$p.adj, digits = 2)})'"))) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    legend.position = c(0.88, 0.3),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = glue("dprime in BG all.jpg"),
       path = save_folder,
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)
