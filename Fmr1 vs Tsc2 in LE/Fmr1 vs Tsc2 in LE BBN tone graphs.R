
# Variables ---------------------------------------------------------------

drop_TP3 = TRUE

# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}

n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}

# Reaction times ----------------------------------------------------------

Rxn_Plotter <- function(df) {
  ggplot(data = df, aes(x = Intensity, y = Rxn)) +
    # geom_point(aes(group = ID, color = Genotype), alpha = 0.3)+
    stat_summary(aes(color = genotype, linetype = line, group = interaction(line, genotype)),
                 fun = mean,
                 fun.min = function(x) mean(x) - se(x),
                 fun.max = function(x) mean(x) + se(x),
                 geom = "errorbar", width = 1.5, position = position_dodge(1)) +
    stat_summary(aes(color = genotype, shape = line, group = interaction(line, genotype)),
                 fun = mean,
                 geom = "point", position = position_dodge(1), size = 3) +
    stat_summary(aes(color = genotype, linetype = line, group = interaction(line, genotype)), 
                 fun = mean, geom = "line", position = position_dodge(1)) +
    labs(x = "Intensity (dB)",
         y = "Reaction time (ms, mean +/- SE)",
         color = "Genotype",
         caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
    # scale_linetype_manual(values = c("Fmr1-LE" = "solid", "Tsc2-LE" = "longdash")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
    scale_x_continuous(breaks = seq(0, 90, by = 10)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    )
}

Rxn_plot_data = Rxn_table %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  rename(Intensity = `Inten (dB)`) %>%
  mutate(Frequency = str_replace_all(`Freq (kHz)`, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$"))

# Single Frequency Graph
single_Frequency = "BBN"

Rxn_plot_data %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == single_Frequency) %>%
  Rxn_Plotter() +
  facet_wrap( ~ detail, ncol = 3, scales = "fixed") +
  labs(title = single_Frequency) +
  theme(legend.position = c(0.88, 0.8),
        legend.background=element_blank())

ggsave(filename = paste0("Rxn_", single_Frequency,".jpg"),
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

# All Frequencies Graph
Rxn_plot_data %>%
  filter(Intensity > 5) %>%
  filter(detail == "Alone") %>%
  mutate(Frequency = factor(Frequency, levels = c("4", "8", "BBN", "16", "32")),
         Type = if_else(Frequency == "BBN", "BBN", "Tones")) %>%
  Rxn_Plotter() +
  facet_wrap( ~ Frequency, ncol = 3, scales = "free") +
  theme(legend.position = c(0.88, 0.22),
        legend.background=element_blank())

ggsave(filename = "Rxn_overall.jpg",
       plot = last_plot(),
       width = 11, height = 4.8, units = "in", dpi = 300)

# Thresholds --------------------------------------------------------------

# Boxplot
TH_table %>%
  filter(Frequency == 0 & Duration == 50) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(detail != "Rotating") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, fill = genotype, color = line, group = interaction(line, genotype))) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
  scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "lightcoral")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3"),
       fill = "Genotype") +
  facet_wrap( ~ detail, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


# Boxplot
TH_table %>%
  filter(Frequency == 0 & Duration == 50) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(detail != "Rotating") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = detail, y = TH, fill = genotype, color = line, group = detail)) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
  scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "lightcoral")) +
  labs(y = "Threshold (dB, mean +/- SE)",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3"),
       fill = "Genotype") +
  facet_wrap( ~ interaction(line, genotype), ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "TH.jpg",
       plot = last_plot(),
       width = 8, height = 6, units = "in", dpi = 300)

