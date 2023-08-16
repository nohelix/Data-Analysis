# Variables ---------------------------------------------------------------

drop_TP3 = TRUE

# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}

n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}

# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}


# TH Averages -------------------------------------------------------------
Tsc2_TH_averages = TH_table %>%
  filter(line == "Tsc2-LE") %>%
  group_by(genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")


# TH analysis -------------------------------------------------------------
# note that for Tukey everything needs to be a factor and Duration (a number
# value) tries to make itself continuous

Tsc2.TH.aov.data = TH_table %>%
  filter(line == "Tsc2-LE") %>%
  filter(Duration == "50" & detail == "Alone") %>%
  # Getting an NA from TP5 so drop
  filter(! (rat_name == "TP5" & Frequency == "4"))

Tsc2.TH.aov.data$Gaus = LambertW::Gaussianize(Tsc2.TH.aov.data$TH)[, 1]

Tsc2.TH.aov = aov(Gaus ~ genotype * as.factor(Frequency), 
                  data = Tsc2.TH.aov.data)

Parametric_Check(Tsc2.TH.aov)

# Non-Normal
# summary(Tsc2.TH.aov)

# Kruskal Testing - only 
lapply(c("Frequency", "genotype"), 
       function(x) kruskal.test(reformulate(x, "TH"), data = Tsc2.TH.aov.data)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sig = gtools::stars.pval(adj.p.value))

# No effect on thresholds

# TODO: BBN only with duration & detail

Tsc2_TH_averages %>%
  filter(Frequency == 0) %>%
  # Combine Frequency and Duration to create a single key column
  mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
  # Remove redundant columns
  select(-all_of(c("Frequency", "Duration"))) %>%
  spread(key, TH)

# TH Duration --------------------------------------------------------------

Tsc2.TH.BBN.aov.data = TH_table %>%
  filter(line == "Tsc2-LE") %>%
  filter(detail != "Rotating" & Frequency == 0) 

Tsc2.TH.BBN.aov.data$Gaus = LambertW::Gaussianize(Tsc2.TH.BBN.aov.data$TH)[, 1]

Tsc2.TH.BBN.aov = aov(TH ~ detail * Duration * genotype,
                      data = Tsc2.TH.BBN.aov.data)

Parametric_Check(Tsc2.TH.BBN.aov)

# Normal
summary(Tsc2.TH.BBN.aov)

# Only primary effects so no post-Hoc needed

## Graph
Tsc2.TH.BBN.aov.data %>%
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
  facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

# Rxn analysis ------------------------------------------------------------
# Get data
Tsc2_Rxn_over_TH = filter(Rxn_table_over_TH, line == "Tsc2-LE")

Tsc2_Rxn_over_TH$Gaus = LambertW::Gaussianize(Tsc2_Rxn_over_TH$Rxn)[, 1]


## Overall Model
  Tsc2_Rxn_overall_model = aov(Gaus ~ Frequency * Intensity * genotype,
                               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  # Normality testing
  Parametric_Check(Tsc2_Rxn_overall_model)
  
  # Non-normal even with transformation
  kruskal.test(Rxn ~ genotype,
               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  
  kruskal.test(Rxn ~ Frequency,
               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))


## BBN Model
  Tsc2.Rxn.BBN.aov.data = Tsc2_Rxn_over_TH %>%
    filter(Frequency == 0)
  
  Tsc2.Rxn.BBN.aov = aov(Gaus ~ detail * Duration * genotype,
                         data = Tsc2.Rxn.BBN.aov.data)
  
  # Normality testing
  Parametric_Check(Tsc2.Rxn.BBN.aov)
  
  # Not even close to normal
  # summary(Tsc2.Rxn.BBN.aov)
  
  # Kruskal Testing - only 
  lapply(c("detail", "Duration", "genotype"), 
         function(x) kruskal.test(reformulate(x, "Rxn"), data = Tsc2.Rxn.BBN.aov.data)) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "BH"),
           sig = gtools::stars.pval(adj.p.value))

  # Only Genotype is significant



# Graphs ------------------------------------------------------------------

Rxn_table %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(line == "Tsc2-LE") %>%
  # rename(Intensity = `Inten (dB)`) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == single_Frequency) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = genotype, group = genotype)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(1)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
  # scale_linetype_manual(values = c("Tsc2-LE" = "solid", "Tsc2-LE" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) +
  labs(title = single_Frequency,
       caption = "Combines all testing conditions (Alone, Mixed, & Rotating)") +
  theme(legend.position = c(0.88, 0.8),
        legend.background=element_blank())

ggsave(filename = paste0("Tsc2_Rxn_", single_Frequency,".jpg"),
       plot = last_plot(),
       width = 5, height = 6, units = "in", dpi = 300)

