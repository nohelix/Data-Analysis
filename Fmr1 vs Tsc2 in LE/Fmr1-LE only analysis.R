
# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}


# TH Averages -------------------------------------------------------------
Fmr1_TH_averages = TH_table %>%
  filter(line == "Fmr1-LE") %>%
  group_by(genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")


# TH analysis -------------------------------------------------------------
# looking only at BBN mixed vs. alone 
# note that for Tukey everything needs to be a factor and Duration (a number
# value) tries to make itself continuous

Fmr1.TH.aov.data = TH_table %>%
  filter(line == "Fmr1-LE") %>%
  filter(Duration == "50" & detail == "Alone") 

Fmr1.TH.aov.data$Gaus = LambertW::Gaussianize(Fmr1.TH.aov.data$TH)[, 1]

Fmr1.TH.aov = aov(TH ~ genotype * as.factor(Frequency), 
                  data = Fmr1.TH.aov.data)

Parametric_Check(Fmr1.TH.aov)

# Normal
summary(Fmr1.TH.aov)

Fmr1_TH_averages %>%
  filter(Frequency == 0) %>%
  # Combine Frequency and Duration to create a single key column
  mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
  # Remove redundant columns
  select(-all_of(c("Frequency", "Duration"))) %>%
  spread(key, TH)

# TH Duration --------------------------------------------------------------

Fmr1.TH.BBN.aov.data = TH_table %>%
  filter(line == "Fmr1-LE") %>%
  filter(detail != "Rotating" & Frequency == 0) 

Fmr1.TH.BBN.aov.data$Gaus = LambertW::Gaussianize(Fmr1.TH.BBN.aov.data$TH)[, 1]

Fmr1.TH.BBN.aov = aov(TH ~ detail * Duration * genotype,
                      data = Fmr1.TH.BBN.aov.data)

Parametric_Check(Fmr1.TH.BBN.aov)

# Normal
summary(Fmr1.TH.BBN.aov)

# Only primary effects so no post-Hoc needed

# Rxn analysis ------------------------------------------------------------
# Get data
Fmr1_Rxn_over_TH = filter(Rxn_table_over_TH, line == "Fmr1-LE")

Fmr1_Rxn_over_TH$Gaus = LambertW::Gaussianize(Fmr1_Rxn_over_TH$Rxn)[, 1]


## Overall Model
  Fmr1_Rxn_overall_model = aov(Gaus ~ Frequency * Intensity * genotype,
                               data = filter(Fmr1_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  # Normality testing
  Parametric_Check(Fmr1_Rxn_overall_model)
  
  # Non-normal even with transformation
  kruskal.test(Rxn ~ genotype,
               data = filter(Fmr1_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  
  kruskal.test(Rxn ~ Frequency,
               data = filter(Fmr1_Rxn_over_TH, detail == "Alone" & Duration == "50"))


## BBN Model
  Fmr1.Rxn.BBN.aov.data = Fmr1_Rxn_over_TH %>%
    filter(Frequency == 0)
  
  Fmr1.Rxn.BBN.aov = aov(Gaus ~ detail * Duration * genotype,
                         data = Fmr1.Rxn.BBN.aov.data)
  
  # Normality testing
  Parametric_Check(Fmr1.Rxn.BBN.aov)

  # Not even close to normal
  # summary(Fmr1.Rxn.BBN.aov)

  # Kruskal Testing - only 
  lapply(c("detail", "Duration", "genotype"), 
         function(x) kruskal.test(reformulate(x, "Rxn"), data = Fmr1.Rxn.BBN.aov.data)) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "BH"),
           sig = gtools::stars.pval(adj.p.value))
  
  # Only Genotype is significant
  
  # Fmr1.Rxn.BBN.aov.postHoc <- 
  #   FSA::dunnTest(Rxn ~ interaction(detail, Duration, genotype), 
  #                 data = Fmr1.Rxn.BBN.aov.data,
  #                 method = "bonf")
  # 
  # Fmr1.Rxn.BBN.aov.postHoc$res %>% 
  #   as_tibble() %>%
  #   select(-P.unadj) %>%
  #   mutate(Sig = gtools::stars.pval(P.adj),
  #          D1 = str_extract(.$Comparison, '[:alpha:]+?\\.') %>% str_remove("\\."),
  #          D2 = str_extract(.$Comparison, ' - [:alpha:]+?\\.') %>% str_remove(" - ") %>% str_remove("\\."),
  #          T1 = str_extract(.$Comparison, '\\.[:digit:]+?\\.') %>% str_remove_all("\\."),
  #          T2 = str_extract(.$Comparison, ' - [:alpha:]+?\\.[:digit:]+?\\.') %>% str_remove(" - [:alpha:]+?\\.") %>% str_remove("\\."),
  #          l1 = str_extract(.$Comparison, '\\.[:alpha:]+[:digit:]') %>% str_remove("\\."),
  #          l2 = str_extract(.$Comparison, ' - .+-LE') %>% str_remove(" - .+\\..+\\.") %>% str_remove("-LE") ,
  #          G1 = str_extract(.$Comparison, '-LE.+?(\\W)') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
  #          G2 = str_extract(.$Comparison, '-LE\\.(Het|KO|WT)$') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
  #          P.adj = round(P.adj, digits = 3)) %>%
  #   # select(all_of(c("G1", "G2", "l1", "l2", "Z", "P.adj", "Sig"))) %>%
  #   filter(l1 == l2, G1 != G2, D1 == D2, T1 == T2)
  # 
  # Fmr1.Rxn.BBN.alone.postHoc <- 
  #   FSA::dunnTest(Rxn ~ interaction(Duration, line, genotype), 
  #                 data = filter(Fmr1.Rxn.BBN.aov.data, detail == "Alone"),
  #                 method = "bonf")
  # 
  # Fmr1.Rxn.BBN.alone.postHoc$res %>%
  #   as_tibble() %>%
  #   select(-P.unadj) %>%
  #   mutate(Sig = gtools::stars.pval(P.adj),
  #          Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
  #          Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
  #          dur1 = str_split_fixed(Comp1, '\\.', 3)[,1],
  #          line1 = str_split_fixed(Comp1, '\\.', 3)[,2],
  #          geno1 = str_split_fixed(Comp1, '\\.', 3)[,3],
  #          dur2 = str_split_fixed(Comp2, '\\.', 3)[,1],
  #          line2 = str_split_fixed(Comp2, '\\.', 3)[,2],
  #          geno2 = str_split_fixed(Comp2, '\\.', 3)[,3]) %>%
  #   select(-Comparison, -Comp1, -Comp2) %>%
  #   filter(line1 == line2, dur1 == dur2)



# Graphs ------------------------------------------------------------------

  Rxn_table %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
    filter(line == "Fmr1-LE") %>%
    rename(Intensity = `Inten (dB)`) %>%
    mutate(Frequency = str_replace_all(`Freq (kHz)`, pattern = "0", replacement = "BBN")) %>%
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
    # scale_linetype_manual(values = c("Fmr1-LE" = "solid", "Tsc2-LE" = "longdash")) +
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
  
ggsave(filename = paste0("Fmr1_Rxn_", single_Frequency,".jpg"),
       plot = last_plot(),
       width = 5, height = 6, units = "in", dpi = 300)

