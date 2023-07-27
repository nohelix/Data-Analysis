
# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}


# TH Averages -------------------------------------------------------------
TH_averages = TH_table %>%
  group_by(line, genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")


TH_averages %>%
  filter(Frequency == 0) %>%
  # Combine Frequency and Duration to create a single key column
  mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
  # Remove redundant columns
  select(-all_of(c("Frequency", "Duration"))) %>%
  spread(key, TH) 


# TH analysis -------------------------------------------------------------
# looking only at BBN mixed vs. alone 
# note that for Tukey everything needs to be a factor and Duration (a number
# value) tries to make itself continuous

TH.aov.data = TH_table %>%
  filter(detail != "Rotating") %>%
  filter(Frequency == 0)  %>%
  mutate(Duration = as.factor(Duration))

TH.aov.data$Gaus = LambertW::Gaussianize(TH.aov.data$TH)[, 1]

TH.aov = aov(TH ~ detail * Duration * line * genotype,
             data = TH.aov.data)

Parametric_Check(TH.aov)

summary(TH.aov)

# Note only primary effects are significant
TH_stats = broom::tidy(TukeyHSD(TH.aov)) %>% mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  filter(sig != " " & str_detect(term, pattern = ":", negate = TRUE)) 


# TH overall --------------------------------------------------------------

TH.all.aov.data = TH_table %>%
  filter(detail == "Alone" & Duration == 50 & rat_name != "TP5")

TH.all.aov.data$Gaus = LambertW::Gaussianize(TH.all.aov.data$TH)[, 1]

TH.all.aov = aov(TH ~ Frequency * line * genotype,
             data = TH.all.aov.data)

Parametric_Check(TH.all.aov)

summary(TH.all.aov)

broom::tidy(TukeyHSD(TH.all.aov)) %>% mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  filter(sig != " " & str_detect(term, pattern = ":", negate = TRUE)) 

# Rxn analysis ------------------------------------------------------------

Rxn.BBN.aov.data = Rxn_table_over_TH %>%
  filter(detail != "Rotating") %>%
  filter(Frequency == 0) 

Rxn.BBN.aov.data$Gaus = LambertW::Gaussianize(Rxn.BBN.aov.data$Rxn)[, 1]

Rxn.BBN.aov = aov(Gaus ~ detail * Duration * line * genotype,
              data = Rxn.BBN.aov.data)


# Normality testing
Parametric_Check(Rxn.BBN.aov)

# Not even close to normal
# summary(TH.aov)

# Kruskal Testing - only 
lapply(c("detail", "Duration", "line", "genotype"), 
       function(x) kruskal.test(reformulate(x, "Rxn"), data = Rxn.BBN.aov.data)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sig = gtools::stars.pval(adj.p.value))

Rxn.BBN.aov.postHoc <- 
  FSA::dunnTest(Rxn ~ interaction(detail, Duration, line, genotype), 
                data = Rxn.BBN.aov.data,
                method = "bonf")

Rxn.BBN.aov.postHoc$res %>% 
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         D1 = str_extract(.$Comparison, '[:alpha:]+?\\.') %>% str_remove("\\."),
         D2 = str_extract(.$Comparison, ' - [:alpha:]+?\\.') %>% str_remove(" - ") %>% str_remove("\\."),
         T1 = str_extract(.$Comparison, '\\.[:digit:]+?\\.') %>% str_remove_all("\\."),
         T2 = str_extract(.$Comparison, ' - [:alpha:]+?\\.[:digit:]+?\\.') %>% str_remove(" - [:alpha:]+?\\.") %>% str_remove("\\."),
         l1 = str_extract(.$Comparison, '\\.[:alpha:]+[:digit:]') %>% str_remove("\\."),
         l2 = str_extract(.$Comparison, ' - .+-LE') %>% str_remove(" - .+\\..+\\.") %>% str_remove("-LE") ,
         G1 = str_extract(.$Comparison, '-LE.+?(\\W)') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
         G2 = str_extract(.$Comparison, '-LE\\.(Het|KO|WT)$') %>% str_remove("-LE\\.") %>% str_remove("\\W"),
         P.adj = round(P.adj, digits = 3)) %>%
  # select(all_of(c("G1", "G2", "l1", "l2", "Z", "P.adj", "Sig"))) %>%
  filter(l1 == l2, G1 != G2, D1 == D2, T1 == T2)

Rxn.BBN.alone.postHoc <- 
  FSA::dunnTest(Rxn ~ interaction(Duration, line, genotype), 
                data = filter(Rxn.BBN.aov.data, detail == "Alone"),
                method = "bonf")

Rxn.BBN.alone.postHoc$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         dur1 = str_split_fixed(Comp1, '\\.', 3)[,1],
         line1 = str_split_fixed(Comp1, '\\.', 3)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 3)[,3],
         dur2 = str_split_fixed(Comp2, '\\.', 3)[,1],
         line2 = str_split_fixed(Comp2, '\\.', 3)[,2],
         geno2 = str_split_fixed(Comp2, '\\.', 3)[,3]) %>%
  select(-Comparison, -Comp1, -Comp2) %>%
  filter(line1 == line2, dur1 == dur2)



## Overall Model
Rxn_overall_model = aov(Rxn ~ Frequency * Intensity * genotype * line, 
                        data = filter(Rxn_table_over_TH))
# Normality testing
Parametric_Check(Rxn_overall_model)


kruskal.test(Rxn ~ genotype, 
             data = Rxn_table_over_TH %>% 
               filter(! str_detect(Intensity, pattern = "5$")) %>% filter(Intensity < 90))

Rxn_overall_model_postHoc <- 
  FSA::dunnTest(Rxn ~ genotype, 
                data = Rxn_table_over_TH %>% 
                  filter(! str_detect(Intensity, pattern = "5$")) %>% filter(Intensity < 90),
                method = "bonf")

Rxn_overall_model_postHoc$res %>% 
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj))
