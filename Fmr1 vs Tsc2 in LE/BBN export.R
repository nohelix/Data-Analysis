left_join(select(BBN.dprime, -reaction), select(BBN.rxn, date, rat_ID, `Dur (ms)`, `Freq (kHz)`, `Inten (dB)`, Rxn), 
          by = c("date", "rat_ID", "Freq" = "Freq (kHz)", "dB" = "Inten (dB)", "Dur" = "Dur (ms)")) %>% 
  relocate(Rxn, .after = dprime) %>%
  fwrite(file = "Tsc2 and Fmr1 BBN data by day.csv")

rename(Rxn_table, Freq = `Freq (kHz)`, Dur = `Dur (ms)`, dB = `Inten (dB)`) %>%
  filter(Freq == 0) %>%
  full_join(core_data %>%
              # Omit Training & Reset days
              filter(! task %in% c("Training", "Reset")) %>%
              # Omit days with > 45% FA, i.e. guessing
              filter(FA_percent < FA_cutoff) %>%
              # Get dprimes
              unnest(dprime) %>%
              filter(Freq == 0) %>%
              group_by(rat_ID, detail, Freq, Dur, dB) %>%
              transmute(dprime = mean(dprime, na.rm = TRUE)) %>% unique(),
            by = c("rat_ID", "detail", "Freq", "dB", "Dur")) %>%
  relocate(dprime, .before = Rxn) %>%
  fwrite(file = "Tsc2 and Fmr1 BBN data by rat.csv")
  