# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}



Rxn_table %>%
  filter(task %in% c("Rotating", "CNO 3mg/kg", "Saline")) %>%
  filter(rat_ID != 100) %>%
  ggplot(aes(x = Position, y = Rxn, color = task)) +
  # geom_point(aes(group = ID, color = Genotype), alpha = 0.3)+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 1, position = position_dodge(0.1)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.1), size = 3) +
  stat_summary(fun = mean, geom = "line")  +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )

ggsave(filename = "Rxn with CNO.jpg",
       plot = last_plot(),
       width = 7.1, height = 5.75, units = "in", dpi = 300)  
