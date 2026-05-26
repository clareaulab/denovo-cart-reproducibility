library(tidyverse)
library(cowplot)
library(BuenColors)
library(reshape2)

construct_colors <- c(
  "NB"     = "#3B5B7A",
  "Abecma" = "#FFB81C",
  "B5"     = "#D91E18",
  "B5.I0"  = "#1B9E77"
)
construct_levels <- c("NB", "Abecma", "B5", "B5.I0")

read_csv("../data/BCMA_Restim2_LAG3-PD1.csv") %>%
  rename(Construct = Sample) %>%
  reshape2::melt(id.vars = "Construct",
                 variable.name = "Stimulation", value.name = "pct_double_pos") %>%
  mutate(
    Stimulation = as.integer(gsub("Stim", "", Stimulation)),
    Construct   = factor(Construct, levels = construct_levels)
  ) -> lag3_long

lag3pd1_summary %>%
  mutate(
    x_pos = factor(interaction(Construct, Stimulation, sep = "_"),
                   levels = paste(rep(construct_levels, each = 4),
                                  rep(1:4, times = 4), sep = "_"))
  ) -> lag3pd1_plot

lag3_long %>%
  mutate(
    x_pos = factor(interaction(Construct, Stimulation, sep = "_"),
                   levels = paste(rep(construct_levels, each = 4),
                                  rep(1:4, times = 4), sep = "_"))
  ) -> lag3_long_plot

lag3_long %>%
  mutate(
    Construct = factor(Construct, levels = construct_levels),
    x_pos = factor(interaction(Construct, Stimulation, sep = "_"),
                   levels = paste(rep(construct_levels, each = 4),
                                  rep(1:4, times = 4), sep = "_"))
  ) -> lag3_long_plot

ggplot(lag3pd1_plot, aes(x = x_pos, y = mean, fill = Construct)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.25, linewidth = 0.35) +
  geom_point(data = lag3_long_plot, aes(x = x_pos, y = pct_double_pos),
             shape = 21, size = 1.2, color = "black", fill = "black",
             inherit.aes = FALSE) +
  scale_fill_manual(values = construct_colors) +
  scale_x_discrete(labels = rep(paste0("Stim ", 1:4), 4)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "", y = "% LAG3+ PD-1+ of CAR-T") +
  guides(fill = "none") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    text        = element_text(family = "Helvetica")
  ) -> p_lag3pd1

p_lag3pd1
ggsave2(p_lag3pd1, filename = "../plots/DNCT_rev/BCMA_Restim2_LAG3PD1.pdf", width = 3, height = 2)
