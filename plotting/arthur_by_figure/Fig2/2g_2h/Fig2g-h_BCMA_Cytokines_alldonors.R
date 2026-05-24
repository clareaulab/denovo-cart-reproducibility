library(tidyverse)
library(cowplot)
library(BuenColors)
library(reshape2)

construct_colors <- c(
  "UTD"    = "grey70",
  "NB"     = "grey50",
  "B2"     = "#4A90D9",
  "B5"     = "#E03030",
  "Abecma" = "#F5A623"
)
construct_levels <- c("UTD", "NB", "B2", "B5", "Abecma")
target_levels    <- c("CAR alone", "K562 mC","Raji", "K BCMA mC", "RPMI mC")
cytokine_levels  <- c("IFNG", "IL-2", "TNFa", "GranzymeB")

raw_in <- read_csv("../data/ELISA_BCMA_cytokines_alldonors.csv") %>%
  dplyr::select(cytokine, construct, target, rep1, rep2, rep3, donor) %>%
  dplyr::mutate(across(starts_with("rep"), ~ pmax(.x, 0)))

df_long_bcma_rpmi <- raw_in %>%
  filter(target == "RPMI mC") %>%
  pivot_longer(cols = starts_with("rep"), names_to = "rep", values_to = "value") %>%
  group_by(cytokine, construct, target, donor) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(construct = factor(construct, levels = construct_levels))

df_summary_bcma <- df_long_bcma_rpmi %>%
  group_by(cytokine, construct, target) %>%
  summarise(mean = mean(value),
            sem  = sd(value)/sqrt(n()),
            .groups = "drop") %>%
  mutate(construct = factor(construct, levels = construct_levels))

plot_cytokine_single <- function(cyt, ylab, show_strips = FALSE, show_xaxis = FALSE) {
  ggplot(
    filter(df_summary_bcma, cytokine == cyt),
    aes(x = construct, y = mean, fill = construct)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = 0.2, linewidth = 0.35) +
    geom_jitter(
      data        = filter(df_long_bcma, cytokine == cyt),
      aes(x = construct, y = value),
      position = position_jitter(width = 0.1, seed = 42),
      size = 0.8, shape = 21, fill = "black",
      inherit.aes = FALSE) +
    scale_fill_manual(values = construct_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "", y = ylab) +
    pretty_plot(fontsize = 7) + L_border() +
    theme(
      legend.position = "none",
      text            = element_text(family = "Helvetica"),
      axis.text.x     = if (show_xaxis) element_text(angle = 45, hjust = 1, size = 6) else element_blank(),
      axis.ticks.x    = if (show_xaxis) element_line() else element_blank(),
      plot.margin     = margin(2, 4, 2, 4))
}

p_bcma_ifng <- plot_cytokine_single("IFNG", "IFN-g (pg/mL)",         show_xaxis = FALSE)
p_bcma_il2  <- plot_cytokine_single("IL-2", "IL-2 (pg/mL)",          show_xaxis = FALSE)
p_bcma_tnfa <- plot_cytokine_single("TNFa", "TNF-a (pg/mL)",         show_xaxis = FALSE)
p_bcma_gzmb <- plot_cytokine_single("GranzymeB", "Granzyme B (pg/mL)",    show_xaxis = FALSE)

p_bcma_stack <- plot_grid(
  p_bcma_ifng, p_bcma_il2, p_bcma_tnfa, p_bcma_gzmb,
  nrow        = 1,
  align       = "h",
  axis        = "tb",
  rel_widths = c(1, 1, 1, 1)
)

print(p_bcma_stack)

ggsave2(p_bcma_stack,
        filename = "../plots/BCMA_RPMI-8226_allcytokines_stacked.pdf",
        width = 6, height = 1.2)

#compute statistics (might want to try other tests?)
library(rstatix)

ttest_table <- df_long_bcma %>%
  group_by(cytokine) %>%
  t_test(value ~ construct,
         comparisons = list(c("NB", "Abecma"),
                            c("B2", "Abecma"),
                            c("B5", "Abecma")),
         p.adjust.method = "BH") %>%
  dplyr::select(cytokine, group1, group2, n1, n2, statistic, df, p, p.adj, p.adj.signif)

print(ttest_table)

tukey_table <- df_long_bcma %>%
  group_by(cytokine) %>%
  tukey_hsd(value ~ construct) %>%
  filter((group1 == "NB"     & group2 == "Abecma") |
           (group1 == "Abecma" & group2 == "NB")     |
           (group1 == "B2"     & group2 == "Abecma") |
           (group1 == "Abecma" & group2 == "B2")     |
           (group1 == "B5"     & group2 == "Abecma") |
           (group1 == "Abecma" & group2 == "B5")) %>%
  dplyr::select(cytokine, group1, group2, estimate, conf.low, conf.high, p.adj, p.adj.signif)

print(tukey_table)

write_csv(ttest_table, "../tables/BCMA_RPMI-8226_ttests.csv")
write_csv(tukey_table, "../tables/BCMA_RPMI-8226_tukey.csv")

#heatmap for 2
p_heat <- raw_in %>%
  melt(id.vars = c("cytokine", "construct", "target")) %>%
  filter(target != "RPMI mC") %>%
  group_by(cytokine, construct, target) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  group_by(cytokine) %>%
  mutate(pct_max = mean / max(mean) * 100) %>%
  ungroup() %>%
  mutate(target    = factor(target,    levels = rev(target_levels))) %>%
  mutate(construct = factor(construct, levels = construct_levels)) %>%
  mutate(cytokine  = factor(cytokine,  levels = cytokine_levels)) %>%
  ggplot(., aes(x = construct, y = target, fill = pct_max)) +
  geom_tile(color = "black") +
  facet_wrap(~cytokine, ncol = 4) +
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "none") +
  labs(x = "construct", y = "co-culture line") +
  pretty_plot(fontsize = 8) +
  theme(
    panel.border       = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    legend.position    = "none",
    strip.background   = element_blank(),
    axis.text.x        = element_blank(),
  )

print(p_heat)
ggsave2(p_heat, filename = "../plots/BCMA_cytokines_heatmap_nolabels.pdf", width = 6.5, height = 1.2)
