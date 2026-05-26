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

df_mean <- read_csv("../data/BCMA_Restim2_incucyte_tumorgrowth_mean.csv") %>%
  melt(id.vars = c("elapsed_hours", "cellline"),
       variable.name = "Construct", value.name = "mean_value") %>%
  mutate(Construct = factor(Construct, levels = construct_levels))

df_se <- read_csv("../data/BCMA_Restim2_incucyte_tumorgrowth_se.csv") %>%
  melt(id.vars = c("elapsed_hours", "cellline"),
       variable.name = "Construct", value.name = "se_value") %>%
  mutate(Construct = factor(Construct, levels = construct_levels))

df_incucyte <- left_join(df_mean, df_se, by = c("elapsed_hours", "cellline", "Construct")) %>%
  mutate(ymin = mean_value - se_value, ymax = mean_value + se_value)

p_killing <- ggplot(df_incucyte, aes(x = elapsed_hours, y = mean_value,
                        color = Construct, group = Construct)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 1, linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_color_manual(values = construct_colors) +
  labs(x = "Hours elapsed", y = "% Intensity to 0 hr", color = "") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(text = element_text(family = "Helvetica")) +
  theme(legend.position = "none")

print(p_killing)
ggsave2(p_killing, filename = "../plots/DNCT_rev/BCMA_Restim2_Incucyte_killing.pdf", width = 2, height = 2)

