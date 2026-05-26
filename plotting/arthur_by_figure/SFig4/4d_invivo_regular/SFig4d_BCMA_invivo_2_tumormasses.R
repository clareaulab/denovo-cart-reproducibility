library(BuenColors)
library(dplyr)
library(data.table)
library(tidyr)
library(ggbeeswarm)

raw_in <- read_csv("../data/DNCT_rev/EXPT002_invivo_TumorMasses.csv")

tumor_df <- raw_in %>% mutate(Treatment = factor(Treatment, levels = binder_order)) %>% group_by(Treatment)

binder_order <- c("Vehicle", "UTD", "Abecma", "B5.I0")

color_mapping <- c(
  "Vehicle" = "#3B5B7A",
  "UTD" = "#B0A18F",
  "Abecma" = "#FFB81C",  
  "B5.I0" = "#1B9E77"
)

p_tumor <- ggplot(tumor_df, aes(x = Treatment, y = Mass)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_point(aes(fill = Treatment), position = position_jitter(width = 0.12, height = 0, seed = 42),
             size = 3, shape = 21, color = "black", stroke = 0.2) +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
  labs(x = "", y = "Tumor mass (g)") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none")

p_tumor
ggsave2(p_tumor, filename = "../plots/DNCT_rev/invivo_tumormass.pdf", width = 3, height = 2)
