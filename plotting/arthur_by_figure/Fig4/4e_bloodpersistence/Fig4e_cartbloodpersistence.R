library(BuenColors)
library(dplyr)
library(data.table)

df <- fread("../data/car-t_mouse_blood_persistence.csv")
df <- df %>% mutate(group = factor(group, levels = c("Vehicle", "UTD", "C11D5.3", "B5.I0")))


color_mapping <- c(
  "Vehicle" = "#3B5B7A",
  "UTD" = "#B0A18F",
  "C11D5.3" = "#FFB81C",
  "B5.I0" = "#1B9E77"
)

ggplot(df, aes(x = group, y = counts, color = group)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 2.5) +
  scale_color_manual(values = color_mapping) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") +
  labs(x = "CAR Construct", y = "CD3+GFP+ events") -> p1
p1
cowplot::ggsave2(p1, file = "../plots/car-t_blood_persistence.pdf", dpi = 300, width = 2.5, height = 2.5)