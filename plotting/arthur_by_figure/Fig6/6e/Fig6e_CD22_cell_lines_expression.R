library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(BuenColors)   
library(cowplot)

raw_in <- read_csv("../data/CD22expression_celllines.csv")
df_summary <- raw_in %>% 
  group_by(cell_line) %>%
  summarise(mean_MFI = mean(MFI), sem = sd(MFI) / sqrt(n())) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "RPMI8226", "K562_CD22", "Raji")))

print(df_summary)

p1 <- ggplot(df_summary, aes(x= cell_line, y = mean_MFI, fill = cell_line)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD")) +
  geom_errorbar(aes(ymin = mean_MFI - sem, ymax = mean_MFI + sem), 
                width = 0.2, linewidth = 0.3) +
  geom_jitter(
    data = raw_in %>%
      mutate(cell_line = factor(cell_line, levels = c("K562", "RPMI8226", "K562_CD22", "Raji"))),
    aes(x = cell_line, y = MFI),
    width = 0.12, size = 1, shape = 21, fill = "black",
    inherit.aes = FALSE
  ) +
  pretty_plot(fontsize = 8) +
  L_border() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0,0.07))) +
  labs(x = "cell line", y = "CD22 MFI")

p1
cowplot::ggsave2(p1, file = "../plots/DNCT_rev/CD22_cellline_expression.pdf", width = 2, height = 2)
