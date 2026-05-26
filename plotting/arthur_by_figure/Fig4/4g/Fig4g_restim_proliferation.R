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

# 1. CAR-T Proliferation
prolif_summary<-read_csv("../data/BCMA_Restim2_CARTproliferation.csv") %>%
  reshape2::melt(id.vars = c("Construct", "Rep"),
                 variable.name = "Stimulation", value.name = "fold_expansion") %>%
  mutate(
    Stimulation = case_when(
      Stimulation == "Input" ~ 0L,
      startsWith(as.character(Stimulation), "Stim") ~ as.integer(gsub("Stim", "", Stimulation)),
      TRUE ~ NA_integer_
    ),
    Construct = factor(Construct, levels = construct_levels)
  ) %>%
  filter(!is.na(Stimulation)) %>%
  group_by(Construct, Stimulation) %>%
  summarise(
    mean = mean(fold_expansion),
    sem  = sd(fold_expansion) / sqrt(n()),
    .groups = "drop"
  )

ggplot(prolif_summary, aes(x = Stimulation, y = mean, color = Construct, group = Construct)) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15, linewidth = 0.35) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.8) +
  scale_color_manual(values = construct_colors) +
  scale_x_continuous(breaks = 0:4, labels = c("Input", "Stim 1", "Stim 2", "Stim 3", "Stim 4")) +
  labs(x = "", y = "CAR+ T cells\n(normalized to input)", color = "") +
  pretty_plot(fontsize = 8) + L_border() +
  theme(text = element_text(family = "Helvetica")) -> p_prolif

print(p_prolif)
ggsave2(p_prolif, filename = "../plots/DNCT_rev/BCMA_Restim2_CARTproliferation.pdf", width = 3, height = 2)