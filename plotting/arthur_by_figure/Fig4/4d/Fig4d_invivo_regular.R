library(ggplot2)
library(dplyr)
library(BuenColors)

## Load normalized data
df_vivo <- read.csv("/home/llab/users/big_storage/drive_a/chowa/denovo-cart-reproducibility/plotting/data/BCMA-CAR-T_in-vivo_full_ivis_normalized.csv")

## Calculate mean and SEM per treatment per week
summary_df <- df_vivo %>%
  filter(!is.na(pct_growth)) %>%
  group_by(treatment, week) %>%
  summarise(
    mean_pct = mean(pct_growth),
    se_pct = sd(pct_growth) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_pct - se_pct,
    ymax = mean_pct + se_pct
  )

## Set treatment order
summary_df$treatment <- factor(summary_df$treatment,
                               levels = c("Vehicle", "UTD", "Abecma", "B11"))

## Color mapping
color_mapping <- c(
  "Vehicle" = "#3B5B7A",
  "UTD" = "#B0A18F",
  "Abecma" = "#FFB81C",
  "B11" = "#1B9E77"
)

## Shape mapping
shape_mapping <- c(
  "Vehicle" = 16,
  "UTD" = 16,
  "Abecma" = 16,
  "B11" = 16
)

## Plot with pseudo-log y-axis
p <- ggplot(summary_df, aes(x = week, y = mean_pct, color = treatment,
                            shape = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000),
                     labels = c("0", "10", "100", "1,000", "10,000")) +
  scale_x_continuous(breaks = 0:4, labels = paste0("Week ", 0:4)) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = c(0.095, 0.15),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank()) +
  labs(x = "Weeks post-treatment",
       y = "% Tumor growth\n(photon flux, normalized to Week 0)")

p

## Save
cowplot::ggsave2("../plots/DNCT_rev/BCMA-in-vivo_ivis_tumor_growth.pdf", p, dpi = 300, width = 3, height = 2)
