library(ggplot2)
library(dplyr)
library(BuenColors)
library(tidyr)

## Load normalized data
df_vivo <- read.csv("/home/llab/users/big_storage/drive_a/chowa/denovo-cart-reproducibility/plotting/data/DNCT_rev/BCMA_CAR-T_invivo_lowerET_fc.csv") %>% 
  drop_na(Week)

## Calculate mean and SEM per treatment per week
df_long <- df_vivo %>%
  pivot_longer(
    cols = c(Vehicle, UTD, Abecma, CARVYKTI, B11),
    names_to = "treatment",
    values_to = "value"
  )

summary_df <- df_long %>%
  group_by(treatment, Week) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sem = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

## Set treatment order
summary_df$treatment <- factor(summary_df$treatment,
                               levels = c("Vehicle", "UTD", "Abecma","CARVYKTI", "B11"))

## Color mapping
color_mapping <- c(
  "Vehicle" = "#3B5B7A",
  "UTD" = "#B0A18F",
  "Abecma" = "#FFB81C",
  "CARVYKTI" = "#A020F0",  
  "B11" = "#1B9E77"
)

p_lin <- ggplot(summary_df, aes(x = Week, y = mean_value, color = treatment,
                            group = treatment)) +
  geom_errorbar(aes(ymin = mean_value - sem, ymax = mean_value + sem),
                width = 0.15, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    breaks = c(0, 5000, 10000, 15000, 20000, 25000),
    labels = function(x) x / 100
  ) +
  scale_x_continuous(breaks = 0:5, labels = paste0("Week ", 0:5)) +
  scale_color_manual(values = color_mapping) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Weeks post-treatment",
       y = "% Tumor growth\n(photon flux, normalized to Week 0)")

p_lin

cowplot::ggsave2("../plots/DNCT_rev/BCMA_invivo_lowerET_lin.pdf", p_lin, dpi = 300, width = 3.1, height = 2.5)

p_log <- ggplot(summary_df, aes(x = Week, y = mean_value, color = treatment,
                                group = treatment)) +
  geom_errorbar(aes(ymin = mean_value - sem, ymax = mean_value + sem),
                width = 0.15, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_y_log10(
    breaks = c(0, 200, 2000, 20000),
    labels = c("0", "200", "2000", "20,000"),
  ) +
  scale_x_continuous(breaks = 0:5, labels = paste0("Week ", 0:5)) +
  scale_color_manual(values = color_mapping) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Weeks post-treatment",
       y = "% Tumor growth\n(photon flux, normalized to Week 0)")

p_log

cowplot::ggsave2("../plots/DNCT_rev/BCMA_invivo_lowerET_log.pdf", p_log, dpi = 300, width = 3.5, height = 2.5)

#statistics
library(emmeans)

df_endpoint <- df_long %>%
  filter(Week == 5, !is.na(value)) %>%
  mutate(log_value = log10(value + 1),
         treatment = factor(treatment,
                            levels = c("B11", "Vehicle", "UTD", "Abecma", "CARVYKTI")))
emm_endpoint <- emmeans(m_endpoint, ~ treatment)
endpoint_results <- contrast(emm_endpoint, method = "trt.vs.ctrl", ref = 1)
summary(endpoint_results)
