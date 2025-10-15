library(BuenColors)
library(dplyr)
library(data.table)

## Load the mean and standard error values
mean_df <- fread("../data/BCMA_incucyte_mean.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline")) %>% rename(mean_value = value)

se_df <- fread("../data/BCMA_incucyte_se.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline")) %>% rename(se_value = value)

valsdf_combined <- left_join(mean_df, se_df,by = c("elapsed_hours", "cellline", "variable")) %>%
  mutate(ymin_se = mean_value - se_value,ymax_se = mean_value + se_value)

## Manual color palette
color_mapping <- c(
  "NB" = jdb_palette("corona")[1],
  "Abecma" = jdb_palette("corona")[2],
  "X561726_WT" = jdb_palette("corona")[3],
  "int_1525" = jdb_palette("corona")[4],
  "nonint_0366" = jdb_palette("corona")[5]
)
## Alternative color scheme
color_mapping <- c(
  "NB" = "black",
  "Abecma" = "#FFB81C",
  "X561726_WT" = "#D91E18",
  "int_1525" = "#8B0000",
  "nonint_0366" = "#9966CC"
)

p1 <- valsdf_combined %>%
  filter(cellline == "K562") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p1
cowplot::ggsave2("../plots/BCMA_Incucyte_K562.pdf",p1,dpi=300,width = 4.5,height=1.6)

p2 <- valsdf_combined %>%
  filter(cellline == "K562 BCMA") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p2
cowplot::ggsave2("../plots/BCMA_Incucyte_K562_BCMA_OE.pdf",p2,dpi=300,width = 4.5,height=1.6)
cowplot::ggsave2("../plots/BCMA_Incucyte_K562_BCMA_OE_Small.pdf",p2,dpi=300,width = 2.5,height=1.6)


p3 <- valsdf_combined %>%
  filter(cellline == "Raji") %>%
  filter(variable != "nonint_0366") %>%
  filter(variable != "X561726_WT") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p3
cowplot::ggsave2("../plots/BCMA_Incucyte_Raji_best_only.pdf",p3,dpi=300,width = 2,height=2)

p3_all <- valsdf_combined %>%
  filter(cellline == "Raji") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p3_all
cowplot::ggsave2("../plots/BCMA_Incucyte_Raji_all.pdf",p3,dpi=300,width = 2,height=2) 

p4 <- valsdf_combined %>%
  filter(cellline == "RPMI-8226") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p4
cowplot::ggsave2("../plots/BCMA_Incucyte_RPMI_8226.pdf",p4,dpi=300,width = 4.5,height=1.6)
cowplot::ggsave2("../plots/BCMA_Incucyte_RPMI_8226_Small.pdf",p4,dpi=300,width = 2,height=1.5)


raji_best_de_novo = valsdf_combined %>% filter(elapsed_hours == 72, cellline=="Raji", variable == "int_1525")
raji_abecma = valsdf_combined %>% filter(elapsed_hours == 72, cellline=="Raji", variable == "Abecma")

raji_best_de_novo_se2 = (raji_best_de_novo$se_value)**2
raji_abecma_se2 = (raji_abecma$se_value)**2

raji_tstat = (raji_best_de_novo$mean_value-raji_abecma$mean_value)/(sqrt(raji_best_de_novo_se2+raji_abecma_se2))
raji_tstat
pt(raji_tstat,df=5000)
# combined_plot <- cowplot::plot_grid(
#   p1,p2,p3,p4
# )
# combined_plot
# 
# cowplot::ggsave2("../plots/BCMA_Incucyte.pdf",combined_plot,dpi=300,width = 4.5,height=1.6)
library(lme4)

valsdf_combined_subset <- valsdf_combined %>% 
  filter(cellline=="Raji") %>%
  filter(variable == "Abecma" | variable == "int_1525") %>%
  filter(elapsed_hours != 0)

valsdf_prepared <- valsdf_combined_subset %>%
  # 1. Set the reference level (e.g., 'NB' for Uninfected control)
  mutate(variable = as.factor(variable)) %>%
  mutate(variable = relevel(variable, ref = "Abecma")) %>%
  # 2. Scale time (optional but often improves LMM convergence)
  mutate(time_scaled = elapsed_hours / max(elapsed_hours)) %>%
  mutate(log_mean_value = log(mean_value))

wls_fit <- lm(
  log_mean_value ~ variable * time_scaled,
  data = valsdf_prepared,
  weights = elapsed_hours
)

summary(wls_fit)

anova(wls_fit)

# lmm_fit <- lmer(
#   log(mean_value) ~ variable * time_scaled + (1 + time_scaled | cellline),
#   data = valsdf_prepared,
#   control = lmerControl(optimizer = "bobyqa") # Use a robust optimizer
# )

