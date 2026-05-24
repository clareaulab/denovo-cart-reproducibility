library(BuenColors)
library(dplyr)
library(data.table)


## Load the mean and standard error values
mean_df <- fread("../data/BCMA_incucyte_mean_alldonors.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline", "donor")) %>% rename(mean_value = value)

se_df <- fread("../data/BCMA_incucyte_se_alldonors.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline", "donor")) %>% rename(se_value = value)

valsdf_combined <- left_join(mean_df, se_df, by = c("elapsed_hours", "cellline", "donor", "variable")) %>%
  mutate(ymin_se = mean_value - se_value, ymax_se = mean_value + se_value)

## Color scheme
color_mapping <- c(
  "NB" = "#3B5B7A",
  "Abecma" = "#FFB81C",
  "X561726_WT" = "#D91E18",
  "int_1525" = "#1B9E77",
  "nonint_0366" = "#88CCEE"
)

## Donor 1
p1 <- valsdf_combined %>%
  filter(cellline == "K562", donor == "donor1") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p1
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_K562_donor1.pdf",p1,dpi=300,width = 2,height=1.5)

p2 <- valsdf_combined %>%
  filter(cellline == "K562 BCMA", donor == "donor1") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p2
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_K562_BCMA_OE_donor1.pdf",p2,dpi=300,width = 2,height=1.5)

p3 <- valsdf_combined %>%
  filter(cellline == "Raji", donor == "donor1") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p3
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_Raji_donor1.pdf",p3,dpi=300,width = 2,height=1.5)

p4 <- valsdf_combined %>%
  filter(cellline == "RPMI-8226", donor == "donor1") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p4
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_RPMI_8226_donor1.pdf",p4,dpi=300,width = 2,height=1.5)

## Donor 2
p5 <- valsdf_combined %>%
  filter(cellline == "K562", donor == "donor2") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p5
cowplot::ggsave2("../plots/BCMA_Incucyte_K562_donor2_4h.pdf",p5,dpi=300,width = 2,height=1.5)

p6 <- valsdf_combined %>%
  filter(cellline == "K562 BCMA", donor == "donor2") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p6
cowplot::ggsave2("../plots/BCMA_Incucyte_K562_BCMA_OE_donor2_4h.pdf",p6,dpi=300,width = 2,height=1.5)

p7 <- valsdf_combined %>%
  filter(cellline == "Raji", donor == "donor2") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p7
cowplot::ggsave2("../plots/BCMA_Incucyte_Raji_donor2_4h.pdf",p7,dpi=300,width = 2,height=1.5)

p8 <- valsdf_combined %>%
  filter(cellline == "RPMI-8226", donor == "donor2") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p8
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_RPMI_8226_donor2_4h.pdf",p8,dpi=300,width = 2,height=1.5)

## Donor 3
p9 <- valsdf_combined %>%
  filter(cellline == "K562", donor == "donor3") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p9
cowplot::ggsave2("../plots/BCMA_Incucyte_RPMI_8226_donor2_4h.pdf",p9,dpi=300,width = 2,height=1.5)

p10 <- valsdf_combined %>%
  filter(cellline == "K562 BCMA", donor == "donor3") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p10
cowplot::ggsave2("../plots/BCMA_Incucyte_RPMI_8226_donor2_4h.pdf",p10,dpi=300,width = 2,height=1.5)

p11 <- valsdf_combined %>%
  filter(cellline == "Raji", donor == "donor3") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p11
cowplot::ggsave2("../plots/BCMA_Incucyte_RPMI_8226_donor2_4h.pdf",p11,dpi=300,width = 2,height=1.5)

p12 <- valsdf_combined %>%
  filter(cellline == "RPMI-8226", donor == "donor3") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.7) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p12
cowplot::ggsave2("../plots/DNCT_rev/BCMA_Incucyte_RPMI_8226_donor3_4h.pdf",p12,dpi=300,width = 2,height=1.5)

## Combined grid
combined <- plot_grid(
  p1, p2, p3, p4,
  p5, p6, p7, p8,
  p9, p10, p11, p12,
  nrow = 3, ncol = 4,
  labels = c("Donor 1", "", "", "", "Donor 2", "", "", "", "Donor3", "", "", ""),
  label_size = 8
)
cowplot::ggsave2("../plots/BCMA_Incucyte_alldonors.pdf", combined, dpi = 300, width = 8, height = 4.5)

#plot average of donors 1-3
valsdf_mean <- valsdf_combined %>%
  group_by(elapsed_hours, cellline, variable) %>%
  summarize(
    mean_value = mean(mean_value, na.rm = TRUE),
    se_value = sd(mean_value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ymin_se = mean_value - se_value,
    ymax_se = mean_value + se_value
  )

p_mean <- valsdf_mean %>%
  filter(cellline %in% c("K562", "K562 BCMA", "Raji", "RPMI-8226")) %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size = 1.7) +
  geom_line() +
  facet_wrap(~ cellline, nrow = 1) +
  scale_y_continuous() +
  pretty_plot(fontsize = 8) +
  L_border() +
  scale_color_manual(values = color_mapping) +
  labs(x = "Hours elapsed", y = "cell growth *(% intensity)")

p_mean

cowplot::ggsave2(
  "../plots/DNCT_rev/BCMA_Incucyte_mean_3donors_allcelllines.pdf",
  p_mean,
  dpi = 300,
  width = 8,
  height = 2
)
