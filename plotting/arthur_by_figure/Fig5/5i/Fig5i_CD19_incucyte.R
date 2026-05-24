library(BuenColors)
library(dplyr)
library(data.table)

## Load the mean and standard error values
mean_df <- fread("../data/CD19_incucyte_mean_updated.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline")) %>% rename(mean_value = value)

se_df <- fread("../data/CD19_incucyte_se_updated.csv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("elapsed_hours", "cellline")) %>% rename(se_value = value)

valsdf_combined <- left_join(mean_df, se_df,by = c("elapsed_hours", "cellline", "variable")) %>%
  filter(variable != "l115.D2") %>% 
  mutate(ymin_se = mean_value - se_value,ymax_se = mean_value + se_value)

## Manual color palette
color_mapping <- c(
  "NB" = "black",
  "FMC63" = "orange",
  "l112" = "#A63603",
  "l115" = "#08519C"
)

p1 <- valsdf_combined %>%
  filter(cellline == "K562") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.6) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p1
cowplot::ggsave2("../plots/DNCT_rev/CD19_Incucyte_K562_donor3.pdf",p1,dpi=300,width = 1.6,height=1.8)

p2 <- valsdf_combined %>%
  filter(cellline == "K562_CD19") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(size=1.6) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p2
cowplot::ggsave2("../plots/DNCT_rev/CD19_Incucyte_K562_CD19_OE_donor3.pdf",p2,dpi=300,width = 1.6,height=1.8)
cowplot::ggsave2("../plots/CD19_Incucyte_K562_CD19_OE_Small.pdf",p2,dpi=300,width = 2,height=1.5)

p3 <- valsdf_combined %>%
  filter(cellline == "K562_CD19_CD81_KD") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point(1.6) + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="cell growth *(% intensity)")
p3
cowplot::ggsave2("../plots/DNCT_rev/CD19_Incucyte_K562_CD19_CD81KD.pdf",p3,dpi=300,width = 1.6,height=1.8)
cowplot::ggsave2("../plots/CD19_Incucyte_K562_CD19_CD81KD_Small_all.pdf",p3,dpi=300,width = 2,height=1.5)

## Combined grid
combined <- plot_grid(
  p1, p2, p3,
  nrow = 1, ncol = 3
)
combined

cowplot::ggsave2("../plots/CD19_Incucyte_donor2_stacked.pdf", combined, dpi = 300, width = 6.5, height = 1.5)
