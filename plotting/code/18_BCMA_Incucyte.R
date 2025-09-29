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

p3 <- valsdf_combined %>%
  filter(cellline == "Raji") %>%
  ggplot(aes(x = elapsed_hours, y = mean_value, color = variable, group = variable)) +
  geom_errorbar(aes(ymin = ymin_se, ymax = ymax_se), width = 1, linewidth = 0.3) +
  geom_point() + geom_line() +
  scale_y_continuous()+ pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "none") +
  labs(x="Hours elapsed", y="% Intensity to 0 hr")
p3
cowplot::ggsave2("../plots/BCMA_Incucyte_Raji.pdf",p3,dpi=300,width = 4.5,height=1.6)

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

# combined_plot <- cowplot::plot_grid(
#   p1,p2,p3,p4
# )
# combined_plot
# 
# cowplot::ggsave2("../plots/BCMA_Incucyte.pdf",combined_plot,dpi=300,width = 4.5,height=1.6)


