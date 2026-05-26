library(BuenColors)
library(dplyr)
library(ggplot2)
library(tidyr)

cd22_incucyte <- read.csv("../data/CD22_incucyte_mean_alldonors.csv")

cd22_incucyte_long <- cd22_incucyte %>%
  pivot_longer(
    cols = c(NB, m971, X819258_WT, nonint_4434),
    names_to = "condition",
    values_to = "confluence"
  )

# Add cellline to grouping
cd22_incucyte_summary <- cd22_incucyte_long %>%
  group_by(elapsed_hours, cellline, condition) %>%
  summarise(
    mean_conf = mean(confluence, na.rm = TRUE),
    sem_conf  = sd(confluence, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

cond_colors


cond_colors = c(jdb_palette("corona")[1],jdb_palette("corona")[2],jdb_palette("corona")[3],"darkgreen")
names(cond_colors) <- c("NB", "m971", "X819258_WT", "nonint_4434")

k562_df = df_summary %>% filter(cellline == "K562")

k562_plot = ggplot(k562_df, aes(x = elapsed_hours, y = mean_conf,
                                   color = condition, fill = condition)) +
  geom_errorbar(aes(ymin = mean_conf - sem_conf,
                    ymax = mean_conf + sem_conf),
                width = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cond_colors) +
  scale_fill_manual(values  = cond_colors) +
  pretty_plot() + L_border() #+
  #theme(axis.title = element_blank(),legend.position = "none")

k562_cd22_df = df_summary %>% filter(cellline == "K562 CD22")

k562_cd22_plot = ggplot(k562_cd22_df, aes(x = elapsed_hours, y = mean_conf,
                                color = condition, fill = condition)) +
  geom_errorbar(aes(ymin = mean_conf - sem_conf,
                    ymax = mean_conf + sem_conf),
                width = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cond_colors) +
  scale_fill_manual(values  = cond_colors) +
  pretty_plot() + L_border() #+
  #theme(axis.title = element_blank(),legend.position = "none")
k562_cd22_plot

rpmi_df = df_summary %>% filter(cellline == "RPMI-8226")

rpmi_plot = ggplot(rpmi_df, aes(x = elapsed_hours, y = mean_conf,
                                          color = condition, fill = condition)) +
  geom_errorbar(aes(ymin = mean_conf - sem_conf,
                    ymax = mean_conf + sem_conf),
                width = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cond_colors) +
  scale_fill_manual(values  = cond_colors) +
  pretty_plot() + L_border() #+
  #theme(axis.title = element_blank(),legend.position = "none")

rpmi_plot

k562_plot | k562_cd22_plot | rpmi_plot

