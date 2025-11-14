library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

dt_draws <- fread("../data/CD19_l112_l115_draws.csv") %>% data.frame()
dt_draws$mean_pos_activity <- (dt_draws$coculture_DG.75 + dt_draws$coculture_Raji + dt_draws$coculture_Ramos + dt_draws$coculture_K562_CD19)/4
dt_draws$activity_gain = dt_draws$mean_pos_activity - dt_draws$coculture_CAR_Only

## Keep negative control for one of the binders
dt_draws_filtered <- dt_draws %>% distinct(binder_name,.keep_all = TRUE)
dt_draws_filtered <- dt_draws_filtered %>% mutate(binder_label = case_when(
  binder_unique_id == "CD19_l112_FMC63" ~ "FMC63",
  binder_unique_id == "UTD" ~ "Control",
  binder_unique_id == "NB" ~ "Control",
  binder_unique_id == "CD19_l112_int_wildtype" ~ "Parental (l112)",
  binder_unique_id == "CD19_l115_int_wildtype" ~ "Parental (l115)",
  stringr::str_detect(binder_unique_id,"l112") ~ "CD19 l112",
  stringr::str_detect(binder_unique_id,"l115") ~ "CD19 l115",
))

dt_draws_filtered_sorted = dt_draws_filtered %>% 
  filter(binder_name != "FMC63", binder_name != "UTD", binder_name != "NB",
         binder_name != "l115_WT", binder_name != "l112_WT") %>%
  arrange(-activity_gain) %>%
  mutate(
    category = case_when(
      str_detect(binder_unique_id, "nonint") ~ "N",
      str_detect(binder_unique_id, "int") ~ "I",
      TRUE ~ ""
    ))

dt_draws_filtered_sorted_l112 = dt_draws_filtered_sorted %>%
  filter(binder_experiment == "CD19_l112") %>%
  mutate(new_id =paste0("C1.",category,ave(-activity_gain, category, FUN = function(x) rank(-x, ties.method = "min")-1)))

dt_draws_filtered_sorted_l115 = dt_draws_filtered_sorted %>%
  filter(binder_experiment == "CD19_l115") %>%
  mutate(new_id =paste0("C2.",category,ave(-activity_gain, category, FUN = function(x) rank(-x, ties.method = "min")-1)))
  
dt_draws_filtered_sorted_l112[,c("new_id","binder_unique_id")] %>% write.csv("../data/cd19_l112_mapping.csv",quote=FALSE)
dt_draws_filtered_sorted_l115[,c("new_id","binder_unique_id")] %>% write.csv("../data/cd19_l115_mapping.csv",quote=FALSE)

## Plot stimulation via Recombinant CD19 results
dt_activate <- fread("../data/CD19_l112_l115_strep_activation.csv") %>% data.frame()
dt_activate <- dt_activate %>% mutate(binder_label = case_when(
  binder_unique_id == "FMC63" ~ "FMC63",
  binder_unique_id == "UTD" ~ "Control",
  binder_unique_id == "NB" ~ "Control",
  binder_unique_id == "CD19_l112_int_wildtype" ~ "Parental (l112)",
  binder_unique_id == "CD19_l115_int_wildtype" ~ "Parental (l115)",
  stringr::str_detect(binder_unique_id,"l112") ~ "CD19 l112",
  stringr::str_detect(binder_unique_id,"l115") ~ "CD19 l115",
))
dt_activate$activity_gain_100nM <- dt_activate$CD69_100_nM - dt_activate$CD69_0_nM
dt_activate$activity_gain_10nM <- dt_activate$CD69_10_nM - dt_activate$CD69_0_nM
dt_activate$activity_gain_1nM <- dt_activate$CD69_1_nM - dt_activate$CD69_0_nM
dt_activate$activity_gain_0.1nM <- dt_activate$CD69_0.1_nM - dt_activate$CD69_0_nM

## Viz using staining at hr1 instead
dt_activate_merged = merge(dt_activate,dt_draws,by="binder_unique_id",all.x=TRUE)
dt_activate_merged <- dt_activate_merged %>% distinct(binder_unique_id,.keep_all = TRUE)

dt_draws_filtered
s1 = ggplot(dt_draws_filtered , aes(x = staining_10_nM, y = activity_gain,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("FMC63" = "orange",
               "CD19 l115" = "dodgerblue3",
               "Parental (l115)" = "dodgerblue",
               "CD19 l112" = "firebrick3",
               "Parental (l112)" = "firebrick1",
               "Control" == "gray"
               )) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining with 10nM (hr1)",y = "Coculture Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") +
  #geom_text(aes(label=binder_unique_id)) +
  ylim(0,100) + xlim(0,100)

s1_alt = ggplot(dt_draws_filtered , aes(x = staining_10_nM, y = coculture_Ramos-coculture_CAR_Only,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("FMC63" = "orange",
               "CD19 l115" = "dodgerblue3",
               "Parental (l115)" = "dodgerblue",
               "CD19 l112" = "firebrick3",
               "Parental (l112)" = "firebrick1",
               "Control" == "gray"
    )) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining with 10nM (hr1)",y = "Ramos Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") +
  ylim(0,100)
  #geom_text_repel(data=dt_draws_filtered %>% filter(staining_10_nM > 50),aes(label=binder_unique_id),max.overlaps = Inf)
  #geom_text(aes(label=binder_unique_id))
  #ylim(0,100) + xlim(0,100)

## Ramos look good
s1_alt
cowplot::ggsave2("../plots/CD19_10nM_hr1_vs_Ramos_gain_coculture.pdf",s1_alt,dpi=300,height=1.6,width=1.6)


dt_draws_filtered %>% arrange(-activity_gain)
dt_draws_filtered %>% filter(binder_unique_id=="CD19_l115_nonint_3533")

s1
cowplot::ggsave2("../plots/CD19_10nM_hr1_vs_activity_gain_coculture.pdf",s1,dpi=300,height=1.6,width=1.6)

t2 = ggplot(dt_activate_merged, aes(x = strep_10_nM, y = activity_gain_10nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("FMC63" = "orange",
               "CD19 l115" = "dodgerblue3",
               "Parental (l115)" = "dodgerblue",
               "CD19 l112" = "firebrick3",
               "Parental (l112)" = "firebrick1",
               "Control" == "gray"
    )) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining with 10nM (hr24)",y = "Recomb. CD19 (10nM) Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  #geom_text(aes(label=binder_unique_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none") +
  ylim(0,100)
t2
cowplot::ggsave2("../plots/CD19_10nM_hr24_vs_activity_gain_recomb_cd19.pdf",t2,dpi=300,height=1.6,width=1.6)

t2_alt = ggplot(dt_activate_merged, aes(x = strep_10_nM, y = activity_gain_10nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("FMC63" = "orange",
               "CD19 l115" = "dodgerblue3",
               "Parental (l115)" = "dodgerblue",
               "CD19 l112" = "firebrick3",
               "Parental (l112)" = "firebrick1",
               "Control" == "gray"
    )) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining with 10nM (hr24)",y = "Recomb. CD19 (10nM) Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  #geom_text(aes(label=binder_unique_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none") +
  ylim(0,100)
t2_alt


s1 | t2

#dt_activate_merged %>% filter(binder_name=="l115_WT")

dt_activate_merged

## Plot CD19 antigen dose response
dt_activate_merged$activity_gain_0nM = 0
dt_activate_merged_subset = dt_activate_merged %>% filter(
  binder_unique_id == "CD19_l115_int_wildtype" |
  binder_unique_id == "FMC63" |
  binder_unique_id == "CD19_l115_nonint_3533__CD19" |
  binder_unique_id == "NB"
) %>% 
  dplyr::select(binder_unique_id,activity_gain_0nM,activity_gain_0.1nM,activity_gain_1nM,activity_gain_10nM,activity_gain_100nM) %>%
  pivot_longer(cols = starts_with("activity_gain_"),names_to = "concentration_nM",values_to = "CD69_value") %>% mutate(
   concentration_nM = as.numeric(str_replace(str_replace(concentration_nM, "activity_gain_", ""),"nM",""))
  )
dt_activate_merged_subset <- dt_activate_merged_subset %>%
  mutate(
    concentration_nM_label = paste0(concentration_nM, "nM"),
    concentration_nM_factor = factor(concentration_nM, levels = c(0, 0.1, 1, 10, 100))
  )

dose_plot <- ggplot(dt_activate_merged_subset, aes(x = concentration_nM_factor, y = CD69_value, color = binder_unique_id)) +
  geom_point() +
  geom_line(aes(group = binder_unique_id)) + # Group is needed for geom_line when x is a factor
  labs(x = "Concentration (nM)") + # Relabel the axis
  scale_x_discrete(labels = unique(dt_activate_merged_subset$concentration_nM_factor)) + # Use the correct nM labels
  pretty_plot(fontsize = 8) +
  L_border() +
  scale_color_manual(values = c(
    "CD19_l115_int_wildtype"="dodgerblue",
    "CD19_l115_nonint_3533__CD19"="dodgerblue3",
    "FMC63"="orange",
    "NB"="black"
  )) + theme(legend.position = "none") +
  labs(x="CD19 Antigen (nM)",y="%CD69 (normalized to 0nM)")
dose_plot
cowplot::ggsave2("../plots/CD19_antigen_dose_plot.pdf",dose_plot,dpi=300,height=1.6,width=1.6)


## Plot Coculture results
p1 <- ggplot(dt_draws_filtered, aes(x = coculture_CAR_Only, y = mean_pos_activity,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%CD69+ with CAR Only",y = "%CD69+ with CD19+ Cells") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") +
  xlim(0,100) + ylim(0,100)
p1
cowplot::ggsave2("../plots/CD19_MPNN_Cocultures.pdf",p1,dpi=300,width=1.8,height=1.4)

## Plot at different nM again
p1a <- ggplot(dt_draws_filtered, aes(x = staining_0.1_nM, y = activity_gain,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining at 0.1nM",y = "Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
p1b <- ggplot(dt_draws_filtered, aes(x = staining_1_nM, y = activity_gain,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining at 1nM",y = "Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
p1c <- ggplot(dt_draws_filtered, aes(x = staining_10_nM, y = activity_gain,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining at 10nM",y = "Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
p1d <- ggplot(dt_draws_filtered, aes(x = staining_100_nM, y = activity_gain,color = binder_label)) +
  #geom_line(aes(group = binder_id), alpha = 0.6) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "Staining at 100nM",y = "Activity Gain") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")

test_p1p4 = cowplot::plot_grid(p1a,p1b,p1c,p1d,ncol=4)
test_p1p4
cowplot::ggsave2("../plots/CD19_MPNN_Coculture_nM.pdf",test_p1p4,dpi=300,width=1.8*4,height=1.4)


p2 <- ggplot(dt_activate, aes(x = CD69_0_nM, y = CD69_100_nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%CD69+ Unstim",y = "%CD69+ 100nM CD19") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") +
  xlim(0,100) + ylim(0,100)
p2
cowplot::ggsave2("../plots/CD19_MPNN_Recombinant_hr24.pdf",p2,dpi=300,width=1.8,height=1.4)

p2a <- ggplot(dt_activate, aes(x = strep_100_nM, y = activity_gain_100nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Strep at 100nM (hr24)",y = "%CD69+ (100nM-0nM)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") #+
  #xlim(0,100) + ylim(0,100)
p2a

p2b <- ggplot(dt_activate, aes(x = strep_10_nM, y = activity_gain_10nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Strep at 10nM (hr24)",y = "%CD69+ (10nM-0nM)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") #+
#xlim(0,100) + ylim(0,100)
p2b

p2c <- ggplot(dt_activate, aes(x = strep_1_nM, y = activity_gain_1nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Strep at 1nM (hr24)",y = "%CD69+ (1nM-0nM)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") #+
#xlim(0,100) + ylim(0,100)
p2b

p2d <- ggplot(dt_activate, aes(x = strep_0.1_nM, y = activity_gain_0.1nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Strep at 0.1nM (hr24)",y = "%CD69+ (0.1nM-0nM)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") #+
#xlim(0,100) + ylim(0,100)
p2b

test_p2p4 = cowplot::plot_grid(p2a,p2b,p2c,p2d,ncol=4)
cowplot::ggsave2("../plots/CD19_MPNN_Recombinant_hr24_nM.pdf",test_p2p4,dpi=300,width=1.8*4,height=1.4)

## Viz using staining at hr1 instead
dt_activate_merged = merge(dt_activate,dt_draws,by.x="binder_id",by.y="binder_name",all.x=TRUE)
dt_activate_merged <- dt_activate_merged %>% distinct(binder_id,.keep_all = TRUE)


p2e <- ggplot(dt_activate_merged, aes(x = staining_10_nM, y = activity_gain_100nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Strep at 0.1nM (hr24)",y = "%CD69+ (0.1nM-alone)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") 
p2e


dt_activate_merged %>% arrange(-strep_100_nM)

p3a = ggplot(dt_activate_merged, aes(x = staining_100_nM, y = strep_100_nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Stain at 100nM (hr1)",y = "%Strep at 100nM (hr24)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") 

p3b = ggplot(dt_activate_merged, aes(x = staining_10_nM, y = strep_10_nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Stain at 10nM (hr1)",y = "%Strep at 10nM (hr24)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") 

p3c = ggplot(dt_activate_merged, aes(x = staining_1_nM, y = strep_1_nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Stain at 1nM (hr1)",y = "%Strep at 1nM (hr24)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") 

p3d = ggplot(dt_activate_merged, aes(x = staining_0.1_nM, y = strep_0.1_nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "%Stain at 0.1nM (hr1)",y = "%Strep at 0.1nM (hr24)") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none") 

test_p3p4 = cowplot::plot_grid(p3a,p3b,p3c,p3d,ncol=4)
cowplot::ggsave2("../plots/CD19_MPNN_hr1_hr24_compraison_nM.pdf",test_p3p4,dpi=300,width=1.8*4,height=1.4)

dt_activate_merged$bound_hr1_hr24_100nM = dt_activate_merged$staining_100_nM - dt_activate_merged$strep_100_nM
dt_activate_merged$bound_hr1_hr24_10nM = dt_activate_merged$staining_10_nM - dt_activate_merged$strep_10_nM
dt_activate_merged$bound_hr1_hr24_1nM = dt_activate_merged$staining_1_nM - dt_activate_merged$strep_1_nM
dt_activate_merged$bound_hr1_hr24_0.1nM = dt_activate_merged$staining_0.1_nM - dt_activate_merged$strep_0.1_nM

dt_activate_merged$cell_recomb_gain_100nM = dt_activate_merged$activity_gain_100nM - dt_activate_merged$activity_gain
dt_activate_merged$cell_recomb_gain_10nM = dt_activate_merged$activity_gain_10nM - dt_activate_merged$activity_gain
dt_activate_merged$cell_recomb_gain_1nM = dt_activate_merged$activity_gain_1nM - dt_activate_merged$activity_gain
dt_activate_merged$cell_recomb_gain_0.1nM = dt_activate_merged$activity_gain_0.1nM - dt_activate_merged$activity_gain

t1 = ggplot(dt_activate_merged, aes(x = strep_1_nM, y = activity_gain,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "staining_10_nM",y = "%CD69 coculture") +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_text(aes(label=binder_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none") + 
  ylim(0,100)

t1

s1 | t2

t1 | t2

s1 | t2

q1 = ggplot(dt_activate_merged, aes(x = strep_1_nM, y = activity_gain,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "staining_10_nM",y = "%CD69 coculture") +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_text(aes(label=binder_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none") #+ 
  #ylim(0,100)

q2 = ggplot(dt_activate_merged, aes(x = strep_1_nM, y = cell_recomb_gain_100nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "strep_1_nM",y = "activity_gain_10nM (recomb)") +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_text(aes(label=binder_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none")# +
  #ylim(0,100)

q1 | q2

q3 = ggplot(dt_activate_merged, aes(x = activity_gain_100nM, y = cell_recomb_gain_100nM,color = binder_label)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FMC63" = "orange", "CD19 l115" = "dodgerblue3", "CD19 l112" = "firebrick3","Control" == "gray")) +
  #geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  labs(x = "strep_1_nM",y = "activity_gain_10nM (recomb)") +
  pretty_plot(fontsize = 8) + L_border() + 
  geom_text(aes(label=binder_id)) +
  #geom_hline(mean(dt_activate_merged$activity_gain_10nM)) +
  theme(legend.position = "none")# +
#ylim(0,100)
q3
