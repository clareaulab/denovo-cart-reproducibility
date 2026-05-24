library(BuenColors)
library(dplyr)
library(data.table)

bcma_rfd_tpm <- read.csv("../data/BCMA_e3_fold_conditioned_tpm.csv")
bcma_tpm <- read.csv("../data/BCMA_bindcraft_tpm.csv")
cd19_big_tpm <- read.csv("../data/CD19_big_bindcraft_tpm.csv") %>% dplyr::rename(binder_id = target_id)
cd19_pd_tpm <- read.csv("../data/CD19_PartialDiffusion_tpm.csv")  %>% dplyr::rename(binder_id = target_id)
cd22_tpm <- read.csv("../data/CD22_bindcraft_tpm.csv")
cd19_e3_tpm <- read.csv("../data/CD19_e3_fold_conditioned_tpm.csv")

## Highlight the binders selected for testing
bcma_rfd_tested_in_cars <- c(
  "BCMAE3FoldConditionedAllFolds_1XU2ChainR_R18-R20-R26_50-100_i16_51_dldesign_0",
  "BCMAE3FoldConditionedAllFolds_1XU2ChainR_R18-R20-R26_50-100_i9_54_dldesign_0",
  "BCMAE3FoldConditionedAllFolds_1XU2ChainR_R18-R20-R26_50-100_i1_29_dldesign_0",
  "BCMAE3FoldConditioned_1XU2ChainR_R18-R20-R26_50-100_i55_32_dldesign_0",
  "BCMAE3FoldConditioned_1XU2ChainR_R18-R20-R26_50-100_i68_76_dldesign_0"
)

bcma_tested_in_cars <- c(
  "BCMA_l53_s847868_mpnn2","BCMA_l58_s298275_mpnn2","BCMA_l51_s395146_mpnn6",
  "BCMA_l60_s719112_mpnn1","BCMA_l64_s766805_mpnn3","BCMA_l60_s542430_mpnn1",
  "BCMA_l59_s561726_mpnn1","BCMA_l64_s814818_mpnn1","BCMA_l55_s655034_mpnn2",
  "BCMA_l59_s903443_mpnn2","BCMA_l64_s766805_mpnn1","BCMA_l61_s2373_mpnn5",
  "BCMA_l61_s2373_mpnn8","BCMA_l63_s961592_mpnn9","BCMA_l55_s655034_mpnn1",
  "BCMA_l54_s816947_mpnn2","BCMA_l61_s550568_mpnn3"
)
cd19_big_tested_in_cars <- c(
  "CD19_Big_l186_s727792_mpnn20","CD19_Big_l115_s450328_mpnn6","CD19_Big_l112_s592202_mpnn2"
)
cd19_pd_tested_in_cars <- c(
  "CD19_par186_0217__CD19","CD19_par186_1283__CD19","CD19_par186_1362__CD19"
)
cd22_tested_in_cars <- c(
  "CD22_l61_s819258_mpnn1","CD22_l61_s819258_mpnn2","CD22_l64_s571098_mpnn2",
  "CD22_l52_s770345_mpnn1","CD22_l59_s943680_mpnn1"
)
cd19_e3_tested_in_cars <- c(
  "i47_37_dldesign_0","i40_65_dldesign_0","i89_343_dldesign_0",
  "i55_32_dldesign_7","i18_63_dldesign_7","i15_32_dldesign_0"
)

bcma_rfd_tpm <- bcma_rfd_tpm %>% mutate(is_selected = binder_id %in% bcma_rfd_tested_in_cars)
bcma_tpm <- bcma_tpm %>% mutate(is_selected = binder_id %in% bcma_tested_in_cars)
cd19_big_tpm <- cd19_big_tpm %>% mutate(is_selected = binder_id %in% cd19_big_tested_in_cars)
cd19_pd_tpm <- cd19_pd_tpm %>% mutate(is_selected = binder_id %in% cd19_pd_tested_in_cars)
cd22_tpm <- cd22_tpm %>% mutate(is_selected = binder_id %in% cd22_tested_in_cars)
cd19_e3_tpm <- cd19_e3_tpm %>% mutate(is_selected = binder_id %in% cd19_e3_tested_in_cars)


## Make fold-change plots
plot_fc <- function(df,pcol,curcol) {
  # Log transform tpm and fold change
  df$log1p_tpm <- log10(df %>% pull({{curcol}})+1)
  df$logFc <- log2(df %>% pull({{curcol}}) / df %>% pull({{pcol}}))
  # Cap negative FC at -10
  df$logFc <- pmax(df$logFc,-10)
  print(df)
  # Plot
  fc_plot <- ggplot(df, aes(x = log1p_tpm, y = logFc, color = is_selected)) + 
    geom_point() +
    scale_color_manual(values = c("lightgray","dodgerblue3")) + 
    pretty_plot(fontsize = 8) + L_border() + 
    theme(legend.position = "none") +
    geom_hline(yintercept =0, linetype = 2) +
    labs(x="log10(TPM+1)",y="logFC(TPM/Parental)")
  return(fc_plot)
}

plot_fc_alt <- function(df,pcol,curcol) {
  # Log transform tpm and fold change
  df$log1p_tpm <- log10(df %>% pull({{curcol}})+1)
  df$logFc <- log2(df %>% pull({{curcol}}) / df %>% pull({{pcol}}))
  # Cap negative FC at -10
  df$logFc <- pmax(df$logFc,-10)
  print(df %>% filter(is_selected))
  
  df$logFC_100 = log2(df %>% pull(BCMA_M100_mean) / df %>% pull(BCMA_parental_mean))
  df$logFC_500 = log2(df %>% pull(BCMA_M500_mean) / df %>% pull(BCMA_parental_mean))
  df$logFC_1000 = log2(df %>% pull(BCMA_M1000_mean) / df %>% pull(BCMA_parental_mean))
  
  df$is_subfc_blue = df$binder_id %in% c("BCMA_l54_s816947_mpnn2","BCMA_l61_s550568_mpnn3")
  
  # Plot
  fc_plot <- ggplot(df, aes(x = log1p_tpm, y = logFc, color = is_subfc_blue)) + 
    geom_point() +
    scale_color_manual(values = c("lightgray","dodgerblue3")) + 
    pretty_plot(fontsize = 8) + L_border() + 
    theme(legend.position = "none") +
    geom_hline(yintercept =0, linetype = 2) +
    labs(x="log10(TPM+1)",y="logFC(TPM/Parental)")
  return(fc_plot)
}

## The two blue dots with logFC <0
#"BCMA_l54_s816947_mpnn2"
#"BCMA_l61_s550568_mpnn3"

write.csv(bcma_rfd_tpm,"../data/bcma_rfd_tpm.csv",row.names = FALSE)

plot_fc_alt(bcma_tpm,"BCMA_parental_mean","BCMA_M1000_mean")
plot_fc_alt(bcma_tpm,"BCMA_parental_mean","BCMA_M500_mean")
plot_fc_alt(bcma_tpm,"BCMA_parental_mean","BCMA_M100_mean")

bcma_rfd_tpm$logFC_1000 = log2(bcma_rfd_tpm %>% pull(BCMA_1000_mean) / bcma_rfd_tpm %>% pull(BCMA_parental_mean))
bcma_rfd_tpm$logFC_100 = log2(bcma_rfd_tpm %>% pull(BCMA_100_mean) / bcma_rfd_tpm %>% pull(BCMA_parental_mean))

# ggplot(bcma_rfd_tpm,aes(x=logFC_1000,y=logFC_100,color=is_selected)) + 
#   geom_point() +
#   pretty_plot() + L_border() +
#   geom_hline(yintercept=0) + geom_vline(xintercept=0)

# ggplot(bcma_rfd_tpm,aes(x=log10(BCMA_parental_mean+1),y=logFC_100,color=is_selected)) + 
#   geom_point() +
#   pretty_plot() + L_border() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = log10(1/dim(bcma_rfd_tpm)[1]), linetype="dashed") +
#   scale_color_manual(values = c("lightgray","dodgerblue3"))

bcma_rfd_tpm$BCMA_parental_occupancy = bcma_rfd_tpm$BCMA_parental_mean/sum(bcma_rfd_tpm$BCMA_parental_mean)
bcma_rfd_tpm$BCMA_1000_occupancy = bcma_rfd_tpm$BCMA_1000_mean/sum(bcma_rfd_tpm$BCMA_1000_mean)
bcma_rfd_tpm$BCMA_100_occupancy = bcma_rfd_tpm$BCMA_100_mean/sum(bcma_rfd_tpm$BCMA_100_mean)
bcma_rfd_tpm$BCMA_1000_occupancy_FC = bcma_rfd_tpm$BCMA_1000_occupancy/bcma_rfd_tpm$BCMA_parental_occupancy
bcma_rfd_tpm$BCMA_100_occupancy_FC =bcma_rfd_tpm$BCMA_100_occupancy/bcma_rfd_tpm$BCMA_parental_occupancy

bli_tested = c(
  "BCMAE3FoldConditionedAllFolds_1XU2ChainR_R18-R20-R26_50-100_i9_54_dldesign_0",
  "BCMAE3FoldConditioned_1XU2ChainR_R18-R20-R26_50-100_i55_32_dldesign_0"
)

bcma_rfd_tpm = bcma_rfd_tpm %>% mutate(
  binder_id_short = str_extract(binder_id, "([^_]+_){3}[^_]+$")
)

bcma_rfd_tpm_m100_fc_plot = ggplot(bcma_rfd_tpm,aes(x=BCMA_parental_occupancy,y=logFC_100,color=is_selected,label=binder_id_short)) +
  geom_point() +
  pretty_plot() + L_border() +
  geom_hline(yintercept = 0) +
  scale_x_continuous(labels = label_percent()) +
  geom_vline(xintercept = 1/dim(bcma_rfd_tpm)[1], linetype="dashed") +
  labs(x="%Total CPM in Parental",y="logFC(100nM/Parental)") +
  scale_color_manual(values = c("lightgray","dodgerblue3")) +
  geom_label_repel(data = bcma_rfd_tpm %>% filter(logFC_100 > 0)) +
  theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 5))
bcma_rfd_tpm_m100_fc_plot

cowplot::ggsave2("../plots/tpm_FC_bcma_rfd_percent_parental_100nM.pdf",bcma_rfd_tpm_m100_fc_plot,dpi=300,height = 1,width=1)

## Read in CAR data
rfd5_car_readout <- read.csv("../data/RFD5_CAR-J_CD69_cocultures.csv",check.names = FALSE)
rfd5_car_readout_df <- rfd5_car_readout %>%
  column_to_rownames("Cellline") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("binder_id") %>%
  mutate(across(c(K562, MM1S, Raji, `CAR alone`), as.numeric)) %>%
  mutate(avg_BCMA_pos = (MM1S + Raji) / 2)
rfd5_car_readout_df = rfd5_car_readout_df %>% mutate(
  binder_id_short = case_when(
    binder_id == "No binder" ~ "No binder",
    binder_id == "Abecma" ~ "Abecma",
    TRUE ~ str_extract(binder_id, "([^_]+_){3}[^_]+$")
  )
)

rfd5_car_readout_df$is_selected = rfd5_car_readout_df$binder_id %in% bcma_rfd_tested_in_cars
rfd5_car_readout_df_barplot = ggplot(rfd5_car_readout_df,aes(x=reorder(binder_id_short,-avg_BCMA_pos),y=avg_BCMA_pos)) +
  geom_col(aes(fill=is_selected)) +
  geom_hline(yintercept = rfd5_car_readout_df %>% filter(binder_id_short == "No binder") %>% pull(avg_BCMA_pos), linetype="dashed") +
  scale_y_continuous(expand=c(0,Inf)) +
  pretty_plot() + L_border() +
  scale_fill_manual(values=c("gray","dodgerblue3")) +
  theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text = element_text(size = 5),axis.text.x = element_blank(), axis.ticks.x = element_blank()
        )
rfd5_car_readout_df_barplot
cowplot::ggsave2("../plots/rfd5_car_readout_df_barplot.pdf",rfd5_car_readout_df_barplot,dpi=300,height = 1,width=1)

  

#BCMA_parental_occupancy < (1/dim(bcma_rfd_tpm)[1])
# 
bcma_rfd_tpm_filtered = bcma_rfd_tpm #%>% filter(logFC_100 > 1, logFC_1000 > 1)
bcma_rfd_tpm_m100_fc_plot_alt = ggplot(bcma_rfd_tpm_filtered,aes(x=logFC_100,y=logFC_1000,color=is_selected,label=binder_id_short)) +
  geom_point() +
  pretty_plot() + L_border() +
  geom_hline(yintercept = 0) +
  scale_x_continuous(labels = label_percent()) +
  geom_vline(xintercept = 1/dim(bcma_rfd_tpm)[1], linetype="dashed") +
  labs(x="%Total CPM in Parental",y="logFC(100nM/Parental)") +
  scale_color_manual(values = c("lightgray","dodgerblue3")) +
  #geom_text_repel() +
  #geom_text_repel(data = bcma_rfd_tpm %>% filter(logFC_100 > 0)) +
  theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 5))
bcma_rfd_tpm_m100_fc_plot_alt



ggplot(bcma_rfd_tpm,aes(x=BCMA_100_occupancy,fill=is_selected)) +
geom_histogram()

bcma_rfd_m1000_plot <- plot_fc(bcma_rfd_tpm,"BCMA_parental_mean", "BCMA_1000_mean")
bcma_rfd_m100_plot <- plot_fc(bcma_rfd_tpm,"BCMA_parental_mean", "BCMA_100_mean")
bcma_m1000_plot <- plot_fc(bcma_tpm,"BCMA_parental_mean","BCMA_M1000_mean")
bcma_m100_plot <- plot_fc(bcma_tpm,"BCMA_parental_mean","BCMA_M100_mean")
cd22_m100_plot <- plot_fc(cd22_tpm,"CD22_parental_mean","CD22_100_mean")
cd19_big_m100_plot <- plot_fc(cd19_big_tpm,"CD19_Big_Parental_mean","CD19_Big_MACS_100_mean")
cd19_pd_m50_plot <- plot_fc(cd19_pd_tpm,"CD19_PartialDiffusion_P_mean","CD19_PartialDiffusion_MACS_50_mean")
cd19_e3_m3f_plot <- plot_fc(cd19_e3_tpm,"CD19_parental_mean","CD19_M3F_mean")

cowplot::ggsave2(bcma_rfd_m1000_plot, file = "../plots/tpm_FC_bcma_rfd_1000nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(bcma_rfd_m100_plot, file = "../plots/tpm_FC_bcma_rfd_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(bcma_m1000_plot, file = "../plots/tpm_FC_bcma_1000nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(bcma_m100_plot, file = "../plots/tpm_FC_bcma_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd22_m100_plot, file = "../plots/tpm_FC_cd22_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd19_big_m100_plot, file = "../plots/tpm_FC_cd19_big_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd19_pd_m50_plot, file = "../plots/tpm_FC_cd19_pd_50nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd19_e3_m3f_plot, file = "../plots/tpm_FC_cd19_m3f_1000nM.pdf", width = 1.3, height = 1.3)

bcma_m100_plot_small = bcma_m100_plot + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 5))
cowplot::ggsave2(bcma_m100_plot_small, file = "../plots/tpm_FC_bcma_100nM_small.pdf", width = 1, height = 1)

bcma_tpm %>% filter(is_selected, logFC_100 < 0)
"BCMA_l54_s816947_mpnn2"
"BCMA_l61_s550568_mpnn3"
under_fc100_bcma_binders = c("BCMA_l54_s816947_mpnn2","BCMA_l61_s550568_mpnn3")
bcma_tpm$under_fc100_bcma_binders = bcma_tpm$binder_id %in% under_fc100_bcma_binders
bcma_tpm$logFC_500 = log2(bcma_tpm %>% pull(BCMA_M500_mean) / bcma_tpm %>% pull(BCMA_parental_mean))

ggplot(bcma_tpm,aes(x=log10(BCMA_M1000_mean+1),y=logFC_1000,color=under_fc100_bcma_binders)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_M1000_mean+1),y=logFC_1000,color=under_fc100_bcma_binders)) + geom_point()
ggplot(bcma_tpm,aes(x=logFC_100,y=logFC_1000,color=is_selected)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_parental_mean+1),y=logFC_1000,color=under_fc100_bcma_binders)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_parental_mean+1),y=logFC_1000,color=is_selected)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_M1000_mean+1),y=logFC_1000,color=under_fc100_bcma_binders)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_M1000_mean+1),y=logFC_1000,color=is_selected)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_M500_mean+1),y=logFC_500,color=under_fc100_bcma_binders)) + geom_point()
ggplot(bcma_tpm,aes(x=log10(BCMA_M500_mean+1),y=logFC_500,color=is_selected)) + geom_point()

write.csv(bcma_tpm,"../data/bcma_bc_tpm.csv",row.names = FALSE)

bcma_m1000_plot

bcma_m100_plot
## Make lineplots
plot_binder_tpm <- function(df, condition_labels) {
  plot_cols <- intersect(names(condition_labels), names(df))
  df_long <- df %>%
    pivot_longer(cols = all_of(plot_cols),names_to = "condition",values_to = "TPM") %>%
    arrange(is_selected)
  df_long$condition <- factor(df_long$condition, levels = plot_cols)
  p <- ggplot(df_long, aes(x = condition, y = TPM, group = binder_id, color = is_selected)) +
    geom_line(aes(group = binder_id), alpha = 0.6) +
    geom_point(size = 1) +
    scale_x_discrete(labels = condition_labels[plot_cols]) +
    scale_color_manual(values = c("TRUE" = "dodgerblue3", "FALSE" = "lightgrey"),
                       labels = c("TRUE" = "Selected Binders", "FALSE" = "Other Binders")) +
    labs(x = "library",y = "TPM",color = "Group") +
    pretty_plot(fontsize = 8) + L_border() + 
    theme(legend.position = "none")
  return(p)
}

# Generate the TPM plots
bcma_rfd_condition_map <- c(
  BCMA_parental_mean = "Parental",
  BCMA_1000_mean = "MACS 1000nM",
  BCMA_100_mean = "MACS 100nM"
)
bcma_rfd_all_tpm_plot <- plot_binder_tpm(bcma_rfd_tpm, bcma_rfd_condition_map)
bcma_rfd_all_tpm_plot
cowplot::ggsave2(bcma_rfd_all_tpm_plot, file = "../plots/tpm_bcma_rfd_all_nM.pdf", width = 3.5, height = 1.3)

bcma_rfd_all_tpm_plot


bcma_condition_map <- c(
  BCMA_parental_mean = "Parental",
  BCMA_M1000_mean = "MACS 1000nM",
  BCMA_M500_mean = "MACS 500nM",
  BCMA_M100_mean = "MACS 100nM"
)
bcma_all_tpm_plot <- plot_binder_tpm(bcma_tpm, bcma_condition_map)
bcma_all_tpm_plot
cowplot::ggsave2(bcma_all_tpm_plot, file = "../plots/tpm_bcma_all_nM.pdf", width = 3.5, height = 1.3)

cd19_big_condition_map <- c(
  CD19_Big_Parental_mean = "Parental",
  CD19_Big_MACS_1000_mean = "MACS 1000nM",
  CD19_Big_MACS_100_mean = "MACS 100nM"
)
cd19_big_all_tpm_plot <- plot_binder_tpm(cd19_big_tpm, cd19_big_condition_map)
cd19_big_all_tpm_plot
cowplot::ggsave2(cd19_big_all_tpm_plot, file = "../plots/tpm_CD19_big_all_nM.pdf", width = 3.5, height = 1.3)

cd22_big_condition_map <- c(
  CD22_parental_mean = "Parental",
  CD22_1000_mean = "MACS 1000nM",
  CD22_500_mean = "MACS 500nM",
  CD22_100_mean = "MACS 100nM"
)
cd22_tpm_plot <- plot_binder_tpm(cd22_tpm, cd22_big_condition_map)
cowplot::ggsave2(cd22_tpm_plot, file = "../plots/tpm_CD22_all_nM.pdf", width = 3.5, height = 1.3)

cd19_e3_condition_map <- c(
  CD19_parental_mean = "Parental",
  CD19_M1_mean = "MACS 1000nM",
  CD19_M1F_mean = "FACSx1",
  CD19_M2F_mean = "FACSx2",
  CD19_M3F_mean = "FACSx3"
)
cd19_e3_tpm_plot <- plot_binder_tpm(cd19_e3_tpm, cd19_e3_condition_map)
cd19_e3_tpm_plot
cowplot::ggsave2(cd19_e3_tpm_plot, file = "../plots/tpm_CD19_e3_all_nM.pdf", width = 3.5, height = 1.3)





