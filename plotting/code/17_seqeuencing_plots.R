library(BuenColors)
library(dplyr)
library(data.table)

bcma_tpm <- read.csv("../data/BCMA_bindcraft_tpm.csv")
cd19_big_tpm <- read.csv("../data/CD19_big_bindcraft_tpm.csv") %>% rename(binder_id = target_id)
cd19_pd_tpm <- read.csv("../data/CD19_PartialDiffusion_tpm.csv")  %>% rename(binder_id = target_id)
cd22_tpm <- read.csv("../data/CD22_bindcraft_tpm.csv")

## Highlight the binders selected for testing
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

bcma_tpm <- bcma_tpm %>% mutate(is_selected = binder_id %in% bcma_tested_in_cars)
cd19_big_tpm <- cd19_big_tpm %>% mutate(is_selected = binder_id %in% cd19_big_tested_in_cars)
cd19_pd_tpm <- cd19_pd_tpm %>% mutate(is_selected = binder_id %in% cd19_pd_tested_in_cars)
cd22_tpm <- cd22_tpm %>% mutate(is_selected = binder_id %in% cd22_tested_in_cars)

## Make fold-change plots
plot_fc <- function(df,pcol,curcol) {
  # Log transform tpm and fold change
  df$log1p_tpm <- log10(df %>% pull({{curcol}})+1)
  df$logFc <- log2(df %>% pull({{curcol}}) / df %>% pull({{pcol}}))
  # Cap negative FC at -10
  df$logFc <- pmax(df$logFc,-10)
  # Plot
  fc_plot <- ggplot(df, aes(x = log1p_tpm, y = logFc, color = is_selected)) + 
    geom_point() +
    scale_color_manual(values = c("dodgerblue3","firebrick")) + 
    pretty_plot(fontsize = 8) + L_border() + 
    theme(legend.position = "none") +
    geom_hline(yintercept =0, linetype = 2) +
    labs(x="log10(TPM+1)",y="logFC(TPM/Parental)")
  return(fc_plot)
}

bcma_m1000_plot <- plot_fc(bcma_tpm,"BCMA_parental_mean","BCMA_M1000_mean")
bcma_m100_plot <- plot_fc(bcma_tpm,"BCMA_parental_mean","BCMA_M100_mean")
cd22_m100_plot <- plot_fc(cd22_tpm,"CD22_parental_mean","CD22_100_mean")
cd19_big_m100_plot <- plot_fc(cd19_big_tpm,"CD19_Big_Parental_mean","CD19_Big_MACS_100_mean")
cd19_pd_m50_plot <- plot_fc(cd19_pd_tpm,"CD19_PartialDiffusion_P_mean","CD19_PartialDiffusion_MACS_50_mean")

cowplot::ggsave2(bcma_m1000_plot, file = "../plots/tpm_FC_bcma_1000nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(bcma_m100_plot, file = "../plots/tpm_FC_bcma_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd22_m100_plot, file = "../plots/tpm_FC_cd22_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd19_big_m100_plot, file = "../plots/tpm_FC_cd19_big_100nM.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(cd19_pd_m50_plot, file = "../plots/tpm_FC_cd19_pd_50nM.pdf", width = 1.3, height = 1.3)

## Make lineplots
plot_binder_tpm <- function(df, condition_labels) {
  plot_cols <- intersect(names(condition_labels), names(df))
  df_long <- df %>%
    pivot_longer(cols = all_of(plot_cols),names_to = "condition",values_to = "TPM") %>%
    arrange(is_selected)
  df_long$condition <- factor(df_long$condition, levels = plot_cols)
  p <- ggplot(df_long, aes(x = condition, y = TPM, group = binder_id, color = is_selected)) +
    geom_line(aes(group = binder_id), alpha = 0.6) +
    geom_point(size = 3) +
    scale_x_discrete(labels = condition_labels[plot_cols]) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"),
                       labels = c("TRUE" = "Selected Binders", "FALSE" = "Other Binders")) +
    labs(x = NULL,y = "TPM",color = "Group") +
    pretty_plot(fontsize = 8) + L_border() + 
    theme(legend.position = "none")
  return(p)
}

# Generate the TPM plots
bcma_condition_map <- c(
  BCMA_parental_mean = "Parental",
  BCMA_M1000_mean = "MACS 1000nM",
  BCMA_M500_mean = "MACS 500nM",
  BCMA_M100_mean = "MACS 100nM"
)
bcma_all_tpm_plot <- plot_binder_tpm(bcma_tpm, bcma_condition_map)
bcma_all_tpm_plot
cowplot::ggsave2(bcma_all_tpm_plot, file = "../plots/tpm_bcma_all_nM.pdf", width = 1.3*4, height = 1.3)

cd19_big_condition_map <- c(
  CD19_Big_Parental_mean = "Parental",
  CD19_Big_MACS_1000_mean = "MACS 1000nM",
  CD19_Big_MACS_100_mean = "MACS 100nM"
)
cd19_big_all_tpm_plot <- plot_binder_tpm(cd19_big_tpm, cd19_big_condition_map)
cd19_big_all_tpm_plot
cowplot::ggsave2(cd19_big_all_tpm_plot, file = "../plots/tpm_CD19_big_all_nM.pdf", width = 1.3*3, height = 1.3)

cd22_big_condition_map <- c(
  CD22_parental_mean = "Parental",
  CD22_1000_mean = "MACS 1000nM",
  CD22_500_mean = "MACS 500nM",
  CD22_100_mean = "MACS 100nM"
)
cd22_tpm_plot <- plot_binder_tpm(cd22_tpm, cd22_big_condition_map)
cowplot::ggsave2(cd19_big_all_tpm_plot, file = "../plots/tpm_CD22_all_nM.pdf", width = 1.3*4, height = 1.3)







