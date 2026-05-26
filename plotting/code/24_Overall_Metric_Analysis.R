library(BuenColors)
library(dplyr)
library(data.table)
library(yardstick)
library(stringr)
library(ggbeeswarm)
library(cowplot)

#dt = fread("../data/2025_10_24_all_binders_with_labels.csv") %>% data.frame()
#dt = fread("../data/2025_11_15_all_binders_with_labels.csv") %>% data.frame()
#dt = fread("../data/2025_12_14_all_binders_with_labels.csv") %>% data.frame()
dt = fread("../data/2026_05_19_all_binders_with_labels.csv") %>% data.frame()


dt = dt %>% mutate(
  binder_by_YSD_1000nM = as.factor(binder_by_YSD_1000nM),
  binder_by_YSD_500nM = as.factor(binder_by_YSD_500nM),
  binder_by_YSD_100nM = as.factor(binder_by_YSD_100nM),
  binder_by_estimated_KD_1000nM = as.factor(binder_by_estimated_KD_1000nM),
  binder_by_estimated_KD_500nM = as.factor(binder_by_estimated_KD_500nM),
  binder_by_estimated_KD_100nM = as.factor(binder_by_estimated_KD_100nM),
  binder_by_estimated_KD_10nM = as.factor(binder_by_estimated_KD_10nM),
  
  Average_ipSAE = as.numeric(Average_ipSAE),
  Average_pae_interaction = as.numeric(Average_pae_interaction),
  Average_iptm = as.numeric(Average_iptm),
  Average_complex_plddt = as.numeric(Average_complex_plddt),
  
  below_pae_10 = Average_pae_interaction < 10,
  above_iptm_85 = Average_iptm >= 0.85,
  above_ipSAE_85 = Average_ipSAE >= 0.85,
  
  above_ipSAE_max_85 = Average_ipSAE_max >= 0.85,
  above_ipSAE_min_85 = Average_ipSAE_min >= 0.85,
  above_ipSAE_a_to_b_85 = Average_Chn1_A_to_Chn2_B_ipSAE >= 0.85,
  above_ipSAE_b_to_a_85 = Average_Chn1_B_to_Chn2_A_ipSAE >= 0.85,
  
  below_pde_1 = Average_complex_pde < 1,
  below_ipde_1 = Average_complex_ipde < 1,

  antigen = case_when(
    str_detect(campaign, regex("CD19", ignore_case = TRUE)) ~ "CD19",
    str_detect(campaign, regex("BCMA", ignore_case = TRUE)) ~ "BCMA",
    str_detect(campaign, regex("CD22", ignore_case = TRUE)) ~ "CD22",
    TRUE ~ "Other Target"
  ),
  
  binder_evolution_type = case_when(
    str_detect(binder_name, regex("wildtype", ignore_case = TRUE)) ~ "parental",
    str_detect(binder_name, regex("_int_", ignore_case = TRUE)) ~ "interface",
    str_detect(binder_name, regex("_nonint_", ignore_case = TRUE)) ~ "non-interface",
    TRUE ~ NA
  ),
  
  high_tonic = coculture.CAR.Only > 25
)

dim(dt)

## Check alanine percentage
dt = dt %>%
  mutate(
    sequence_length = str_length(binder_sequence),
    alanine_count = str_count(binder_sequence, 'A'),
    leucine_count = str_count(binder_sequence, 'L'),
    alanine_percentage = (alanine_count / sequence_length) * 100,
    leucine_percentage = (leucine_count / sequence_length) * 100,
    is_bindcraft = str_detect(campaign,"Bind"),
  )

dt$campaign_short_replaced = str_replace_all(dt$campaign_short, " \\(", "\n(")

table(dt$campaign_short_replaced)

dt_ysd = dt %>% filter(!is.na(binder_by_YSD_1000nM))

dt$campaign_short_replaced = factor(dt$campaign_short_replaced,levels=c(
  "RFD 1\n(BCMA)","RFD 2\n(BCMA)","RFD 3\n(BCMA)","RFD 4\n(BCMA)","RFD 5\n(BCMA)","BC 1\n(BCMA)",
  "RFD 1\n(CD19)","RFD 2\n(CD19)","RFD 3\n(CD19)","RFD 4\n(CD19)","BC 1\n(CD19)","BC 2\n(CD19)",
  #"RFD3 1\n(CD19)", "BG VHH\n(CD19)",
  "RFD 1\n(CD22)", "BC 1\n(CD22)"
))

table(dt$campaign)

dt_ysd = dt %>% filter(!is.na(binder_by_YSD_1000nM))
## Alanine plot
ala_plot = ggplot(dt_ysd %>% arrange(binder_by_YSD_1000nM), aes(x = campaign_short_replaced, y=alanine_percentage, color=binder_by_YSD_1000nM)) +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  L_border() +
  scale_color_manual(values = c("gray", "firebrick")) + 
  labs(x="Yeast Campaign",y="% Alanine in Sequence") +
  theme(legend.position = "none")
#scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink", "forestgreen","palegreen")) +
#labs(y="NetMHC Class I Burden",x="Antigen")
ala_plot

cowplot::ggsave2("../plots/alanine_percentages_by_campaign.pdf",ala_plot,dpi=300,width=6.5,height=1.8,unit="in")

ggplot(data=dt %>% filter(alanine_count < 10, Average_monomer_lDDT > 80),aes(y=campaign,x=Average_ipSAE,color=binder_by_YSD_100nM)) +
  geom_jitter() +
  geom_boxplot() +
  geom_vline(xintercept=0.85)

ggplot(data=dt %>% filter(Average_monomer_lDDT > 80, alanine_percentage < 20),aes(y=campaign,x=Average_ipSAE,color=binder_by_YSD_1000nM)) +
  geom_jitter() +
  geom_vline(xintercept=0.85)

ggplot(dt %>% arrange(binder_by_YSD_1000nM),aes(x=alanine_percentage,y=Average_Binder_BetaSheet.,color=binder_by_YSD_1000nM)) + geom_point()

# find_significant_diffs <- function(df, group_col, alpha = 0.05) {
#   # Get numeric columns (excluding the grouping column)
#   num_cols <- names(df)[sapply(df, is.numeric) & names(df) != group_col]
#   
#   # Test each numeric column
#   results <- lapply(num_cols, function(col) {
#     tryCatch({
#       test <- t.test(df[[col]] ~ df[[group_col]])
#       data.frame(
#         variable = col,
#         p_value = test$p.value,
#         mean_true = mean(df[[col]][df[[group_col]] == TRUE], na.rm = TRUE),
#         mean_false = mean(df[[col]][df[[group_col]] == FALSE], na.rm = TRUE),
#         significant = test$p.value < alpha,
#         stringsAsFactors = FALSE
#       )
#     }, error = function(e) {
#       # Handle cases with no variance or other errors
#       data.frame(
#         variable = col,
#         p_value = NA,
#         mean_true = mean(df[[col]][df[[group_col]] == TRUE], na.rm = TRUE),
#         mean_false = mean(df[[col]][df[[group_col]] == FALSE], na.rm = TRUE),
#         significant = FALSE,
#         stringsAsFactors = FALSE
#       )
#     })
#   })
#   
#   # Combine into single dataframe
#   do.call(rbind, results)
# }

find_significant_corrs <- function(df, outcome_col, alpha = 0.05) {
  # Get numeric columns (excluding the outcome column)
  num_cols <- names(df)[sapply(df, is.numeric) & names(df) != outcome_col]
  
  # Test each numeric column
  results <- lapply(num_cols, function(col) {
    tryCatch({
      test <- cor.test(df[[col]], df[[outcome_col]], method = "spearman")
      data.frame(
        variable = col,
        p_value = test$p.value,
        rho = test$estimate,
        significant = test$p.value < alpha,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      # Handle cases with no variance or other errors
      data.frame(
        variable = col,
        p_value = NA,
        rho = NA,
        significant = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })
  
  # Combine into single dataframe
  do.call(rbind, results)
}

bcma_mpnn_dt = dt %>% filter(campaign=="BCMA_l59_MPNN")
bcma_mpnn_dt_filtered <- bcma_mpnn_dt[, !grepl("^X|^max|CD22|CD19|BCMA", names(bcma_mpnn_dt))]
#bcma_mpnn_dt_filtered_diffs = find_significant_diffs(bcma_mpnn_dt_filtered,"high_tonic")
bcma_mpnn_dt_filtered_diffs = find_significant_corrs(bcma_mpnn_dt_filtered,"coculture.CAR.Only")

bcma_mpnn_dt_filtered_diffs = bcma_mpnn_dt_filtered_diffs %>% mutate(target="BCMA")
bcma_mpnn_dt_filtered_diffs %>% arrange(p_value)

cd19_l112_mpnn_dt = dt %>% filter(campaign=="CD19_l112_MPNN")
cd19_l112_mpnn_dt_filtered <- cd19_l112_mpnn_dt[, !grepl("^X|^max|CD22|CD19|BCMA", names(cd19_l112_mpnn_dt))]
#cd19_l112_mpnn_dt_filtered_diffs = find_significant_diffs(cd19_l112_mpnn_dt_filtered,"high_tonic")
cd19_l112_mpnn_dt_filtered_diffs = find_significant_corrs(cd19_l112_mpnn_dt_filtered,"coculture.CAR.Only")

cd19_l112_mpnn_dt_filtered_diffs = cd19_l112_mpnn_dt_filtered_diffs %>% mutate(target="CD19_l112")
cd19_l112_mpnn_dt_filtered_diffs %>% arrange(p_value)

cd19_l115_mpnn_dt = dt %>% filter(campaign=="CD19_l115_MPNN")
cd19_l115_mpnn_dt_filtered <- cd19_l115_mpnn_dt[, !grepl("^X|^max|CD22|CD19|BCMA", names(cd19_l115_mpnn_dt))]
cd19_l115_mpnn_dt_filtered_diffs = find_significant_corrs(cd19_l115_mpnn_dt_filtered,"coculture.CAR.Only")
cd19_l115_mpnn_dt_filtered_diffs = cd19_l115_mpnn_dt_filtered_diffs %>% mutate(target="CD19_l115")
cd19_l115_mpnn_dt_filtered_diffs %>% arrange(p_value)

cd22_mpnn_dt = dt %>% filter(campaign=="CD22_l61_MPNN")
cd22_mpnn_dt_filtered <- cd22_mpnn_dt[, !grepl("^X|^max|CD22|CD19|BCMA", names(cd22_mpnn_dt))]
cd22_mpnn_dt_filtered_diffs = find_significant_corrs(cd22_mpnn_dt_filtered,"coculture.CAR.Only")
cd22_mpnn_dt_filtered_diffs = cd22_mpnn_dt_filtered_diffs %>% mutate(target="CD22")
cd22_mpnn_dt_filtered_diffs %>% arrange(p_value)

combined_diffs = rbind(bcma_mpnn_dt_filtered_diffs,cd19_l112_mpnn_dt_filtered_diffs,cd19_l115_mpnn_dt_filtered_diffs,cd22_mpnn_dt_filtered_diffs)

combined_diffs_modified <- combined_diffs %>%
  # 1. Calculate Fold Change (Mean_True / Mean_False)
  # Adding a small epsilon to the denominator to prevent division by zero if mean_false is 0
  mutate(
    spearman = rho,
    #fold_change = mean_true / (mean_false + 1e-9),
    # If mean_false is 0 and mean_true is non-zero, fold_change will be very large.
    # Consider replacing with log2(mean_true / mean_false) for symmetric colors if appropriate,
    # but based on the prompt, we use the raw ratio.
    
    # 2. Define significance label (Star annotation)
    star_label = case_when(
      significant == TRUE ~ "*",
      is.na(significant) ~ "", # No star for NA/No Test
      TRUE ~ "" # No star for Not Significant
    ),
    
    # Simplify variable names for better display
    variable_short = gsub("Average_|_", " ", variable) %>% trimws()
  )

# --- Plotting ---

# Define a diverging color palette for the fold change (e.g., blue to red)
# You may want to use a colorblind-safe palette like 'RdBu' or 'Spectral'
fold_change_colors <- c("#2166AC", "white", "#B2182B") 

# Create the heatmap
ggplot(combined_diffs_modified, aes(x = factor(variable_short), y = factor(target))) +
  # 1. Tile layer: color based on the continuous 'fold_change'
  geom_tile(aes(fill = spearman), color = "black", linewidth = 0.5) +
  
  # 2. Annotation layer: Display the formatted fold change value
  geom_text(
    aes(label = sprintf("%.2f", spearman)), # Format to 2 decimal places
    color = "black", 
    size = 2.5, # Reduced size for clarity
    vjust = 0.5 # Center the text vertically
  ) +
  
  # 2. Text layer: Add stars for significant results
  geom_text(aes(label = star_label), color = "black", size = 4, fontface = "bold") +
  
  # # 3. Scale Fill: Use a continuous, diverging color scale for fold change
  # scale_fill_gradient2(
  #   low = fold_change_colors[1],   # Blue for low FC
  #   mid = fold_change_colors[2],   # White for FC near 1
  #   high = fold_change_colors[3],  # Red for high FC
  #   midpoint = 1,                  # Midpoint is 1 (no change)
  #   name = "Fold Change (Tonic/Non-Tonic)"
  # ) +
  
  # --- Theme and Labels (Kept from original) ---
  labs(
    x = "Feature Variable",
    y = "Target"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
    # Removed legend.position = "none" to display the fold change legend
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_b()

## Calculate AUROC and AUPRC of iptm/pae/ipSAE


target_cols = c("pred_1", "pred_2")

calculate_metrics <- function(data, label_col, pred_cols) {
  y <- as.numeric(as.logical(data[[label_col]]))
  results <- lapply(pred_cols, function(col) {
    scores <- data[[col]]
    
    # PRROC expects scores separated by class
    pos <- scores[y == 1]
    neg <- scores[y == 0]
    
    # Calculate metrics
    roc <- PRROC::roc.curve(scores.class0 = pos, scores.class1 = neg)
    pr  <- PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg)
    
    data.frame(
      predictor = col,
      auroc     = roc$auc,
      auprc     = pr$auc.integral
    )
  })
  do.call(rbind, results)
}

dt_ysd$Average_pae_interaction_inv = 1/dt_ysd$Average_pae_interaction
dt_ysd$Average_pae_all_inv = 1/dt_ysd$Average_pae_all
dt_ysd$Average_complex_ipde_inv = 1/dt_ysd$Average_complex_ipde
dt_ysd$Average_complex_pde_inv = 1/dt_ysd$Average_complex_pde
dt_ysd$Average_dG_neg = -dt_ysd$Average_dG


dl_metrics = c(
  "Average_ipSAE_max", "Average_ipSAE_min","Average_iptm","Average_pae_interaction_inv",
  "Average_i_pLDDT","Average_complex_plddt","Average_complex_ipde_inv",
  "Average_pDockQ","Average_pDockQ2","Average_LIS","Average_complex_pde_inv","Average_pae_all_inv","Average_ss_pLDDT"
)
biophysical_metrics = c(
  "Average_InterfaceHbondsPercentage","Average_Interface_Hydrophobicity","Average_dSASA","Average_Unrelaxed_Clashes",
  "Average_Interface_SASA_.","Average_Surface_Hydrophobicity","Average_dG_neg",
  "Average_n_InterfaceHbonds","Average_ShapeComplementarity"
)

metric_cols = c(dl_metrics,biophysical_metrics)


metrics_summary_overall_1000nM = calculate_metrics(
  data = dt_ysd, 
  label_col = "binder_by_YSD_1000nM", 
  pred_cols = metric_cols
)
metrics_summary_overall_1000nM = metrics_summary_overall_1000nM %>% mutate(
  metric_type = case_when(
    predictor %in% dl_metrics ~ "ML-based metric",
    predictor %in% biophysical_metrics ~ "Physics-based metric",
  )
)
metrics_summary_overall_100nM = calculate_metrics(
  data = dt_ysd, 
  label_col = "binder_by_YSD_100nM", 
  pred_cols = metric_cols
)
metrics_summary_overall_100nM = metrics_summary_overall_100nM %>% mutate(
  metric_type = case_when(
    predictor %in% dl_metrics ~ "ML-based metric",
    predictor %in% biophysical_metrics ~ "Physics-based metric",
  )
)


overall_auroc_plot_1000nM = ggplot(metrics_summary_overall_1000nM, aes(x = auroc, y = reorder(predictor, auroc),fill=metric_type)) +
  geom_col() + pretty_plot() + L_border() +
  scale_x_continuous(expand=c(0,Inf)) +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  scale_fill_manual(values=c("ML-based metric"="dodgerblue3","Physics-based metric"="firebrick")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "none",axis.text = element_text(size=5))
overall_auroc_plot_1000nM

cowplot::ggsave2("../plots/overall_ysd_auroc_plot_1000nM.pdf",overall_auroc_plot_1000nM,dpi=300,width=3,height=1.8,unit="in")

overall_auroc_plot_100nM = ggplot(metrics_summary_overall_100nM, aes(x = auroc, y = reorder(predictor, auroc),fill=metric_type)) +
  geom_col() + pretty_plot() + L_border() +
  scale_x_continuous(expand=c(0,Inf)) +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  scale_fill_manual(values=c("ML-based metric"="dodgerblue3","Physics-based metric"="firebrick")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "none",axis.text = element_text(size=5))
overall_auroc_plot_100nM
cowplot::ggsave2("../plots/overall_ysd_auroc_plot_100nM.pdf",overall_auroc_plot_100nM,dpi=300,width=3,height=1.8,unit="in")

# overall_auprc_plot = ggplot(metrics_summary_overall, aes(x = auprc, y = reorder(predictor, auprc),fill=metric_type)) +
#   geom_col() + pretty_plot() + L_border() +
#   scale_x_continuous(expand=c(0,Inf)) +
#   geom_vline(xintercept = 0.5, linetype="dashed") +
#   theme(axis.title.y = element_blank())


antigens_to_test = c("BCMA","CD19","CD22")
metrics_by_antigen_1000nM = lapply(antigens_to_test, function(ag) {
  dt_ysd %>%
    filter(antigen == ag) %>%
    calculate_metrics(label_col = "binder_by_YSD_1000nM", pred_cols = metric_cols) %>%
    mutate(antigen = ag) # Track which antigen this row belongs to
}) %>% bind_rows()

metrics_by_antigen_100nM = lapply(antigens_to_test, function(ag) {
  dt_ysd %>%
    filter(antigen == ag) %>%
    calculate_metrics(label_col = "binder_by_YSD_100nM", pred_cols = metric_cols) %>%
    mutate(antigen = ag) # Track which antigen this row belongs to
}) %>% bind_rows()

metrics_by_antigen_1000nM_auroc_plot = ggplot(metrics_by_antigen_1000nM, aes(x = auroc, y = reorder(predictor, auroc), color = antigen)) +
  geom_point(size = 1, alpha = 0.75) +
  scale_color_manual(values = c(
    "BCMA" = "firebrick",
    "CD19" = "dodgerblue3",
    "CD22" = "orange"
  )) +
  theme_minimal() +
  pretty_plot() + L_border() +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "none",axis.text = element_text(size=5))

cowplot::ggsave2("../plots/ysd_auroc_by_antigen_plot_1000nM.pdf",metrics_by_antigen_1000nM_auroc_plot,dpi=300,width=3,height=1.8,unit="in")

metrics_by_antigen_100nM_auroc_plot = ggplot(metrics_by_antigen_100nM, aes(x = auroc, y = reorder(predictor, auroc), color = antigen)) +
  geom_point(size = 1, alpha = 0.75) +
  scale_color_manual(values = c(
    "BCMA" = "firebrick",
    "CD19" = "dodgerblue3",
    "CD22" = "orange"
  )) +
  theme_minimal() +
  pretty_plot() + L_border() +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "none",axis.text = element_text(size=5))
cowplot::ggsave2("../plots/ysd_auroc_by_antigen_plot_100nM.pdf",metrics_by_antigen_100nM_auroc_plot,dpi=300,width=3,height=1.8,unit="in")


metrics_by_antigen_100nM_auroc_plot

# TRUTH_COLUMN <- "binder_by_YSD_1000nM"
# TARGET_CAMPAIGNS <- c(
#   "BCMA_E3_fold conditioned", 
#   "BCMA_BindCraft_Small", 
#   "CD19_BindCraft_big_1", 
#   "CD22_BindCraft_Domain 7_small"
# )
# 
# TRUTH_COLUMN <- "binder_by_estimated_KD_100nM"
# TARGET_CAMPAIGNS <- c(
#   "BCMA_l59_MPNN", 
#   "CD19_l112_MPNN", 
#   "CD19_l115_MPNN", 
#   "CD22_l61_MPNN"
# )
# 
# combined_auc_pr_df <- dt %>%
#   # 1. Filter campaigns and select relevant columns
#   filter(campaign %in% TARGET_CAMPAIGNS) %>%
#   dplyr::select(campaign, all_of(TRUTH_COLUMN), starts_with("Average")) %>%
#   dplyr::select(-Average_InterfaceAAs) %>%
#   
#   # 2. Prepare columns for analysis
#   mutate(
#     # Convert all Average columns to numeric
#     across(starts_with("Average"), as.numeric),
#     # Create the truth column: use the simple factor conversion.
#     truth_value = factor(!!sym(TRUTH_COLUMN), levels = c(FALSE, TRUE))
#   ) %>%
#   
#   # 3. Pivot predictor metrics to long format
#   pivot_longer(
#     cols = starts_with("Average"),
#     names_to = "feature",
#     values_to = "predictor_score"
#   ) %>%
#   
#   # 4. Group and filter for valid AUPRC calculation
#   group_by(campaign, feature) %>%
#   # Filter must be TRUE/FALSE *after* it's factored.
#   # We assume FALSE/TRUE are present in the filtered truth_value.
#   filter(n_distinct(truth_value) >= 2, !is.na(predictor_score)) %>%
#   
#   # 5. Calculate AUPRC for each group - FINAL, MOST ROBUST FIX
#   summarise(
#     .metric = "pr_auc",
#     .estimator = "binary",
#     # Pass the columns directly from the current group's context
#     .estimate = pr_auc_vec(
#       # The values are already filtered and numeric
#       truth = factor(truth_value),
#       estimate = predictor_score,
#       event_level = "second"
#     ),
#     .groups = "drop"
#   ) %>%
#   # Rename and select columns
#   rename(AUPRC_Value = .estimate) %>%
#   dplyr::select(feature, campaign, AUPRC_Value)
# 
# heatmap_plot <- combined_auc_pr_df %>%
#   # Ensure feature names are clean for better plotting
#   mutate(feature = gsub("Average_", "", feature)) %>%
#   
#   ggplot(aes(x = campaign, y = feature, fill = AUPRC_Value)) +
#   geom_tile(color = "white") + # Add tile boundaries for separation
#   
#   # Add the AUPRC value as text on each tile
#   geom_text(aes(label = sprintf("%.2f", AUPRC_Value)), color = "black", size = 3) +
#   
#   # Define the color scale (using viridis for better visual distinction)
#   scale_fill_gradient(
#     low = "white", 
#     high = "darkred", 
#     name = "AUPRC Value",
#     limits = c(0, 1) # AUPRC values range from 0 to 1
#   ) +
#   
#   # Customize labels and theme
#   labs(
#     title = "AUPRC Heatmap: Feature Performance by Campaign",
#     x = "Campaign",
#     y = "Predictor Feature"
#   ) +
#   # Rotate x-axis labels for readability and apply a clean theme
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title.x = element_blank(), # Campaign title is redundant with text
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# # Display the plot
# print(heatmap_plot)

## Plot Percent of Binders above threshold
## Get success rate  of CAR staining based campaigns
iptm_success_rate_mpnn = dt %>% group_by(antigen,above_iptm_85) %>%
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) %>%
  mutate(metric="iptm>=0.85") %>% filter(above_iptm_85)
pae_success_rate_mpnn = dt %>% group_by(antigen,below_pae_10) %>%
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) %>%
  mutate(metric="pae_interaction<10") %>% filter(below_pae_10)
ipSAE_success_rate_mpnn = dt %>% group_by(antigen,above_ipSAE_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) %>%
  mutate(metric="ipSAE>=0.85") %>% filter(above_ipSAE_85)
# pde_success_rate_mpnn = dt %>% group_by(antigen,below_pde_1) %>% 
#   summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) %>%
#   mutate(metric="pde<1") %>% filter(below_pde_1)
# ipde_success_rate_mpnn = dt %>% group_by(antigen,below_ipde_1) %>% 
#   summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) %>%
#   mutate(metric="ipde<1") %>% filter(below_ipde_1)

# campaign_success_rates_mpnn_combined = rbind(iptm_success_rate_mpnn,pae_success_rate_mpnn,ipSAE_success_rate_mpnn,pde_success_rate_mpnn,ipde_success_rate_mpnn)
# campaign_success_rates_mpnn_combined$metric = factor(
#   campaign_success_rates_mpnn_combined$metric, levels = c("iptm>=0.85","pae_interaction<10","ipSAE>=0.85","pde<1","ipde<1")
# )

campaign_success_rates_mpnn_combined = rbind(iptm_success_rate_mpnn,pae_success_rate_mpnn,ipSAE_success_rate_mpnn)
campaign_success_rates_mpnn_combined$metric = factor(
  campaign_success_rates_mpnn_combined$metric, levels = c("iptm>=0.85","pae_interaction<10","ipSAE>=0.85")
)

mpnn_metric_barplot = ggplot(campaign_success_rates_mpnn_combined ,aes(x=antigen,y=success_rate*100,fill=metric)) +
  geom_bar(stat="identity",position="dodge", color="black") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = c(
    "iptm>=0.85" = "dodgerblue3",
    "pae_interaction<10" = "orange",
    "ipSAE>=0.85" = "firebrick"
    #"pde<1"="green",
    #"ipde<1"="purple"
    ),name = "metric") +
  theme(
    #legend.position = "none",
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.line = element_line(colour = 'black', size = 0.5),
    ) +
  labs(x="Antigen",y="% Binders <1000nM")
mpnn_metric_barplot
cowplot::ggsave2("../plots/mpnn_metric_comparisons.pdf",mpnn_metric_barplot,height=1.6,width=1.8)


## Show percentage satisfying this metric per campaign
ipSAE_success_comparison = dt %>% group_by(antigen,above_ipSAE_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) 

ipSAE_min_success_comparison = dt %>% group_by(antigen,above_ipSAE_min_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE))

ipSAE_max_success_comparison = dt %>% group_by(antigen,above_ipSAE_max_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) 

ipSAE_a_to_b_success_comparison = dt %>% group_by(antigen,above_ipSAE_a_to_b_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE))

ipSAE_b_to_a_success_comparison = dt %>% group_by(antigen,above_ipSAE_b_to_a_85) %>% 
  summarise(success_rate=mean(binder_by_estimated_KD_1000nM==TRUE, na.rm=TRUE)) 

## Fisher test to see if the proportions are different
ipSAE_fisher_results = dt %>%
  group_by(antigen) %>%
  summarise(fisher_pval = fisher.test(table(above_ipSAE_85, binder_by_estimated_KD_1000nM))$p.value)

ipSAE_success_rate_plot = ggplot(ipSAE_success_comparison, aes(x = antigen, y = success_rate * 100, fill = above_ipSAE_85)) +
  geom_bar(stat = "identity", position = "dodge",color="black") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values=c("gray","firebrick3")) +
  theme(legend.position = "none") + labs(x="Antigen",y="% Binders <1000nM")

ipSAE_success_rate_plot

## The fisher pvalues are added in illustrator
ipSAE_fisher_results
cowplot::ggsave2("../plots/mpnn_ipSAE_prop_test.pdf",ipSAE_success_rate_plot,height=1.6,width=1.8)

ipSAE_success_rate_plot_short = ipSAE_success_rate_plot +
  theme(axis.text.x = element_blank(),axis.title = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
ipSAE_success_rate_plot_short
cowplot::ggsave2("../plots/mpnn_ipSAE_prop_small.pdf",ipSAE_success_rate_plot_short,height=0.4,width=1.6)

## Do fisher for YSD
ysd_short_filtered = dt %>% 
  filter(campaign_short %in% c("BC 1 (BCMA)", "RFD 5 (BCMA)","BC 2 (CD19)","BC 1 (CD22)"))
ipSAE_ysd_success_comparison = ysd_short_filtered %>%
  mutate(campaign_short = factor(ysd_short_filtered$campaign_short,levels=c("BC 1 (BCMA)", "RFD 5 (BCMA)","BC 2 (CD19)","BC 1 (CD22)"))) %>%
  group_by(campaign_short,above_ipSAE_85) %>% 
  summarise(success_rate=mean(binder_by_YSD_1000nM==TRUE, na.rm=TRUE))

ipSAE_ysd_success_rate_plot = ggplot(ipSAE_ysd_success_comparison, aes(x = campaign_short, y = success_rate * 100, fill = above_ipSAE_85)) +
  geom_bar(stat = "identity", position = "dodge",color="black") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values=c("gray","firebrick3")) +
  theme(legend.position = "none") + labs(x="Campaign",y="% Enriched (<1000nM)")

ipSAE_ysd_success_rate_plot_short = ipSAE_ysd_success_rate_plot +
  theme(axis.text.x = element_blank(),axis.title = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
ipSAE_ysd_success_rate_plot_short

cowplot::ggsave2("../plots/ysd_ipSAE_prop_small.pdf",ipSAE_ysd_success_rate_plot_short,height=0.4,width=1.6)

ipSAE_ysd_fisher_results = dt %>%
  filter(campaign_short %in% c("BC 1 (BCMA)","BC 2 (CD19)","BC 1 (CD22)", "RFD 5 (BCMA)")) %>%
  group_by(campaign_short) %>%
  summarise(fisher_pval = fisher.test(table(above_ipSAE_85, binder_by_YSD_1000nM))$p.value)

ipSAE_ysd_fisher_results

# ## Check the correlation between Kd and CAR activity gain
# ## Show BCMA and CD22 since those were actual coculture hits
mpnn_bcma = dt %>% filter(campaign == "BCMA_l59_MPNN", !is.na(overall_activity_gain))
bcma_stain_pred <- predict(lm(overall_activity_gain ~ staining_1_nM, mpnn_bcma),
                   se.fit = TRUE, interval = "confidence")
bcma_stain_lims <- as.data.frame(bcma_stain_pred$fit)
bcma_stain_kd_plot <- ggplot(mpnn_bcma, aes(x = staining_1_nM, y = overall_activity_gain, color = binder_evolution_type)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Staining 1nM BCMA", y = "Activity Gain") +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
  geom_line(aes(x = staining_1_nM, y = bcma_stain_lims$lwr),
            linetype = 2, color = "grey") +
  geom_line(aes(x = staining_1_nM, y = bcma_stain_lims$upr),
            linetype = 2, color = "grey") +
  scale_color_manual(values = c(
    "interface" = "dodgerblue3",
    "non-interface" = "firebrick",
    "parental" = "black"
  ),name = "evolution") +
  #scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
  theme(legend.position = "none")#+ ylim(0,90)

bcma_stain_kd_plot
bcma_stain_kd_cor <- cor.test(mpnn_bcma$staining_1_nM, mpnn_bcma$overall_activity_gain, method = "spearman")
bcma_stain_kd_cor

## Do it for CD22
mpnn_cd22 = dt %>% filter(campaign == "CD22_l61_MPNN", !is.na(overall_activity_gain))
cd22_stain_pred <- predict(lm(overall_activity_gain ~ staining_0.1_nM, mpnn_cd22),
                        se.fit = TRUE, interval = "confidence")
cd22_stain_lims <- as.data.frame(cd22_stain_pred$fit)
cd22_stain_kd_plot <- ggplot(mpnn_cd22, aes(x = staining_0.1_nM, y = overall_activity_gain, color = binder_evolution_type)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Staining 0.1nM CD22", y = "Activity Gain") +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
  geom_line(aes(x = staining_0.1_nM, y = cd22_stain_lims$lwr),
            linetype = 2, color = "grey") +
  geom_line(aes(x = staining_0.1_nM, y = cd22_stain_lims$upr),
            linetype = 2, color = "grey") +
  scale_color_manual(values = c(
    "interface" = "dodgerblue3",
    "non-interface" = "firebrick",
    "parental" = "black"
  ),name = "evolution") +
  #scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
  theme(legend.position = "none")

cd22_stain_kd_plot
cd22_stain_kd_cor <- cor.test(mpnn_cd22$staining_0.1_nM, mpnn_cd22$overall_activity_gain, method = "spearman")
cd22_stain_kd_cor

bcma_stain_kd_plot | cd22_stain_kd_plot
print(bcma_stain_kd_cor)
print(cd22_stain_kd_cor)

bcma_cd22_stain_kd_plot = cowplot::plot_grid(bcma_stain_kd_plot,cd22_stain_kd_plot,ncol=2)
bcma_cd22_stain_kd_plot
cowplot::ggsave2("../plots/bcma_cd22_stain_activity_gain.pdf",bcma_cd22_stain_kd_plot,height=1.6,width=3.6)

## Show ipSAE distribution of YSD campaigns that had binders

# ysd_dt = dt %>% filter(campaign %in% c("BCMA_E3_fold conditioned","BCMA_BindCraft_Small","CD19_BindCraft_big_1","CD22_BindCraft_Domain 7_small"))
# ysd_dt = ysd_dt %>% mutate(
#   campaign_renamed = case_when(
#     campaign == "BCMA_E3_fold conditioned" ~ "BCMA\n(RFD)",
#     campaign == "BCMA_BindCraft_Small" ~ "BCMA\n(BC)",
#     campaign == "CD19_BindCraft_big_1" ~ "CD19\n(BC)",
#     campaign == "CD22_BindCraft_Domain 7_small" ~ "CD22\n(BC)",
#   )
# )
ysd_dt = dt %>% filter(campaign %in% c("BCMA RFDiffusion Campaign 5","BCMA BindCraft Campaign 1","CD19 BindCraft Campaign 2","CD22 BindCraft Campaign 1"))

ysd_campaign_ipsae_plot = ggplot(ysd_dt,aes(x=campaign,y=Average_ipSAE,color=binder_by_YSD_1000nM)) +
  #geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipSAE") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_ipsae_plot
cowplot::ggsave2("../plots/ysd_ipsae_cutoff.pdf",ysd_campaign_ipsae_plot,height=1.6,width=1.8)

ysd_campaign_ipsae_plot

## Add ipSAE max/min/a->b/b->a

ysd_campaign_ipsae_max_plot = ggplot(ysd_dt,aes(x=campaign,y=Average_ipSAE_max,color=binder_by_YSD_1000nM)) +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipSAE max") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_ipsae_max_plot
cowplot::ggsave2("../plots/ysd_ipsae_max.pdf",ysd_campaign_ipsae_max_plot,height=1.5,width=1.75)

ysd_campaign_ipsae_min_plot = ggplot(ysd_dt,aes(x=campaign,y=Average_ipSAE_min,color=binder_by_YSD_1000nM)) +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  #geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipSAE min") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_ipsae_min_plot
cowplot::ggsave2("../plots/ysd_ipsae_min.pdf",ysd_campaign_ipsae_min_plot,height=1.5,width=1.75)

ysd_campaign_ipsae_a_to_b_plot = ggplot(ysd_dt,aes(x=campaign,y=Average_Chn1_A_to_Chn2_B_ipSAE,color=binder_by_YSD_1000nM)) +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  #geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipSAE (Chn1:binder -> Chn2:target)") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_ipsae_a_to_b_plot
cowplot::ggsave2("../plots/ysd_ipsae_a_to_b.pdf",ysd_campaign_ipsae_a_to_b_plot,height=1.5,width=1.75)

ysd_campaign_ipsae_b_to_a_plot = ggplot(ysd_dt,aes(x=campaign,y=Average_Chn1_B_to_Chn2_A_ipSAE,color=binder_by_YSD_1000nM)) +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  #geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipSAE (Chn1:target -> Chn2:binder)") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_ipsae_b_to_a_plot
cowplot::ggsave2("../plots/ysd_ipsae_b_to_a.pdf",ysd_campaign_ipsae_b_to_a_plot,height=1.5,width=1.75)


ysd_campaign_pae_plot = ggplot(ysd_dt,aes(x=campaign_renamed,y=Average_pae_interaction,color=binder_by_YSD_1000nM)) +
  #geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  geom_hline(yintercept = 10,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipAE") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_pae_plot
cowplot::ggsave2("../plots/ysd_ipsae_cutoff.pdf",ysd_campaign_ipsae_plot,height=1.6,width=1.8)

ysd_campaign_iptm_plot = ggplot(ysd_dt,aes(x=campaign_renamed,y=Average_iptm,color=binder_by_YSD_1000nM)) +
  geom_point(position = position_jitterdodge()) +
  pretty_plot() + L_border() +
  scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
  geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
  labs(x="Campaign",y="Average ipTM") +
  theme(legend.position = "none",text = element_text(size = 8))
ysd_campaign_iptm_plot
cowplot::ggsave2("../plots/ysd_iptm_cutoff.pdf",ysd_campaign_iptm_plot,height=1.6,width=1.8)

ysd_cutoff_plots_combined = cowplot::plot_grid(ysd_campaign_pae_plot,ysd_campaign_iptm_plot,ysd_campaign_ipsae_plot,ncol=3)
cowplot::ggsave2("../plots/ysd_cutoff_plots_combined.pdf",ysd_cutoff_plots_combined,height=1.6,width=5)



## Check plDDT across campaigns
# ysd_campaign_lddt_plot = ggplot(dt,aes(x=Average_monomer_lDDT,y=campaign,color=binder_by_YSD_1000nM)) +
#   geom_boxplot() +
#   geom_point(position = position_jitterdodge()) +
#   pretty_plot() + L_border() +
#   scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
#   #geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
#   labs(x="Campaign",y="Average ipSAE") #+
#   #theme(legend.position = "none")
# ysd_campaign_lddt_plot

# ggplot(dt,aes(x=Average_monomer_lDDT,y=Average_ipSAE,color=binder_by_YSD_1000nM)) +
#     geom_point() +
#     pretty_plot() + L_border() +
#     scale_color_manual(values=c("FALSE"="gray","TRUE"="dodgerblue3")) +
#     #geom_hline(yintercept = 0.85,linetype="dashed",color="black") +
#     labs(x="Campaign",y="Average ipSAE") #+
#     #theme(legend.position = "none")
# ysd_campaign_lddt_plot

## Get success rate of YSD based campaigns
iptm_ysd_success_rate = dt %>% group_by(campaign,above_iptm_85) %>%
  summarise(success_rate=mean(binder_by_YSD_100nM==TRUE)) %>%
  mutate(metric="iptm>=0.85") %>% filter(above_iptm_85)
pae_ysd_success_rate = dt %>% group_by(campaign,below_pae_10) %>%
  summarise(success_rate=mean(binder_by_YSD_100nM==TRUE)) %>%
  mutate(metric="pae_interaction<10") %>% filter(below_pae_10)
ipSAE_ysd_success_rate = dt %>% group_by(campaign,above_ipSAE_85) %>% 
  summarise(success_rate=mean(binder_by_YSD_100nM==TRUE)) %>%
  mutate(metric="ipSAE>=0.85") %>% filter(above_ipSAE_85)

# ipSAE_min_ysd_success_rate = dt %>% group_by(campaign,above_ipSAE_85) %>% 
#   summarise(success_rate=mean(binder_by_YSD_100nM==TRUE)) %>%
#   mutate(metric="ipSAE>=0.85") %>% filter(above_ipSAE_85)
# 
# ipSAE_max_ysd_success_rate = dt %>% group_by(campaign,above_ipSAE_85) %>% 
#   summarise(success_rate=mean(binder_by_YSD_100nM==TRUE)) %>%
#   mutate(metric="ipSAE>=0.85") %>% filter(above_ipSAE_85)

# pde_ysd_success_rate = dt %>% group_by(campaign,below_pde_1) %>% 
#   summarise(success_rate=mean(binder_by_YSD_1000nM==TRUE, na.rm=TRUE)) %>%
#   mutate(metric="pde<1") %>% filter(below_pde_1)
# ipde_ysd_success_rate = dt %>% group_by(campaign,below_ipde_1) %>% 
#   summarise(success_rate=mean(binder_by_YSD_1000nM==TRUE, na.rm=TRUE)) %>%
#   mutate(metric="ipde<1") %>% filter(below_ipde_1)

campaign_ysd_success_rates_combined = rbind(iptm_ysd_success_rate,pae_ysd_success_rate,ipSAE_ysd_success_rate)
                                        #pde_ysd_success_rate,ipde_ysd_success_rate)

campaign_ysd_success_rates_combined


campaign_ysd_success_rates_combined = campaign_ysd_success_rates_combined %>% mutate(
  campaign_renamed = case_when(
    campaign == "BCMA_E3_fold conditioned" ~ "BCMA\n(RFD)",
    campaign == "BCMA_BindCraft_Small" ~ "BCMA\n(BC)",
    campaign == "CD19_BindCraft_big_1" ~ "CD19\n(BC)",
    campaign == "CD22_BindCraft_Domain 7_small" ~ "CD22\n(BC)",
  ),
  metric = factor(metric,levels=c("iptm>=0.85","pae_interaction<10","ipSAE>=0.85"))
)
ysd_mpnn_metric_barplot = ggplot(campaign_ysd_success_rates_combined %>% filter(success_rate > 0),aes(x=campaign_renamed,y=success_rate*100,fill=metric)) +
  geom_bar(stat="identity",position="dodge", color="black") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = c(
    "iptm>=0.85" = "dodgerblue3",
    "pae_interaction<10" = "orange",
    "ipSAE>=0.85" = "firebrick"
    #"pde<1"="green",
    #"ipde<1"="purple"
  ),name = "metric") +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.line = element_line(colour = 'black', size = 0.5),
  ) +
  labs(x="Antigen",y="% YSD Binders <1000nM") + 
  #ylim(0,35) +
  scale_y_continuous (expand = c(0,0))
ysd_mpnn_metric_barplot
cowplot::ggsave2("../plots/ysd_mpnn_metric_comparisons.pdf",ysd_mpnn_metric_barplot,height=1.6,width=1.6)

## Check the Kd metric relationship
## ipSAe -> Kd

## MPNN Campaign PR Curves
mpnn_dt = dt %>%
  filter(campaign %in% c("BCMA_l59_MPNN","CD19_l112_MPNN","CD19_l115_MPNN","CD22_l61_MPNN"))

mpnn_campaign_colors = c(
  "BCMA_l59_MPNN" = "firebrick",
  "CD19_l112_MPNN" = "dodgerblue3",
  "CD19_l115_MPNN" = "purple",
  "CD22_l61_MPNN" = "orange"
)

plddt_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_complex_plddt,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average complex pLDDT", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

pae_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_pae_interaction,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average pAE interaction", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

iptm_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_iptm,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipTM", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

ipsae_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_ipSAE,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipSAE", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

ipsae_kd

ipsae_min_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_ipSAE_min,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipSAE min", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

ipsae_max_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_ipSAE_max,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipSAE max", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

ipsae_chn1_a_to_chn2_b_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_Chn1_A_to_Chn2_B_ipSAE,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipSAE (Chn1:binder -> Chn2:target)", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

ipsae_chn1_b_to_chn1_a_kd <- mpnn_dt %>% 
  ggplot(aes(x = Average_Chn1_B_to_Chn2_A_ipSAE,y=neg_log10_Kd_M, color = campaign)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=mpnn_campaign_colors) +
  theme(legend.position = "none") +
  labs(x="Average ipSAE (Chn1:target -> Chn2:binder)", y="Estimated Kd (M)") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^-", x)))

combined_ipsae_kd_plots = cowplot::plot_grid(ipsae_min_kd,ipsae_max_kd,ipsae_chn1_a_to_chn2_b_kd,ipsae_chn1_b_to_chn1_a_kd,nrow=2)
cowplot::ggsave2("../plots/combined_ipsae_kd_plots.pdf",combined_ipsae_kd_plots,dpi=300,height=3.2,width=3.5)

ipsae_comparison_plot = ggplot(mpnn_dt,aes(y=Average_ipSAE_max,x=Average_Chn1_B_to_Chn2_A_ipSAE,color=campaign)) +
  geom_point() + 
  geom_abline(slope=1, intercept=0,linetype="dashed") +
  pretty_plot() + L_border() +
  theme(legend.position = "none") +
  scale_color_manual(values=mpnn_campaign_colors)
ipsae_comparison_plot

cowplot::ggsave2("../plots/ipsae_comparison_plot.pdf",ipsae_comparison_plot,dpi=300,height=1.5,width=1.75)


ipsae_min_comparison_plot = ggplot(mpnn_dt,aes(y=Average_ipSAE_min,x=Average_Chn1_B_to_Chn2_A_ipSAE,color=campaign)) +
  geom_point() + 
  geom_abline(slope=1, intercept=0,linetype="dashed") +
  pretty_plot() + L_border() +
  theme(legend.position = "none") +
  scale_color_manual(values=mpnn_campaign_colors)
ipsae_min_comparison_plot

cowplot::ggsave2("../plots/ipsae_min_comparison_plot.pdf",ipsae_min_comparison_plot,dpi=300,height=1.5,width=1.75)


combined_metric_kd = cowplot::plot_grid(pae_kd,iptm_kd,ipsae_kd,ncol=3)

cowplot::ggsave2("../plots/combined_metric_kd_plot.pdf",combined_metric_kd,dpi=300,height=1.6,width=5)

ipsae_kd_zoomed = ipsae_kd + 
  xlim(0.80,0.9) +
  geom_vline(xintercept = 0.85, linetype="dashed",color="gray")
ipsae_kd_zoomed
cowplot::ggsave2("../plots/ipSAE_kd_plot_zoomed.pdf",ipsae_kd_zoomed,dpi=300,height=1.6,width=1.6)


ipsae_max_binder_kd_zoomed = ipsae_kd + 
  xlim(0.80,0.9) +
  geom_vline(xintercept = 0.85, linetype="dashed",color="gray") +
  theme(legend.position = "none")
ipsae_max_binder_kd_zoomed
cowplot::ggsave2("../plots/ipSAE_max_binder_kd_plot_zoomed.pdf",ipsae_max_binder_kd_zoomed,dpi=300,height=1.6,width=1.75)

ipsae_target_binder_kd_zoomed = ipsae_chn2_a_to_chn1_b_kd + 
  xlim(0.80,0.9) +
  geom_vline(xintercept = 0.85, linetype="dashed",color="gray") +
  theme(legend.position = "none")
ipsae_target_binder_kd_zoomed
cowplot::ggsave2("../plots/ipSAE_target_binder_kd_plot_zoomed.pdf",ipsae_target_binder_kd_zoomed,dpi=300,height=1.6,width=1.75)



## plot charges
ysd_campaigns = c("BCMA_E3_fold conditioned","BCMA_BindCraft_Small","CD19_BindCraft_big_1","CD22_BindCraft_Domain 7_small")
mpnn_campaigns = c("BCMA_l59_MPNN","CD19_l112_MPNN","CD19_l115_MPNN","CD22_l61_MPNN")
dt_subset = dt %>% filter(campaign %in% c(ysd_campaigns,mpnn_campaigns))
#dt_subset = dt %>% filter(!is.na(binder_by_estimated_KD_1000nM)|(!is.na(binder_by_YSD_1000nM)))
ggplot(dt_subset,aes(x=binder_seq_net_charge,fill=campaign)) +
  geom_histogram() + pretty_plot() + L_border()

dt_subset_mpnn = dt %>% filter(campaign %in% c(mpnn_campaigns))

ggplot(dt_subset_mpnn,aes(x=binder_seq_net_charge,y=coculture.CAR.Only)) +
  geom_point() + pretty_plot() + L_border() + facet_wrap(~campaign)

dt_subset_bc_bcma = dt %>% filter(campaign %in% c("BCMA_BindCraft_Small"))

ggplot(dt_subset_bc_bcma,aes(x=binder_seq_net_charge,y=coculture.CAR.Only)) +
  geom_point()# + pretty_plot() + L_border()

  #labs(x = "D1 Average ipSAE", y = "D1.N0 Average ipSAE") +
  #scale_color_manual(values = c("black","firebrick")) +
  #geom_abline(slope=1, intercept=0, color = "gray", linetype="dashed") +
  #theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank()) + ylim(0,1) + xlim(0,1)

p1

ggplot(mpnn_dt %>% filter(kd_fit_quality == "Good (Curve Fit)", neg_log10_Kd_M >= 6),aes(x=Average_ipSAE, y=neg_log10_Kd_M, color=campaign)) +
  geom_point() +
  facet_wrap(~campaign)

ggplot(mpnn_dt %>% filter(kd_fit_quality == "Good (Curve Fit)", neg_log10_Kd_M >= 6),aes(x=Average_ipSAE, y=overall_activity_gain, color=campaign)) +
  geom_point() +
  facet_wrap(~campaign)



# mpnn_cd19_l112 = dt %>% filter(campaign=="CD19_l112_MPNN", !is.na(overall_activity_gain))
# cd19_l112_stain_kd_plot <- ggplot(mpnn_cd19_l112, aes(x = staining_1_nM, y = overall_activity_gain, color = binder_evolution_type)) + 
#   geom_point() + pretty_plot(fontsize = 8) + L_border() + 
#   labs(x = "Staining 1nM CD19", y = "Activity Gain") +
#   geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
#   # geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$lwr), 
#   #           linetype = 2, color = "grey") +
#   # geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$upr), 
#   #           linetype = 2, color = "grey") +
#   scale_color_manual(values = c(
#     "interface" = "dodgerblue3",
#     "non-interface" = "firebrick",
#     "parental" = "black"
#   ),name = "evolution") +
#   #scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
#   theme(legend.position = "none")
# cd19_l112_stain_kd_plot
# 
# cd19_l112_stain_kd_plot
# cd19_l112_stain_kd_cor <- cor.test(mpnn_cd19_l112$staining_1_nM, mpnn_cd19_l112$overall_activity_gain, method = "spearman")
# cd19_l112_stain_kd_cor
# 
# mpnn_cd19_l115 = dt %>% filter(campaign=="CD19_l115_MPNN", !is.na(overall_activity_gain))
# cd19_l115_stain_kd_plot <- ggplot(mpnn_cd19_l115, aes(x = staining_1_nM, y = overall_activity_gain, color = binder_evolution_type)) + 
#   geom_point() + pretty_plot(fontsize = 8) + L_border() + 
#   labs(x = "Staining 1nM CD19", y = "Activity Gain") +
#   geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
#   # geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$lwr), 
#   #           linetype = 2, color = "grey") +
#   # geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$upr), 
#   #           linetype = 2, color = "grey") +
#   scale_color_manual(values = c(
#     "interface" = "dodgerblue3",
#     "non-interface" = "firebrick",
#     "parental" = "black"
#   ),name = "evolution") +
#   #scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
#   theme(legend.position = "none")
# cd19_l115_stain_kd_plot
# 
# cd19_l115_stain_kd_plot
# cd19_l115_stain_kd_cor <- cor.test(mpnn_cd19_l115$staining_1_nM, mpnn_cd19_l115$overall_activity_gain, method = "spearman")
# cd19_l115_stain_kd_cor
# 
# cd19_l112_stain_kd_plot | cd19_l115_stain_kd_plot
# ggplot(dt %>% filter(MPNN_campaign,!is.na(binder_by_estimated_KD_1000nM)),aes(x=Average_ipSAE,fill=binder_by_estimated_KD_1000nM)) + 
#   geom_histogram() + facet_wrap(~campaign,ncol=1)


# generate_pr_curve_df = function(data, campaign_name, target_col) {
#   feature_cols = c("Average_ipSAE","Average_pae_interaction","Average_iptm","Average_complex_plddt")
#   filtered_dt = data %>% filter(campaign == campaign_name)
#   pr_dfs = lapply(feature_cols, function(feature) {
#     pr_curve(filtered_dt, target_col, sym(feature)) %>%
#       mutate(feature = feature, campaign = campaign_name)
#   })
#   pr_aucs = lapply(feature_cols, function(feature) {
#     pr_auc(filtered_dt, target_col, sym(feature)) %>%
#       mutate(feature = feature, campaign = campaign_name)
#   })
#   top_k_prec = map_dfr(c(5, 10, 20), function(k) {
#     filtered_dt %>%
#       arrange(desc(sym(feature))) %>%
#       slice_head(n = k) %>%
#       summarise(
#         precision = sum(target_col == levels(target_col)[2]) / k,
#         metric = paste0("Precision@", k),
#         feature = feature,
#         campaign = campaign_name
#       )
#   })
#   pr_df = bind_rows(pr_dfs)
#   pr_auc =bind_rows(pr_aucs)
#   return(list(pr_df=pr_df,pr_auc=pr_auc,top_k_prec=top_k_prec))
# }

generate_pr_curve_df = function(data, campaign_name, target_col) {
  feature_cols = c("Average_ipSAE", "Average_pae_interaction", "Average_iptm", "Average_complex_plddt")
  
  # Prepare data: filter and ensure target is a factor for yardstick
  filtered_dt = data %>%
    filter(campaign == campaign_name) %>%
    mutate(!!sym(target_col) := as.factor(!!sym(target_col)))
  
  # Lists to store results
  pr_dfs = list()
  pr_aucs = list()
  top_k_precs = list()
  
  # Loop over features to calculate all metrics
  for (feature in feature_cols) {
    feature_sym = sym(feature)
    target_sym = sym(target_col)
    
    # Calculate PR curve points
    pr_dfs[[feature]] = pr_curve(filtered_dt, !!target_sym, !!feature_sym) %>%
      mutate(feature = feature, campaign = campaign_name)
    
    # Calculate AUPRC
    pr_aucs[[feature]] = pr_auc(filtered_dt, !!target_sym, !!feature_sym) %>%
      mutate(feature = feature, campaign = campaign_name)
    
    # Calculate Precision at K
    top_k_precs[[feature]] = map_dfr(c(5, 10, 15, 20, 25), function(k) {
      filtered_dt %>%
        # Arrange by the current feature in descending order
        arrange(desc(!!feature_sym)) %>%
        # Select top k rows
        slice_head(n = k) %>%
        # Calculate precision
        summarise(
          # Check for the positive class (assumes it's the second factor level)
          precision = sum(!!target_sym == levels(!!target_sym)[2]) / k,
          metric = paste0("Precision@", k),
          top_n = k,
          feature = feature,
          campaign = campaign_name
        )
    })
  }
  
  # Combine all results from the lists
  pr_df = bind_rows(pr_dfs)
  pr_auc = bind_rows(pr_aucs)
  top_k_prec = bind_rows(top_k_precs)
  
  # Return the list of data frames
  return(list(pr_df = pr_df, pr_auc = pr_auc, top_k_prec = top_k_prec))
}

pr_feature_col_map = c(
  "Average_iptm" = "dodgerblue3",
  "Average_pae_interaction" = "orange",
  "Average_ipSAE" = "firebrick",
  "Average_complex_plddt"="gray"
)

## YSD Campaigns PR Curves
bcma_e3_pr_res = generate_pr_curve_df(dt %>% mutate(Average_pae_interaction = -Average_pae_interaction),"BCMA_E3_fold conditioned","binder_by_YSD_100nM")
bcma_e3_pr_curve = ggplot(bcma_e3_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
bcma_e3_pr_curve
bcma_e3_pr_res$pr_auc
bcma_e3_pr_res$top_k_prec

bcma_bc_pr_res = generate_pr_curve_df(dt,"BCMA_BindCraft_Small","binder_by_YSD_1000nM")
bcma_bc_pr_curve = ggplot(bcma_bc_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
bcma_bc_pr_curve
bcma_bc_pr_res$pr_auc
bcma_e3_pr_res$top_k_prec

cd19_bc_pr_res = generate_pr_curve_df(dt,"CD19_BindCraft_big_1","binder_by_YSD_1000nM")
cd19_bc_pr_curve = ggplot(cd19_bc_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
cd19_bc_pr_curve
cd19_bc_pr_res$pr_auc
cd19_bc_pr_res$top_k_prec


cd22_d7_pr_res = generate_pr_curve_df(dt,"CD22_BindCraft_Domain 7_small","binder_by_YSD_1000nM")
cd22_d7_pr_curve = ggplot(cd22_d7_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
cd22_d7_pr_curve
cd22_d7_pr_res$pr_auc


bcma_mpnn_pr_res = generate_pr_curve_df(dt,"BCMA_l59_MPNN","binder_by_estimated_KD_100nM")
bcma_mpnn_pr_curve = ggplot(bcma_mpnn_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
bcma_mpnn_pr_curve
bcma_mpnn_pr_res$pr_auc
bcma_mpnn_pr_res$top_k_prec


cd19_l112_mpnn_pr_res = generate_pr_curve_df(dt,"CD19_l112_MPNN","binder_by_estimated_KD_100nM")
cd19_l112_mpnn_pr_curve = ggplot(cd19_l112_mpnn_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
cd19_l112_mpnn_pr_curve
cd19_l112_mpnn_pr_res$pr_auc
cd19_l112_mpnn_pr_res$top_k_prec


cd19_l115_mpnn_pr_res = generate_pr_curve_df(dt,"CD19_l115_MPNN","binder_by_estimated_KD_100nM")
cd19_l115_mpnn_pr_curve = ggplot(bcma_mpnn_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
cd19_l115_mpnn_pr_curve
cd19_l115_mpnn_pr_res$pr_auc
cd19_l115_mpnn_pr_res$top_k_prec


cd22_l61_mpnn_pr_res = generate_pr_curve_df(dt,"CD22_l61_MPNN","binder_by_estimated_KD_100nM")
cd22_l61_mpnn_pr_curve = ggplot(bcma_mpnn_pr_res$pr_df,aes(x=recall,y=precision,color=feature)) +
  geom_path() + coord_equal() + pretty_plot() + L_border() + 
  scale_color_manual(values=pr_feature_col_map) + theme(legend.position = "none") +
  labs(x="Recall",y="Precision")
cd22_l61_mpnn_pr_curve
cd22_l61_mpnn_pr_res$pr_auc
cd22_l61_mpnn_pr_res$top_k_prec


test_dt = dt %>% filter(campaign=="CD19_BindCraft_big_1")
test_dt = test_dt %>% mutate(
  pae_rank=rank(test_dt$Average_pae_interaction),
  iptm_rank=rank(-test_dt$Average_iptm)
)

t1 = ggplot(test_dt %>% arrange(binder_by_YSD_100nM),aes(x=pae_rank,y=Average_pae_interaction,color=binder_by_YSD_1000nM)) + 
  geom_point()

t2 = ggplot(test_dt %>% arrange(binder_by_YSD_100nM),aes(x=iptm_rank,y=Average_iptm,color=binder_by_YSD_1000nM)) + 
  geom_point()

t1 | t2




# t1 = ggplot(test_dt,aes(x=pae_rank,y=Average_pae_interaction,color=binder_by_estimated_KD_100nM)) + 
#   geom_point()
# 
# t2 = ggplot(test_dt,aes(x=iptm_rank,y=Average_iptm,color=binder_by_estimated_KD_100nM)) + 
#   geom_point()
# 
# t1 | t2


g1 = ggplot(mpnn_dt, aes(x=Average_pae_interaction, y=neg_log10_Kd_M, color = campaign)) +
  geom_point() + theme(legend.position = "none")

g2 = ggplot(mpnn_dt, aes(x=Average_ipSAE, y=neg_log10_Kd_M, color = campaign)) +
  geom_point() +
  scale_x_log10() +
  geom_vline(xintercept =  0.85, linestyle="dashed",color="gray")+ theme(legend.position = "none")

g3 = ggplot(mpnn_dt, aes(x=Average_complex_pde, y=neg_log10_Kd_M, color = campaign)) +
  geom_point()+ theme(legend.position = "none")

g4 = ggplot(mpnn_dt, aes(x=Average_complex_ipde, y=neg_log10_Kd_M, color = campaign)) +
  geom_point()+ theme(legend.position = "none")

cowplot::plot_grid(g1,g2,g3,g4,ncol=2)


g1 | g2 | g3 | g4

g3 = ggplot(mpnn_dt, aes(x=Average_complex_plddt, y=neg_log10_Kd_M, color = campaign)) +
  geom_point()

g3


g4 = ggplot(mpnn_dt, aes(x=Average_ipSAE, y=neg_log10_Kd_M, color = campaign)) +
  geom_point() +
  xlim(0.75,1) +
  #scale_x_log10() +
  geom_vline(xintercept =  0.85, linestyle="dashed",color="gray")
g4


## Kd -> Activtiy Gain
ggplot(mpnn_dt %>% filter(kd_fit_quality == "Good (Curve Fit)", neg_log10_Kd_M >= 6),aes(x=neg_log10_Kd_M, y=overall_activity_gain, color=campaign)) +
  geom_point() +
  facet_wrap(~campaign)

ggplot(mpnn_dt %>% filter(kd_fit_quality == "Good (Curve Fit)", neg_log10_Kd_M >= 6),aes(x=staining_1_nM, y=overall_activity_gain, color=campaign)) +
  geom_point() +
  facet_wrap(~campaign)

cd22_t1 = mpnn_dt %>% filter(kd_fit_quality == "Good (Curve Fit)", neg_log10_Kd_M >= 6, campaign == "CD22_l61_MPNN")
cor.test(cd22_t1$neg_log10_Kd_M,cd22_t1$overall_activity_gain,method = "spearman")

# %>%
  

mpnn_dt

ggplot(dt %>% filter(campaign %in% c("BCMA_E3_fold conditioned","BCMA_BindCraft_Small","CD19_BindCraft_big_1","CD22_BindCraft_Domain 7_small")),
       aes(x=campaign,y=Average_pae_interaction,fill=binder_by_YSD_100nM)) +
  geom_boxplot()

ggplot(dt %>% filter(campaign %in% c("BCMA_E3_fold conditioned","BCMA_BindCraft_Small","CD19_BindCraft_big_1","CD22_BindCraft_Domain 7_small")),
       aes(x=campaign,y=Average_pae_interaction,fill=binder_by_YSD_100nM)) +
  geom_boxplot()

# bcma_kd_pred <- predict(lm(overall_activity_gain ~ neg_log10_Kd_M, mpnn_bcma), 
#                    se.fit = TRUE, interval = "confidence")
# bcma_kd_lims <- as.data.frame(bcma_kd_pred$fit)
# 
# bcma_kd_plot <- ggplot(mpnn_bcma, aes(x = neg_log10_Kd_M, y = overall_activity_gain, color = binder_evolution_type)) + 
#   geom_point() + pretty_plot(fontsize = 8) + L_border() + 
#   labs(x = "BCMA Kd (M)", y = "Activity Gain") +
#   geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
#   geom_line(aes(x = neg_log10_Kd_M, y = bcma_kd_lims$lwr), 
#             linetype = 2, color = "grey") +
#   geom_line(aes(x = neg_log10_Kd_M, y = bcma_kd_lims$upr), 
#             linetype = 2, color = "grey") +
#   scale_color_manual(values = c(
#     "interface" = "dodgerblue3",
#     "non-interface" = "firebrick",
#     "parental" = "black"
#   ),name = "evolution") +
#   scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
#   theme(legend.position = "none")#+ ylim(0,90)
# bcma_kd_plot
# bcma_mpnn_kd_cor <- cor.test(mpnn_bcma$neg_log10_Kd_M, mpnn_bcma$overall_activity_gain, method = "spearman")
# bcma_mpnn_kd_cor
# ## Make a similar plot for CD22
# mpnn_cd22 = dt %>% filter(campaign == "CD22_l61_MPNN", !is.na(overall_activity_gain))
# cd22_kd_pred <- predict(lm(overall_activity_gain ~ neg_log10_Kd_M, mpnn_cd22), 
#                         se.fit = TRUE, interval = "confidence")
# cd22_kd_lims <- as.data.frame(cd22_kd_pred$fit)
# cd22_kd_plot <- ggplot(mpnn_cd22, aes(x = neg_log10_Kd_M, y = overall_activity_gain, color = binder_evolution_type)) + 
#   geom_point() + pretty_plot(fontsize = 8) + L_border() + 
#   labs(x = "CD22 Kd (M)", y = "Activity Gain") +
#   geom_smooth(method = "lm", se = FALSE, aes(group = 1))  +
#   geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$lwr), 
#             linetype = 2, color = "grey") +
#   geom_line(aes(x = neg_log10_Kd_M, y = cd22_kd_lims$upr), 
#             linetype = 2, color = "grey") +
#   scale_color_manual(values = c(
#     "interface" = "dodgerblue3",
#     "non-interface" = "firebrick",
#     "parental" = "black"
#   ),name = "evolution") +
#   scale_x_continuous(labels = function(x) parse(text = paste0("10^-", x))) +
#   theme(legend.position = "none")
# cd22_kd_plot
# cd22_mpnn_kd_cor <- cor.test(mpnn_cd22$neg_log10_Kd_M, mpnn_cd22$overall_activity_gain, method = "spearman")
# cd22_mpnn_kd_cor
# 
# bcma_kd_plot | cd22_kd_plot
# 
# bcma_cd22_kd_plot = cowplot::plot_grid(bcma_kd_plot,cd22_kd_plot,ncol=2)
# bcma_cd22_kd_plot
# cowplot::ggsave2("../plots/bcma_cd22_kd_activity_gain.pdf",bcma_cd22_kd_plot,height=1.6,width=3.6)

