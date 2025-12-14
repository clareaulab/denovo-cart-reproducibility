library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)
library(data.table)
library(purrr)

## Read data
bcma_car_df <- fread("../data/BCMA_BindCraft_Campaign_CAR_Activity.csv") %>%
  data.frame() %>%
  rowwise() %>% # Treat each row as a separate group
  mutate(
    mean_pos_activity = mean(c_across(K562.BCMA:RPMI.8226), na.rm = TRUE),
    activity_gain = mean_pos_activity - CAR.alone
  ) %>%
  ungroup()

## Extract metrics columns
average_cols <- bcma_car_df %>%
  dplyr::select(starts_with("Average")) %>%
  names()
other_insilico_cols <- c("Length","MPNN_score","MPNN_seq_recovery")#,"hydrophobic_count","net_charge")
cols_to_cor <- c(average_cols,other_insilico_cols)

bcma_car_df <- bcma_car_df %>%
  mutate(across(all_of(cols_to_cor), as.numeric))

## Helper function to calculate correlations
calculate_spearman <- function(data, x_col) {
  # 1. Create a temporary subset with only the two columns needed and remove NAs
  temp_data <- data %>% 
    dplyr::select(all_of(x_col), activity_gain) %>%
    na.omit()
  n_obs <- nrow(temp_data)
  # 2. Check for sufficient finite observations (cor.test requires N >= 4)
  if (n_obs < 4) {
    # Return NA values and a status note if sample size is too small
    return(data.frame(
      Feature = x_col,
      Spearman_Rho = NA,
      P_Value = NA,
      N_Observations = n_obs,
      Status = "Insufficient Data (N < 4)"
    ))
  }
  # 3. Run the Spearman test on the complete observations
  test_result <- cor.test(
    x = temp_data[[x_col]], 
    y = temp_data$activity_gain, 
    method = "spearman"
  )
  # 4. Extract and return the key statistics
  return(data.frame(
    Feature = x_col,
    Spearman_Rho = test_result$estimate[["rho"]],
    P_Value = test_result$p.value,
    N_Observations = n_obs,
    Status = "OK"
  ))
}

correlation_results <- map_dfr(cols_to_cor, ~calculate_spearman(bcma_car_df, .x))

dl_metrics = correlation_results$Feature[grepl("pLDDT|pAE|pTM|RMSD|MPNN|Length",correlation_results$Feature)]
rosetta_metrics = correlation_results$Feature[grepl("Clashes|Energy|Surface|Interface|SASA|dG|Loop|Helix|BetaSheet|Shape|PackStat",correlation_results$Feature)]


correlation_results = correlation_results %>% mutate(
  metric_type = case_when(
    Feature %in% dl_metrics ~ "AlphaFold",
    Feature %in% rosetta_metrics ~ "Rosetta",
    TRUE ~ "Other"
  )
) %>% mutate(
    rho = round(Spearman_Rho, 4),
    pval = format.pval(P_Value, digits = 4)
)

correlation_results = correlation_results %>% mutate(padj = p.adjust(P_Value, method = "BH"))
filtered_results <- correlation_results %>%
  dplyr::filter(!is.na(Feature) & !is.na(Spearman_Rho)) %>%
  mutate(rho_rank = rank(-Spearman_Rho))

bcma_bc_feature_plot = ggplot(filtered_results, # Use the filtered data frame
                              aes(x = rho_rank, 
                                  y = Spearman_Rho, 
                                  fill=P_Value<0.05)) + 
  #geom_bar(stat=aes(Feature)) +
  geom_bar(stat = "identity",color = "black") + 
  labs(x = "Rank Ordered Features", y = "Activity Gain Correlation") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "dodgerblue3"),name = "Direction") +
  theme(
    legend.position = "none",
    #axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y=element_blank(),
    axis.line = element_line(colour = 'black', size = 0.5),
    ) +
  scale_x_continuous(breaks = c(1,35),expand = c(0, 0)) #+
  #labs()
  #scale_x_discrete(expand = c(0, 0))


bcma_bc_feature_plot
cowplot::ggsave2(bcma_bc_feature_plot, file = "../plots/bcma_bc_feature_plot.pdf", width = 1.8, height = 1.8)
cowplot::ggsave2(bcma_bc_feature_plot, file = "../plots/bcma_bc_feature_plot.pdf", width = 1.6, height = 1.3, units = "in")

bcma_bc_feature_dot_plot = ggplot(filtered_results, # Use the filtered data frame
                              aes(x = rho_rank, 
                                  y = Spearman_Rho, 
                                  color=P_Value<0.05)) + 
  #geom_bar(stat=aes(Feature)) +
  geom_point() +
  labs(x = "Rank Ordered Features", y = "Activity Gain Correlation") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "dodgerblue3"),name = "Direction") +
  theme(
    legend.position = "none",
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.line = element_line(colour = 'black', size = 0.5),
  ) +
  scale_x_continuous(breaks = c(1,35)) 
bcma_bc_feature_dot_plot

bcma_bc_feature_dot_plot
cowplot::ggsave2(bcma_bc_feature_dot_plot, file = "../plots/bcma_bc_feature_dot_plot.pdf", width = 1.8, height = 1.8)

## Now focus on the dSASA association
bcma_car_df_filtered = bcma_car_df %>% filter(grepl("BCMA",binder_name))
pred <- predict(lm(activity_gain ~ Average_dSASA, data = bcma_car_df_filtered),
                se.fit = TRUE, interval = "confidence")
plot_data <- cbind(bcma_car_df_filtered, pred$fit)
plot_data <- plot_data[order(plot_data$Average_dSASA), ]
plot_data$to_highlight = plot_data$binder_name_short %in% c("B4","B5")

dSASA_activity_gain_plot <- plot_data %>%
  ggplot(aes(x = Average_dSASA, y = activity_gain, color=to_highlight)) +
  # Use geom_smooth for the main regression line (simpler than geom_line)
  geom_smooth(aes(group = 1), method = "lm", se = FALSE) + 
  # Use geom_line for the confidence intervals, mapping to the new columns
  geom_line(aes(y = lwr), linetype = 2, color = "grey") +
  geom_line(aes(y = upr), linetype = 2, color = "grey") +
  geom_point() +
  # ADDED: segment.color to draw the connector line
  #geom_text_repel(size = 2, color = "firebrick") + 
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_manual(values=c("black","firebrick")) +
  labs(x = "Average dSASA", y = "Activity Gain") +
  theme(legend.position = "none")
  

dSASA_activity_gain_plot
cowplot::ggsave2(dSASA_activity_gain_plot, file = "../plots/bcma_dSASA_gain_cor.pdf", width = 2.2, height = 1.4)
#cor.test(bcma_df$neg_log10_Kd_M, bcma_df$overall_activity_diff, method = "spearman")

# # Constants for Charge Calculation
# pKa_values <- list(
#   D=3.9, E=4.1, C=8.3, Y=10.1, H=6.0, K=10.5, R=12.5, N_Term=8.0, C_Term=3.1
# )
# # Hydrophobic Residues (A, I, L, M, F, W, V, P, G)
# HYDROPHOBIC_RESIDUES <- "A|I|L|M|F|W|V|P|G" 
# 
# count_hydrophobic_residues <- function(sequence) {
#   sum(str_count(sequence, HYDROPHOBIC_RESIDUES))
# }
# 
# net_charge <- function(sequence, pH = 7.4) {
#   counts <- as.data.frame(table(str_split(sequence, "")))
#   colnames(counts) <- c("AA", "Count")
#   
#   pos_aas <- c("H", "K", "R")
#   neg_aas <- c("D", "E", "C", "Y")
#   
#   pos_charge <- sum(sapply(pos_aas, function(aa) {
#     count <- counts$Count[counts$AA == aa]
#     ifelse(length(count) == 0, 0, count * (1 / (1 + 10^(pH - pKa_values[[aa]]))))
#   }))
#   
#   neg_charge <- sum(sapply(neg_aas, function(aa) {
#     count <- counts$Count[counts$AA == aa]
#     ifelse(length(count) == 0, 0, count * (-1 / (1 + 10^(pKa_values[[aa]] - pH))))
#   }))
#   
#   pos_charge <- pos_charge + (1 / (1 + 10^(pH - pKa_values$N_Term)))
#   neg_charge <- neg_charge + (-1 / (1 + 10^(pKa_values$C_Term - pH)))
#   
#   pos_charge + neg_charge
# }
# 
# calculate_binder_properties <- function(sequence, pH = 7.4) {
#   if (!is.character(sequence) || length(sequence) == 0 || sequence == "" || is.na(sequence)) {
#     return(list(hydrophobic_count = 0, net_charge = 0))
#   }
#   
#   list(
#     hydrophobic_count = count_hydrophobic_residues(sequence),
#     net_charge = net_charge(sequence, pH)
#   )
# }
# 
# bcma_car_df <- bcma_car_df %>%
#   # Apply the analysis function to the 'seq' column row-wise
#   mutate(
#     properties = map(seq, calculate_binder_properties, pH = 7.4),
#     
#     # Extract results into the final target columns
#     binder_n_hydrophobic_residues = map_dbl(properties, "hydrophobic_count"),
#     binder_seq_net_charge = map_dbl(properties, "net_charge")
#   ) %>%
#   # Remove the temporary list column
#   dplyr::select(-properties)

