library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(data.table)

# --- Data Preparation ---
dt = fread("../data/2025_12_14_all_binders_with_labels.csv") %>% data.frame()
dt_bcma = dt %>% filter(campaign == "BCMA_l59_MPNN")

bcma_mpnn_parental_rmsd = fread("../data/bcma_mpnn_parental_rmsd.csv")
dt_bcma = merge(dt_bcma, bcma_mpnn_parental_rmsd, by="binder_name")

bcma_parental_seq_vec_str = dt_bcma %>% filter(binder_id == "BCMA_l59_int_wildtype") %>% pull(binder_sequence)
bcma_parental_activity_gain = dt_bcma %>% filter(binder_id == "BCMA_l59_int_wildtype") %>% pull(overall_activity_gain)
parental_seq_vec <- strsplit(bcma_parental_seq_vec_str, "")[[1]]
seq_len <- length(parental_seq_vec)

evolved_df <- dt_bcma %>% 
  filter(binder_id != "BCMA_l59_int_wildtype") %>%
  dplyr::select(binder_sequence, mean_BCMA_positive_activity, overall_activity_gain, Average_binder_parental_RMSD) %>%
  mutate(RowID = as.factor(row_number())) %>%
  drop_na()

# --- Color and Group Mapping ---
aa_to_group <- c(
  "G"="Small Non-polar", "A"="Small Non-polar", "S"="Small Non-polar", "T"="Small Non-polar",
  "V"="Hydrophobic", "L"="Hydrophobic", "I"="Hydrophobic", "M"="Hydrophobic", "P"="Hydrophobic", "F"="Hydrophobic", "W"="Hydrophobic",
  "N"="Polar", "Q"="Polar", "Y"="Polar", "C"="Polar",
  "D"="Negative", "E"="Negative",
  "K"="Positive", "R"="Positive", "H"="Positive",
  "match"="Parental Match"
)

group_colors <- c(
  "Small Non-polar" = "#E69F00", "Hydrophobic" = "#009E73", 
  "Polar" = "#56B4E9", "Negative" = "#D55E00", "Positive" = "#CC79A7",
  "Parental Match" = "white"
)

heatmap_plot_data <- evolved_df %>%
  separate_rows(binder_sequence, sep = "") %>%
  filter(binder_sequence != "") %>%
  group_by(RowID) %>%
  mutate(Position = row_number()) %>%
  ungroup() %>%
  mutate(
    Parental_AA = parental_seq_vec[Position],
    is_mutated = ifelse(binder_sequence != Parental_AA, "1", "0"),
    plot_fill = ifelse(is_mutated == "1", binder_sequence, "match"),
    Chemical_Group = factor(aa_to_group[plot_fill])
  )

# --- Plots ---

p_left <- ggplot(evolved_df, aes(y = reorder(RowID, -as.numeric(RowID)), x = Average_binder_parental_RMSD)) +
  geom_col(fill = "dodgerblue3", width = 0.7) +
  scale_x_reverse(expand = expansion(mult = c(0.1, 0))) + 
  theme_minimal() + labs(x = "Parental RMSD", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank()) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, color = "black", linewidth = 0.8) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, color = "black", linewidth = 0.8)

p_mid <- ggplot(heatmap_plot_data, aes(x = Position, y = reorder(RowID, -as.numeric(RowID)), fill = Chemical_Group)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is_mutated == "0", "-", "")), color = "gray60", size = 3) +
  geom_text(aes(label = ifelse(is_mutated == "1", binder_sequence, "")), size = 2.5, fontface = "bold") +
  scale_x_continuous(breaks = 1:seq_len, labels = parental_seq_vec, expand = c(0,0), position = "top") +
  # name = NULL removes the title at the scale level
  scale_fill_manual(values = group_colors, 
                    breaks = c("Small Non-polar", "Hydrophobic", "Polar", "Negative", "Positive"),
                    name = NULL) + 
  theme_minimal() + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_blank(), panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) # Force NULL title in guide too

p_right <- ggplot(evolved_df, aes(y = reorder(RowID, -as.numeric(RowID)), x = overall_activity_gain)) +
  geom_col(fill = "firebrick", width = 0.7) +
  geom_vline(xintercept = bcma_parental_activity_gain, linetype="dashed", color = "gray40") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() + labs(x = "Activity Gain", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank()) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, color = "black", linewidth = 0.6) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, color = "black", linewidth = 0.6)

# Assembly WITH legend
final_plot_with_legend <- (p_left + p_mid + p_right) + 
  plot_layout(widths = c(1, 4, 1), guides = "collect") & 
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    # Negative top margin pulls the legend UP toward the plot
    # Positive bottom margin ensures it's not cut off
    legend.margin = margin(t = -15, b = 5), 
    legend.box.margin = margin(t = -10, b = 5),
    plot.margin = margin(t = 5, r = 2, b = 2, l = 2),
    legend.key.size = unit(0.4, "cm")
  )


ggsave("../plots/BCMA_evolved_final_with_legend.pdf", plot = final_plot_with_legend, width = 10, height = 4)

# Assembly WITHOUT legend
final_plot_no_legend <- (p_left + p_mid + p_right) + 
  plot_layout(widths = c(1, 4, 1)) & 
  theme(legend.position = "none", plot.margin = margin(5, 2, 5, 2))

print(final_plot_no_legend)
ggsave("../plots/BCMA_evolved_final_with_no_legend.pdf", plot = final_plot_no_legend, width = 10, height = 4)

