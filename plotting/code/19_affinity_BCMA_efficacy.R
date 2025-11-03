library(BuenColors)
library(dplyr)
library(readxl)

binder_attributes <- read_xlsx("../data/bcma_mpnn_binders_fullattributes.xlsx") %>% data.frame()

bcma_df <- binder_attributes %>% filter(binder_experiment=="BCMA", binder_name != "BCMA_Abecma")
bcma_df <- bcma_df %>% mutate(
  n_hydrophobic_residues = as.numeric(n_hydrophobic_residues),
  n_InterfaceHbonds = as.numeric(n_InterfaceHbonds),
  Average_dSASA = as.numeric(Average_dSASA)
)

kd_to_nM_label <- function(x) {
  Kd_M <- 10^(-x)
  Kd_nM <- Kd_M * 1e9
  labels <- ifelse(
    abs(x) > 10,
    paste0(scientific(Kd_nM), " nM"),
    paste0(round(Kd_nM,0), " nM")
  )
  return(labels)
}

## Remove low quality kd before inferrence
bcma_high_conf_kd <- bcma_df %>% filter(neg_log10_Kd_M > 5)
pred <- predict(lm(overall_activity_diff ~ neg_log10_Kd_M, bcma_high_conf_kd), 
                se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)

p1 <- bcma_high_conf_kd %>% 
  ggplot(aes(x = neg_log10_Kd_M,y=overall_activity_diff)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Estimated Kd (nM)", y = "%CD69+ (activated - alone)") +
  geom_smooth(method = "lm", se = FALSE)  +
  geom_line(aes(x = neg_log10_Kd_M, y = limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = neg_log10_Kd_M, y = limits$upr), 
            linetype = 2,  color = "grey") +
  scale_x_continuous(labels=kd_to_nM_label) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" ))

p1

cowplot::ggsave2(p1 + theme(legend.position = "none"), file = "../plots/bcma_affinity_activity.pdf", width = 1.8, height = 1.8)
cor.test(bcma_df$neg_log10_Kd_M, bcma_df$overall_activity_diff, method = "spearman")







library(dplyr)
library(corrr)

bcma_df_filtered <- bcma_df %>%
  filter(!is.na(mpnn_redesign_type)) %>%
  filter(binder_name != "BCMA_Abecma") %>%
  filter(mpnn_redesign_type == "non-interface")

# Convert all columns that are NOT already numeric (e.g., 'character' or 'factor')
# to numeric (float). This is the key step.
numeric_df <- bcma_df_filtered %>%
  mutate(
    across(
      .cols = !where(is.numeric), # Select all columns that are NOT numeric
      .fns = as.numeric            # Apply as.numeric to those columns
    )
  )

target_column_name <- "overall_activity_diff"

# 2. Calculate the correlation vector
# cor(x, y) calculates the correlation of every column in x against every column in y.
# If y is a single column, it returns a vector of correlations.
correlation_vector <- cor(
  x = numeric_df,
  y = numeric_df[[target_column_name]],
  # 'pairwise.complete.obs' ensures correlations are calculated for each
  # pair using all available non-NA observations for that pair.
  use = "pairwise.complete.obs"
)

# 3. Convert the results to a clean, readable dataframe (using the tidyverse style)
library(dplyr)
library(tibble)

correlation_results <- as.data.frame(correlation_vector) %>%
  rownames_to_column(var = "Feature")

# 2. Get the actual name of the numeric column in the new data frame
# It's always the second column after 'rownames_to_column'
old_col_name <- names(correlation_results)[2]

# 3. Rename the column using the variable name (old_col_name) and the new name ("Correlation")
correlation_results <- correlation_results %>%
  rename(Correlation = all_of(old_col_name)) %>%  # Use all_of() for safe column selection/renaming with variables
  # Remove the row where the target column is correlated against itself (always 1)
  filter(Feature != target_column_name) %>%
  # Sort the results by the magnitude of the correlation (absolute value)
  arrange(desc(abs(Correlation)))

# 4. Print the final results
print(correlation_results)


correlation_results %>% arrange(-Correlation) %>% head(100)

p2 <- bcma_df %>% filter(mpnn_redesign_type == "non-interface") %>%
  ggplot(aes(x = as.numeric(Average_Interface_Hydrophobicity),y=overall_activity_diff)) +
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border()

p2

cd22_df = binder_attributes %>% filter(binder_experiment=="CD22", neg_log10_Kd_M > 5)# binder_name != "CD22_m971")

cd22_df_filtered <- cd22_df %>%
  filter(!is.na(mpnn_redesign_type)) %>%
  filter(binder_name != "BCMA_Abecma") %>%
  filter(mpnn_redesign_type == "non-interface")

# Convert all columns that are NOT already numeric (e.g., 'character' or 'factor')
# to numeric (float). This is the key step.
numeric_df <- bcma_df_filtered %>%
  mutate(
    across(
      .cols = !where(is.numeric), # Select all columns that are NOT numeric
      .fns = as.numeric            # Apply as.numeric to those columns
    )
  )

target_column_name <- "overall_activity_diff"

# 2. Calculate the correlation vector
# cor(x, y) calculates the correlation of every column in x against every column in y.
# If y is a single column, it returns a vector of correlations.
correlation_vector <- cor(
  x = numeric_df,
  y = numeric_df[[target_column_name]],
  # 'pairwise.complete.obs' ensures correlations are calculated for each
  # pair using all available non-NA observations for that pair.
  use = "pairwise.complete.obs"
)

# 3. Convert the results to a clean, readable dataframe (using the tidyverse style)
library(dplyr)
library(tibble)

correlation_results <- as.data.frame(correlation_vector) %>%
  rownames_to_column(var = "Feature")

# 2. Get the actual name of the numeric column in the new data frame
# It's always the second column after 'rownames_to_column'
old_col_name <- names(correlation_results)[2]

# 3. Rename the column using the variable name (old_col_name) and the new name ("Correlation")
correlation_results <- correlation_results %>%
  rename(Correlation = all_of(old_col_name)) %>%  # Use all_of() for safe column selection/renaming with variables
  # Remove the row where the target column is correlated against itself (always 1)
  filter(Feature != target_column_name) %>%
  # Sort the results by the magnitude of the correlation (absolute value)
  arrange(desc(abs(Correlation)))

# 4. Print the final results
print(correlation_results)


correlation_results %>% arrange(-Correlation) %>% head(100)

p2 <- bcma_df %>% filter(mpnn_redesign_type == "non-interface") %>%
  ggplot(aes(x = as.numeric(Average_Interface_Hydrophobicity),y=overall_activity_diff)) +
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border()

p2