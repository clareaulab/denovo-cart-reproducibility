library(BuenColors)
library(dplyr)
library(data.table)

map_df = fread("../data/2025_10_31_denovo_CAR_T_ID_conversions.csv") %>% as.data.frame()
bcma_df = fread("../data/BCMA_l59_mpnn_NetMHC.csv") %>% as.data.frame()

bcma_df_merged = bcma_df %>% merge(map_df,by.x="Original_ID",by.y="order_binders",all.x=TRUE)

bcma_df_merged_filtered = bcma_df_merged %>%
  reshape2::melt(id.vars = c("new_name_final", "Strong Binder", "Allele")) %>% drop_na()

## Define the allele level to plot
allele_levels = bcma_df_merged_filtered %>%
  distinct(Allele) %>%
  mutate(
    allele_class = str_extract(Allele, "^HLA-([ABC])"),
    allele_number = str_extract(Allele, "[0-9]+$")
  ) %>%
  # Define the explicit order for the class
  mutate(
    class_order = case_when(
      allele_class == "HLA-A" ~ 1,
      allele_class == "HLA-B" ~ 2,
      allele_class == "HLA-C" ~ 3,
      TRUE ~ 99 # For any unexpected values
    )
  ) %>%
  arrange(class_order, allele_number) %>%
  pull(Allele)


bcma_df_merged_filtered <- bcma_df_merged_filtered %>%
  mutate(
    Allele = factor(Allele, levels = allele_levels),
  ) %>%
  mutate(Allele_Formatted = sub("^HLA-([A-Z])", "\\1:", Allele)) %>%
  mutate(Allele_Formatted = sub("([0-9]{2})([0-9]{2})$", "\\1*\\2", Allele_Formatted))

bcma_df_merged_filtered$Allele_Formatted

bcma_hla_heatmap = bcma_df_merged_filtered %>%
  ggplot(aes(y = new_name_final, x = Allele_Formatted, fill = `Strong Binder`)) +
  geom_tile(aes(fill = `Strong Binder`),color = "black") +
  #facet_wrap(~what) +
  scale_fill_gradientn(
    colors = jdb_palette("solar_blues"),
    #limits = c(0, 100),
    guide = "none" # Use guide = "none" to mimic theme(legend.position = "none") without conflict
  ) +
  theme_bw(base_size = 8) +
  theme(
    #axis.text.x = element_text(),
    #axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 5),
    axis.text.y = element_text(size = 5),
  ) +
  #geom_text(aes(label=`Strong Binder`),size=1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL)

bcma_hla_heatmap

cowplot::ggsave2(bcma_hla_heatmap, file = "../plots/BCMA_HLA_heatmap.pdf", width = 6.8, height = 2)


