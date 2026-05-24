library(BuenColors)
library(dplyr)
library(cowplot)
library(ggplot2)

vhh_df = read.csv("../data/combined_boltz_df_cd19_vhh_bc_style_merged_subset.csv")

vhh_df <- vhh_df %>%
  mutate(
    nb_number = as.integer(str_extract(prediction_name, "(?<=Nb_)\\d+")),
    nb_category = case_when(
      nb_number >= 1  & nb_number <= 39 ~ "BoltzGen Ranked",
      nb_number >= 40 & nb_number <= 68 ~ "ipSAE Ranked",
      TRUE ~ NA_character_
    )
  )

vhh_df

write.csv(vhh_df,"../data/combined_boltz_df_cd19_vhh_bc_style_merged_subset_nb_assigned.csv",row.names = FALSE)

selection_ipSAE_plot = ggplot(vhh_df,aes(x=Average_ipSAE_max,fill=nb_category)) + 
  geom_histogram(color = "black") +
  pretty_plot() + L_border() +
  scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values = c("BoltzGen Ranked" = "firebrick", "ipSAE Ranked" = "dodgerblue3")) +
  theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())

selection_ipSAE_plot

ggsave("../plots/CD19_Boltzgen_VHH_selection_ipSAE_plot.pdf",selection_ipSAE_plot,dpi=300,height=1.4,width=2)

rfd3_df = read.csv("../data/2026_05_19_all_binders_with_labels.csv")
rfd3_df_subset = rfd3_df %>% filter(campaign_short == "RFD3 1 (CD19)")

rfd3_ipsae_plot = ggplot(rfd3_df_subset,aes(x=Average_ipSAE_max)) + 
  geom_histogram(color="black",fill = "dodgerblue3") +
  pretty_plot() + L_border() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
  theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
rfd3_ipsae_plot
ggsave("../plots/CD19_RFD3_minibinder_selection_ipSAE_plot.pdf",rfd3_ipsae_plot,dpi=300,height=1.4,width=2)

