library(BuenColors)
library(dplyr)
library(data.table)

# Load summarized NetMHC outputs
bcma_netmhc <- read.csv("../data/BCMA_l59_mpnn_NetMHC.csv") %>% mutate(antigen="BCMA")
cd19_netmhc <- read.csv("../data/CD19_l186_mpnn_NetMHC.csv") %>% mutate(antigen="CD19")
cd22_netmhc <- read.csv("../data/CD22_l61_2_mpnn_NetMHC.csv") %>% mutate(antigen="CD22")
binder_netmhc <- rbind(bcma_netmhc,cd19_netmhc,cd22_netmhc)

# Add more detailed annotations
binder_netmhc = binder_netmhc %>% mutate(
  binder_type_detailed = case_when(
    Original_ID == "BCMA_Abecma" ~ "scFV (ABECMA)",
    Original_ID == "CD19_FMC63" ~ "scFV (FMC63)",
    Original_ID == "CD22_m971" ~ "scFV (m971)",
    TRUE ~ "de novo"
  )
)

# Sum across all alleles
binder_netmhc_summed <- binder_netmhc %>% 
  group_by(Original_ID) %>%
  summarise(
    Total_Strong_Binders = sum(Strong.Binder, na.rm = TRUE),
    binder_type = first(binder_type),
    binder_type_detailed = first(binder_type_detailed),
    antigen = first(antigen)
  ) %>% ungroup()

binder_netmhc_summed
p0 <- ggplot(shuf(binder_netmhc_summed), aes(x = antigen, y=Total_Strong_Binders, color=binder_type)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  geom_quasirandom(size = 0.5) +
  pretty_plot(fontsize = 8) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  #scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink", "forestgreen","palegreen")) +
  labs(y="NetMHC Class I Burden",x="Antigen")
cowplot::ggsave2(p0, file = "../plots/de_novo_NetMHC.pdf", width = 1.65, height = 1.8)


# data.frame(
#   campaign = c("RF1", "RF2", "RF3", "RF4", "RF5", "zBC1"),
#   Q2 = c(0.4, 0.5, 0.3, 0.4, 1.9, 4.9),
#   Q3 = c(52.2, 61.3, 65.5, 63.4, 77.2, 40.7)
# ) %>% mutate(rat = Q2/(Q2+Q3)) %>%
#   ggplot(aes(x = campaign, y = rat*100)) + 
#   geom_bar(stat = "identity", fill = "lightgrey", color = "black") + labs(x = "", y = "") +
#   pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10)) -> p1 
# cowplot::ggsave2(p1, file = "../plots/fig1c_yeast.pdf", width = 2.1, height = 1.8)
