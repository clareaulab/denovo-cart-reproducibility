library(BuenColors)
library(dplyr)
library(data.table)

# Load summarized NetMHC outputs
bcma_netmhc <- read.csv("../data/BCMA_l59_mpnn_NetMHC.csv") %>% mutate(antigen="BCMA")# %>% select(-is_top_freq_allele)
cd19_l112_netmhc <- read.csv("../data/CD19_l112_mpnn_NetMHC.csv") %>% mutate(antigen="CD19")
cd19_l115_netmhc <- read.csv("../data/CD19_l115_mpnn_NetMHC.csv") %>% mutate(antigen="CD19")
cd22_netmhc <- read.csv("../data/CD22_l61_2_mpnn_NetMHC.csv") %>% mutate(antigen="CD22")#%>% select(-is_top_freq_allele)
binder_netmhc <- bind_rows(bcma_netmhc,cd19_l112_netmhc,cd19_l115_netmhc,cd22_netmhc) %>%
  distinct(Original_ID, Allele, .keep_all = TRUE)

## Get scFV lengths
abecma_length =  nchar("DIVLTQSPPSLAMSLGKRATISCRASESVTILGSHLIHWYQQKPGQPPTLLIQLASNVQTGVPARFSGSGSRTDFTLTIDPVEEDDVAVYYCLQSRTIPRTFGGGTKLEIKGSTSGSGKPGSGEGSTKGQIQLVQSGPELKKPGETVKISCKASGYTFTDYSINWVKRAPGKGLKWMGWINTETREPAYAYDFRGRFAFSLETSASTAYLQINNLKYEDTATYFCALDYSYAMDYWGQGTSVTVSSAAA")
fmc63_length = nchar("DIQMTQTTSSLSASLGDRVTISCRASQDISKYLNWYQQKPDGTVKLLIYHTSRLHSGVPSRFSGSGSGTDYSLTISNLEQEDIATYFCQQGNTLPYTFGGGTKLEITGSTSGSGKPGSGEGSTKGEVKLQESGPGLVAPSQSLSVTCTVSGVSLPDYGVSWIRQPPRKGLEWLGVIWGSETTYYNSALKSRLTIIKDNSKSQVFLKMNSLQTDDTAIYYCAKHYYYGGSYAMDYWGQGTSVTVSS")
m971_length = nchar("QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCAREVTGDLEDAFDIWGQGTMVTVSSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQTIWSYLNWYQQRPGKAPNLLIYAASSLQSGVPSRFSGRGSGTDFTLTISSLQAEDFATYYCQQSYSIPQTFGQGTKLEIKAAA")

# Add more detailed annotations
binder_netmhc = binder_netmhc %>% mutate(
  binder_type_detailed = case_when(
    Original_ID == "BCMA_Abecma" ~ "scFV (ABECMA)",
    Original_ID == "CD19_FMC63" ~ "scFV (FMC63)",
    Original_ID == "CD22_m971" ~ "scFV (m971)",
    TRUE ~ "de novo"
  ),
  binder_length = case_when(
    Original_ID == "BCMA_Abecma" ~ abecma_length,
    Original_ID == "CD19_FMC63" ~ fmc63_length,
    Original_ID == "CD22_m971" ~ m971_length,
    startsWith(Original_ID, "BCMA_l59") ~ 59,
    startsWith(Original_ID, "CD19_l112") ~ 112,
    startsWith(Original_ID, "CD19_l115") ~ 115,
    startsWith(Original_ID, "CD22_l61") ~ 61,
    TRUE ~ NA
  ),
  total_9mers = binder_length - 8,
  strong_binder_normalized = Strong.Binder / total_9mers,
)

# Sum across all alleles
binder_netmhc_summed <- binder_netmhc %>% 
  group_by(Original_ID) %>%
  summarise(
    Total_Strong_Binders = sum(Strong.Binder, na.rm = TRUE),
    mean_strong_binder_fraction = mean(strong_binder_normalized, na.rm=TRUE),
    max_strong_binder_fraction = max(strong_binder_normalized, na.rm=TRUE),
    binder_type = first(binder_type),
    binder_type_detailed = first(binder_type_detailed),
    antigen = first(antigen)
  ) %>% ungroup()

binder_netmhc_summed

## Test to see wald test coefficient on binder type
#burden_lm_test = lm(Total_Strong_Binders ~ antigen + binder_type, binder_netmhc_summed)
burden_lm_test = lm(mean_strong_binder_fraction ~ antigen + binder_type, binder_netmhc_summed)

lm_table = summary(burden_lm_test)$coefficients
lm_table

binder_netmhc_summed
p0 <- ggplot(shuf(binder_netmhc_summed), aes(x = antigen, y=mean_strong_binder_fraction, color=binder_type)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  geom_quasirandom(size = 0.5) +
  pretty_plot(fontsize = 8) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(y="mean MHC-I Burden",x="Antigen")
  #scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink", "forestgreen","palegreen")) +
  #labs(y="NetMHC Class I Burden",x="Antigen")
p0
#cowplot::ggsave2(p0, file = "../plots/de_novo_NetMHC.pdf", width = 1.65, height = 1.5)
cowplot::ggsave2(p0, file = "../plots/de_novo_scFV_NetMHC_mean_9mer_fraction.pdf", width = 1.65, height = 1.5)


# 
# 
# binder_netmhc_normalized = binder_netmhc_summed %>% mutate(
#   Total_Strong_Binders_Normalized = case_when(
#     Original_ID == "BCMA_Abecma" ~ Total_Strong_Binders / (abecma_length-8),
#     Original_ID == "CD19_FMC63" ~ Total_Strong_Binders / fmc63_length,
#     Original_ID == "CD22_m971" ~ Total_Strong_Binders / m971_length,
#     startsWith(Original_ID, "BCMA_l59") ~ Total_Strong_Binders / 59,
#     startsWith(Original_ID, "CD19_l112") ~ Total_Strong_Binders / 112,
#     startsWith(Original_ID, "CD19_l115") ~ Total_Strong_Binders / 115,
#     startsWith(Original_ID, "CD22_l61") ~ Total_Strong_Binders / 61,
#     TRUE ~ 10000
#   )
# )
# 
# binder_netmhc_normalized
# 
# p1 <- ggplot(shuf(binder_netmhc_normalized), aes(x = antigen, y=Total_Strong_Binders_Normalized, color=binder_type)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
#   geom_quasirandom(size = 0.5) +
#   pretty_plot(fontsize = 8) + 
#   L_border() + theme(legend.position = "none") +
#   scale_color_manual(values = c("dodgerblue3", "firebrick")) +
#   #scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink", "forestgreen","palegreen")) +
#   labs(y="Normalized MHC-I Burden",x="Antigen")
# p1

# data.frame(
#   campaign = c("RF1", "RF2", "RF3", "RF4", "RF5", "zBC1"),
#   Q2 = c(0.4, 0.5, 0.3, 0.4, 1.9, 4.9),
#   Q3 = c(52.2, 61.3, 65.5, 63.4, 77.2, 40.7)
# ) %>% mutate(rat = Q2/(Q2+Q3)) %>%
#   ggplot(aes(x = campaign, y = rat*100)) + 
#   geom_bar(stat = "identity", fill = "lightgrey", color = "black") + labs(x = "", y = "") +
#   pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10)) -> p1 
# cowplot::ggsave2(p1, file = "../plots/fig1c_yeast.pdf", width = 2.1, height = 1.8)
