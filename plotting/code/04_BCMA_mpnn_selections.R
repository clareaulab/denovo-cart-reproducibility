library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

dt <- fread("../data/BCMA_MPNN_filter_passed_unselected_included.csv") %>% data.frame() %>%
  mutate(cl_redesign = ifelse(mpnn_redesign_selection %in% c("non_interface_selected", "interface_selected"), "selected","not_selected")) %>%
  mutate(design = ifelse(mpnn_redesign_selection %in% c("interface_selected", "non_interface_unselected"), "interface", "non_interface"))
parental <- data.frame(what = colnames(dt), t(dt[dt$mpnn_redesign_selection == "parental",])); parental <- parental[nchar(parental[[2]]) < 150,]
rownames(parental) <- 1:length(rownames(parental))
dt <- dt[dt$mpnn_redesign_selection != "parental",]
dt <- dt %>% mutate(size = (cl_redesign == "selected")*1 + 1)

# Make diversity plots
parental %>% filter(what %in% c("Average_n_InterfaceHbonds", "Average_dSASA"))
p0 <- ggplot(dt %>% arrange(cl_redesign), aes(x = Average_dSASA, y = Average_n_InterfaceHbonds, color = mpnn_redesign_selection)) + 
  geom_point() +
  geom_point(data = data.frame(Average_dSASA = 1473.928, Average_n_InterfaceHbonds = 8), color = "black") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" )) + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")

parental %>% filter(what %in% c("binder_seq_net_charge", "n_hydrophobic_residues"))

p1 <- ggplot(shuf(dt) %>% arrange(cl_redesign), aes(x = binder_seq_net_charge, y = n_hydrophobic_residues, )) + 
  geom_point(aes(color = mpnn_redesign_selection)) +
  geom_point(data = data.frame(binder_seq_net_charge = 0.7949567, n_hydrophobic_residues = 18), color = "black") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" )) + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(p0, p1, nrow = 1), 
                 file = "../plots/bcma_mpnn_diversity.pdf", width = 3.3, height = 1.6)

# make plots verifiying quality
parental %>% filter(what %in% c("Average_protein_iptm", "Average_ipSAE","Average_complex_plddt", "Average_pae_interaction", ""))

p2 <- ggplot(shuf(dt), aes(x = cl_redesign, y = Average_protein_iptm,color = mpnn_redesign_selection)) +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  geom_hline(yintercept =0.9224972, linetype = 2) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" ))

p3 <- ggplot(shuf(dt), aes(x = cl_redesign, y = Average_ipSAE, color = mpnn_redesign_selection)) +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  geom_hline(yintercept =0.8822802, linetype = 2) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" ))

p4 <- ggplot(shuf(dt), aes(x = cl_redesign, y = Average_complex_plddt, color = mpnn_redesign_selection)) +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  geom_hline(yintercept =0.8126818, linetype = 2) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" ))

p5 <- ggplot(shuf(dt), aes(x = cl_redesign, y = Average_pae_interaction*31, color = mpnn_redesign_selection)) +
  geom_quasirandom(size = 0.5) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  geom_hline(yintercept = 0.08541517*31, linetype = 2) + 
  L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue3", "lightblue", "firebrick", "lightpink" ))

cowplot::ggsave2(cowplot::plot_grid(p2, p3, p4, p5, nrow = 1), 
                 file = "../plots/bcma_mpnn_metrics_verify.pdf", width = 6.6, height = 1.4)



