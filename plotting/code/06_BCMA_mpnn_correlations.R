library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

binder_attributes <- read_xlsx("../data/bcma_mpnn_binders_fullattributes.xlsx") %>% data.frame()

bcma_df <- binder_attributes %>% filter(binder_experiment=="BCMA", binder_name != "BCMA_Abecma")
bcma_df <- bcma_df %>% mutate(
  binder_seq_net_charge = as.numeric(binder_seq_net_charge),
  n_hydrophobic_residues = as.numeric(n_hydrophobic_residues),
  n_InterfaceHbonds = as.numeric(n_InterfaceHbonds),
  Average_dSASA = as.numeric(Average_dSASA)
)

p1_pred <- predict(lm(coculture.CAR.Only ~ binder_seq_net_charge, bcma_df), 
                se.fit = TRUE, interval = "confidence")
p1_limits <- as.data.frame(p1_pred$fit)

## Show associaiton between net charge and tonic activity
p1 <- df %>%
  ggplot(aes(x = binder_seq_net_charge,y=coculture.CAR.Only)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Binder Net Charge", y = "%CD69+ CAR alone") +
  geom_smooth(method = "lm", se = FALSE)  +
  geom_line(aes(x = binder_seq_net_charge, y = p1_limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = binder_seq_net_charge, y = p1_limits$upr), 
            linetype = 2,  color = "grey") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" )) +
  theme(legend.position = "none")

p1
cowplot::ggsave2(p1, file = "../plots/fig3d_bcma_tonic.pdf", width = 1.8, height = 1.8)
cor.test(df$binder_seq_net_charge, df$coculture.CAR.Only, method = "spearman")

## Showing no association between hydrophobic residues and tonic signaling
bcma_df_filtered = bcma_df %>% filter(!is.na(n_hydrophobic_residues))
p2_pred <- predict(lm(coculture.CAR.Only ~ n_hydrophobic_residues, bcma_df_filtered), 
                se.fit = TRUE, interval = "confidence")
p2_limits <- as.data.frame(p2_pred$fit)
p2 <- bcma_df_filtered %>%
  ggplot(aes(x = n_hydrophobic_residues, y=coculture.CAR.Only)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "# Hydrophobic Residues", y = "%CD69+ CAR alone") +
  geom_smooth(method = "lm", se = FALSE, color = "gray")  +
  geom_line(aes(x = n_hydrophobic_residues, y = p2_limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = n_hydrophobic_residues, y = p2_limits$upr), 
            linetype = 2,  color = "grey") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" )) +
  theme(legend.position = "none")

p2
cowplot::ggsave2(p2, file = "../plots/bcma_hydrophobic_tonic_activity.pdf", width = 1.8, height = 1.8)
cor.test(bcma_df_filtered$n_hydrophobic_residues, bcma_df_filtered$coculture.CAR.Only, method = "spearman")

## Showing no association between H-bond and activity gain
p3_pred <- predict(lm(overall_activity_diff ~ n_InterfaceHbonds, bcma_df_filtered), 
                se.fit = TRUE, interval = "confidence")
p3_limits <- as.data.frame(p3_pred$fit)
p3 <- bcma_df_filtered %>%
  ggplot(aes(x = n_InterfaceHbonds, y=overall_activity_diff)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Average # Hbonds", y = "%CD69+ (activated - alone)") +
  geom_smooth(method = "lm", se = FALSE, color = "gray")  +
  geom_line(aes(x = n_InterfaceHbonds, y = p3_limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = n_InterfaceHbonds, y = p3_limits$upr), 
            linetype = 2,  color = "grey") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" )) +
  theme(legend.position = "none")

p3
cowplot::ggsave2(p3, file = "../plots/bcma_hbond_activity.pdf", width = 1.8, height = 1.8)
cor.test(bcma_df_filtered$n_InterfaceHbonds, bcma_df_filtered$overall_activity_diff, method = "spearman")

## Showing no association between dSASA and activity gain
p4_pred <- predict(lm(overall_activity_diff ~ Average_dSASA, bcma_df_filtered), 
                se.fit = TRUE, interval = "confidence")
p4_limits <- as.data.frame(p4_pred$fit)
p4 <- bcma_df_filtered %>%
  ggplot(aes(x = Average_dSASA, y=overall_activity_diff)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Average dSASA", y = "%CD69+ (activated - alone)") +
  geom_smooth(method = "lm", se = FALSE, color = "gray")  +
  geom_line(aes(x = Average_dSASA, y = p4_limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = Average_dSASA, y = p4_limits$upr), 
            linetype = 2,  color = "grey") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" )) +
  theme(legend.position = "none")

p4
cowplot::ggsave2(p4, file = "../plots/bcma_dSASA_activity.pdf", width = 1.8, height = 1.8)
cor.test(bcma_df_filtered$Average_dSASA, bcma_df_filtered$overall_activity_diff, method = "spearman")

bcma_correlations = cowplot::plot_grid(p1,p2,p3,p4,ncol=4)

cowplot::ggsave2(bcma_correlations, file = "../plots/bcma_mpnn_correlations.pdf", width = 1.7*4, height = 1.6)


