library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

df <- fread("../data/BCMA_MPNN_netcharge-tonic.csv") %>%
  data.frame() 

pred <- predict(lm(coculture.CAR.Only ~ binder_seq_net_charge, df), 
                se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)

p1 <- df %>%
  ggplot(aes(x = binder_seq_net_charge,y=coculture.CAR.Only)) + 
  geom_point(aes(color = mpnn_redesign_type)) + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Binder net charge", y = "%CD69+ CAR alone") +
  geom_smooth(method = "lm", se = FALSE)  +
  geom_line(aes(x = binder_seq_net_charge, y = limits$lwr), 
            linetype = 2, color = "grey") +
  geom_line(aes(x = binder_seq_net_charge, y = limits$upr), 
            linetype = 2,  color = "grey") +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "black" ))

p1
cowplot::ggsave2(p1 + theme(legend.position = "none"), file = "../plots/fig3d_bcma_tonic.pdf", width = 1.8, height = 1.8)
cor.test(df$binder_seq_net_charge, df$coculture.CAR.Only, method = "spearman")
