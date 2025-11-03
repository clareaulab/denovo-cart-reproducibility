library(BuenColors)
library(dplyr)
library(data.table)

dt = fread("../data/CD22_off_target_parental_vs_evolved.csv") %>% data.frame()

dt$to_highlight = dt$target_name %in% c("CD22_P20273","CXCR4_P61073")

p1 <- dt %>% 
  ggplot(aes(x = Average_ipSAE_parental,y=Average_ipSAE_evolved, color = to_highlight)) + 
  geom_point() + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "D1 Average ipSAE", y = "D1.N0 Average ipSAE") +
  scale_color_manual(values = c("black","firebrick")) +
  geom_abline(slope=1, intercept=0, color = "gray", linetype="dashed") +
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank()) + ylim(0,1) + xlim(0,1)
  
p1

cowplot::ggsave2(p1, file = "../plots/CD22_offtarget.pdf", width = 1.3, height = 1.3)

