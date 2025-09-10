library(BuenColors)
library(dplyr)

data.frame(
  campaign = c("RF1", "RF2", "RF3", "RF4", "RF5", "zBC1"),
  Q2 = c(0.4, 0.5, 0.3, 0.4, 1.9, 4.9),
  Q3 = c(52.2, 61.3, 65.5, 63.4, 77.2, 40.7)
) %>% mutate(rat = Q2/(Q2+Q3)) %>%
  ggplot(aes(x = campaign, y = rat*100)) + 
  geom_bar(stat = "identity", fill = "lightgrey", color = "black") + labs(x = "", y = "") +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10)) -> p1 
cowplot::ggsave2(p1, file = "../plots/fig1c_yeast.pdf", width = 2.1, height = 1.8)
