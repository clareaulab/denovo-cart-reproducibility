library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

data.frame(
  binder =  c("aUTD", "bNB", "cFMC63", "dl112", "el115", "fl186"),
  Yenrich = c(NA, NA, NA, 3.77, 6.2, 9.74),
  Y100nM = c(NA, NA, NA, 80.75, 89.0, 18.35),
  Y1000nM = c(NA, NA, NA, 81.89, 91.45, 56.1)
) %>%
  reshape2::melt(id.vars = "binder") %>%
  mutate(variable = factor(as.character(variable), levels = rev(c("Yenrich", "Y100nM", "Y1000nM")))) %>%
  ggplot(.,aes(x = binder, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100))  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none")  -> p1
cowplot::ggsave2(p1, file = "../plots/CD19_parental_binding_yeast.pdf", width= 1.4, height = 0.8)



data.frame(
  binder = c("aUTD", "bNB", "cFMC63", "dl112", "el115", "fl186"),
  MYC = c(1.16, 90.3, 93.0, 79.2, 54.3, 54.6),
  CAR_100nM = c(0.5, 0.35, 98, 91.3, 53.4, 0),
  CAR_1000nM = c(0.38, 0.25, 97.3, 93.8, 67.1, 1.59)
) %>% 
  reshape2::melt(id.vars = "binder") %>%
  mutate(variable = factor(as.character(variable), levels = rev(c("MYC", "CAR_100nM", "CAR_1000nM")))) %>%
  ggplot(.,aes(x = binder, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("brewer_purple"), limits = c(0, 100))  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none")  -> p1
cowplot::ggsave2(p1, file = "../plots/CD19_parental_binding_CAR-purple.pdf", width= 1.4, height = 0.8)
