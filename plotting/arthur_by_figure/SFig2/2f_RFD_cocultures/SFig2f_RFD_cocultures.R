library(BuenColors)
library(dplyr)
library(readxl)
library(reshape2)

raw_in <- read_csv("../data/DNCT_rev/RFD5_CAR-J_CD69_cocultures.csv") %>% data.frame()

melt_df <- raw_in[, c(
  "binder_name",
  "coculture.CAR.Only",
  "coculture.K562",
  "coculture.MM1S",
  "coculture.Raji"
)] %>%
  reshape2::melt(id.vars = c("binder_name")) %>%
  mutate(variable = gsub("coculture.", "", variable))

order_cl <- c("CAR.Only", "K562","Raji", "MM1S")
order_binder_controls <- c("No binder", "i16_51", "i1_29", "i68_76", "i9_54", "i55_32","Abecma")

p1 <- melt_df %>%
  mutate(variable = factor(variable, levels = rev(order_cl))) %>%
  mutate(binder_name = factor(binder_name, levels = order_binder_controls)) %>%
  ggplot(aes(x = binder_name, y = variable, fill = value)) +
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") + labs(x = "BCMA MPNN Draw", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

cowplot::ggsave2(p1, file = "../plots/BCMA_RFD_co-culture.pdf", width = 2.7, height = 2.6)
