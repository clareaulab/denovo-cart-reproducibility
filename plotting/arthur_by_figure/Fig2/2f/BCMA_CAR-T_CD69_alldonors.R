library(BuenColors)
library(dplyr)
library(ggplot2)
library(reshape2)

raw_in <- read.csv("../data/BCMA_CD69_alldonors_master.csv", check.names = TRUE) %>% data.frame()

melt_df <- raw_in[c("binder_name", "coculture.CAR.Only", "coculture.K562",
                      "coculture.K562.BCMA.OE", "coculture.Raji",
                      "coculture.RPMI.8226", "Donor")] %>%
  reshape2::melt(id.vars = c("binder_name", "Donor")) %>%
  mutate(variable = gsub("coculture\\.", "", variable))

order_cl <- c("CAR.Only", "K562", "RPMI.8226", "Raji", "K562.BCMA.OE")

order_binders <- c("UTD", "NB", "B2", "B5", "B4", "Abecma")

p1 <- melt_df %>%
  filter(binder_name %in% order_binders) %>%
  mutate(variable = factor(as.character(variable), levels = rev(order_cl))) %>%
  mutate(binder_name = factor(as.character(binder_name), levels = order_binders)) %>%
  mutate(Donor = paste0("Donor ", Donor)) %>%
  ggplot(., aes(x = binder_name, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~Donor, nrow = 1) +
  theme(legend.position = "none") +
  labs(x = "BCMA binder", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p1

cowplot::ggsave2("../plots/CD69_heatmap_BCMA_withB4_alldonors.pdf",
                 p1, width = 5, height = 2.5)
