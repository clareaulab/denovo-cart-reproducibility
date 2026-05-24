library(BuenColors)
library(dplyr)
library(readxl)
library(reshape2)

raw_in <- read_xlsx("../data/CD19_CD81-KD_coculture.xlsx") %>% data.frame()
colnames(raw_in)

melt_df <- raw_in[, c(
  "binder_name",
  "coculture.CAR.Only",
  "coculture.K562",
  "coculture.K562.CD19.OE",
  "coculture.K562.CD19.OE.CD81.KD"
)] %>%
  reshape2::melt(id.vars = c("binder_name")) %>%
  mutate(variable = gsub("coculture.", "", variable))

order_cl <- c("CAR.Only", "K562","K562.CD19.OE", "K562.CD19.OE.CD81.KD")
order_binder_controls <- c("UTD", "NB", "FMC63", "C1", "C2", "C2.N2")

p1 <- melt_df %>%
  mutate(variable = factor(variable, levels = rev(order_cl))) %>%
  mutate(binder_name = factor(binder_name,
                              levels = c("UTD", "NB", "FMC63", "C1", "C2", "C2.N2"))) %>%
  ggplot(aes(x = binder_name, y = variable, fill = value)) +
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") + labs(x = "BCMA MPNN Draw", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

cowplot::ggsave2(p1, file = "../plots/full_parental_CAR_CD19.pdf", width = 4, height = 1.5)

