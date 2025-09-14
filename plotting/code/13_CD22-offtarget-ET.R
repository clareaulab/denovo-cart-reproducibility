library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

valsdf <- fread("../data/CD22_ET_many.tsv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("ET_Ratio", "Cellline")) %>%
  filter(ET_Ratio != "i1:8")

binders_show <- c("No_Binder", "m971", "l61_819258_mpnn2")

p0 <- valsdf %>%
  filter(Cellline == "K562") %>%
  filter(variable %in% binders_show) %>%
  ggplot(aes(x = ET_Ratio, y = value, color = variable, group = variable)) + 
  geom_point() + geom_line() +
  scale_y_continuous(limits = c(0,100))+ pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + scale_color_manual(values = jdb_palette("corona"))

p1 <- valsdf %>%
  filter(Cellline == "K562CD22") %>%
  filter(variable %in% binders_show) %>%
  ggplot(aes(x = ET_Ratio, y = value, color = variable, group = variable)) + 
  geom_point() + geom_line() +
  scale_y_continuous(limits = c(0,100))+ pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + scale_color_manual(values = jdb_palette("corona"))

p2 <- valsdf %>%
  filter(Cellline == "Raji") %>%
  filter(variable %in% binders_show) %>%
  ggplot(aes(x = ET_Ratio, y = value, color = variable, group = variable)) + 
  geom_point() + geom_line() +
  scale_y_continuous(limits = c(0,100)) + pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + scale_color_manual(values = jdb_palette("corona"))

p3 <- valsdf %>%
  filter(Cellline == "RPMI8266") %>%
  filter(variable %in% binders_show) %>%
  ggplot(aes(x = ET_Ratio, y = value, color = variable, group = variable)) + 
  geom_point() + geom_line() +
  scale_y_continuous(limits = c(0,100)) + pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + scale_color_manual(values = jdb_palette("corona"))

cowplot::ggsave2(cowplot::plot_grid(p0, p1, p2, p3, nrow = 1), file = "../plots/CD22ET_grid_4_coculture.pdf", width = 7.6, height = 1.8)

valsdf %>%
  filter(Cellline == "RPMI8266") %>%
  filter(variable %in% c("l61_819258_mpnn2", "m971")) %>% 
  lm(data=.,value ~ variable + ET_Ratio) %>% summary()

valsdf %>%
  filter(Cellline == "Raji") %>%
  filter(variable %in% c("l61_819258_mpnn2", "m971")) %>% 
  lm(data=.,value ~ variable + ET_Ratio) %>% summary()
