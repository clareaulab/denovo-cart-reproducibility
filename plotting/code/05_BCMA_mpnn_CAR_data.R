library(BuenColors)
library(dplyr)
library(readxl)

raw_in <- read_xlsx("../data/bcma_mpnn_binders_fullattributes.xlsx") %>% data.frame()

melt_df <- raw_in[,c("binder_name","coculture.CAR.Only","coculture.K562","coculture.K562.BCMA.OE","coculture.Raji","coculture.RPMI.8226")] %>%
  reshape2::melt(id.vars = c("binder_name")) %>% mutate(variable = gsub("coculture.", "", variable))

order_cl <- c("CAR.Only", "K562", "RPMI.8226", "Raji", "K562.BCMA.OE")
order_binder_controls <- c("BCMA_UTD","BCMA_NB", "BCMA_Abecma")

reorder_df <- raw_in[,c("coculture.K562.BCMA.OE","coculture.CAR.Only", "coculture.Raji", "coculture.RPMI.8226", "binder_name")] %>%
  mutate(class = case_when(
    coculture.CAR.Only > 50 ~ "Tonic",
    TRUE ~ "zOK")) %>%
  mutate(metric = (coculture.K562.BCMA.OE + coculture.Raji + coculture.RPMI.8226)/3 -coculture.CAR.Only) %>%
  arrange(desc(class), desc(metric)) 

order_vec <- reorder_df %>% pull(binder_name)

order_binders <- c(order_binder_controls, order_vec[!(order_vec %in% order_binder_controls)])
new_name <- c("UTD", "NB", "Abecma", paste0("B5.", toupper(substr(order_binders, 10,10))[4:25]))
new_name[order_binders == "BCMA_l59_int_wildtype"] <- "B5"
new_name_final <- make.unique(new_name, sep = "")

data.frame(
  new_name_final, order_binders
)

p1 <- melt_df %>%
  mutate(variable = factor(as.character(variable), levels = rev(order_cl))) %>% 
  mutate(binder_name = factor(as.character(binder_name), levels = (order_binders))) %>% 
  arrange(binder_name) %>% mutate(new_name_final1=rep(new_name_final, each = 5)) %>% 
  mutate(new_name_final1 = factor(as.character(new_name_final1), levels = (new_name_final))) %>% 
  ggplot(.,aes(x = new_name_final1, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") + labs(x = "BCMA MPNN Draw", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p2 <- raw_in[,c("binder_name","staining_0_nM","staining_1_nM","staining_10_nM","staining_100_nM")] %>%
  reshape2::melt(id.vars = c("binder_name")) %>%  mutate(variable = gsub("staining_", "", variable)) %>%
  mutate(variable = factor(as.character(variable), levels = rev(c("0_nM","1_nM","10_nM","100_nM")))) %>% 
  mutate(binder_name = factor(as.character(binder_name), levels = (order_binders))) %>% 
  arrange(binder_name) %>% mutate(new_name_final1=rep(new_name_final, each = 4)) %>% 
  mutate(new_name_final1 = factor(as.character(new_name_final1), levels = (new_name_final))) %>%
  ggplot(.,aes(x = new_name_final1, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100))  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pX <- cowplot::plot_grid(p2, p1, nrow = 2, align = "v", axis = "l", rel_heights = c(0.7, 1.3))
cowplot::ggsave2(pX, file = "../plots/MPNN_CAR_BCMA.pdf", width = 4, height = 1.9)
