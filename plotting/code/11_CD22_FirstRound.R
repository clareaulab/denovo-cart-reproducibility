library(BuenColors)
library(dplyr)
library(data.table)

raw_in <- fread("../data/CD22-early-screen-data.tsv") %>% data.frame() %>% arrange(desc(Binder))

melt_df <- raw_in[,c("Binder", "CARs.Only","K562.Parental","Nalm6","RPMI.8226", "Raji")] %>%
  reshape2::melt(id.vars = c("Binder")) 

order_cl <- c("CARs.Only","K562.Parental","Nalm6","RPMI.8226", "Raji")
order_binder_controls <- c("UTD", "m971")

reorder_df <- raw_in[,c("Binder",order_cl)] %>%
  mutate(class = case_when(
    RPMI.8226 > 100 ~ "RPMIhigh",
    TRUE ~ "zOK")) %>%
  mutate(metric = Raji - CARs.Only) %>%
  arrange(desc(class), desc(metric)) 

order_vec <- reorder_df %>% pull(Binder)

order_binders <- c(order_binder_controls, order_vec[!(order_vec %in% order_binder_controls)])
new_name <- c("UTD", "m971", "C1", "C2", "C3", "C4")

data.frame(
  new_name, order_binders
)

p1 <- melt_df %>%
  mutate(variable = factor(as.character(variable), levels = rev(order_cl))) %>% 
  mutate(binder_name = factor(as.character(Binder), levels = (order_binders))) %>% 
  ggplot(.,aes(x = binder_name, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") + labs(x = "BCMA MPNN Draw", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

cowplot::ggsave2(p1, file = "../plots/CD22_initial_screen.pdf", width = 2, height = 1.2)

pS <- raw_in[,c("ID","X0nM", "X0.1nM", "X1nM", "X10nM")] %>%
  reshape2::melt(id.vars = c("ID")) %>%
  filter(ID %in% c("No_Binder", "m971","819258_2_wiltype", "CD22_l61_nonint_2261_CD22")) %>%
  mutate(variable = factor(as.character(variable), levels = rev(c("X0nM", "X0.1nM", "X1nM","X10nM")))) %>% 
  mutate(ID = factor(as.character(ID), levels = c("No_Binder", "m971","819258_2_wiltype", "CD22_l61_nonint_2261_CD22"))) %>% 
  ggplot(.,aes(x = ID, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100))  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y = "")

cowplot::ggsave2(pS, file = "../plots/concentrations3_CD22.pdf", width = 1.3, height = 1)

