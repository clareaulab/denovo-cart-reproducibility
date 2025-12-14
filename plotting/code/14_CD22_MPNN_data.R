library(BuenColors)
library(dplyr)
library(data.table)

raw_in <- fread("../data/CD22_MPNN_data.txt") %>% data.frame() %>% arrange(desc(ID))

# append nameing
raw_in$what <- c("m971", "UTD", "NB", rep("nonint", 12), rep("int", 11), "WT")

p1 <- ggplot(raw_in, aes(x = RPMI_8226, y = K562_CD22, color = what)) +
  geom_point(size = 1) + scale_color_manual(values = c("dodgerblue3", "purple", "grey", "firebrick","grey", "red" ))+
  geom_abline(slope=1, intercept=0, linetype = 2) +
  scale_y_continuous(limits = c(0, 80)) + scale_x_continuous(limits = c(0, 80)) +
  pretty_plot(fontsize = 8) + L_border() +
  geom_text(aes(label=ID)) + 
  theme(legend.position = "none")
p1
cowplot::ggsave2(p1, file = "../plots/scatter_CD22_RPMI.pdf", width = 1.3, height = 1.3)
  
melt_df <- raw_in[,c("ID","CAR_Only","RPMI_8226","K562_Parental","Raji", "K562_CD22")] %>%
  reshape2::melt(id.vars = c("ID")) 

order_cl <- c("CAR_Only","RPMI_8226","K562_Parental","Raji")
order_binder_controls <- c("UTD","No_Binder", "m971")

reorder_df <- raw_in[,c("ID","CAR_Only","RPMI_8226","K562_Parental","Raji","K562_CD22")] %>%
  mutate(class = case_when(
    RPMI_8226 > 20 ~ "RPMIhigh",
    TRUE ~ "zOK")) %>%
  mutate(metric = (Raji + K562_CD22)/2 -CAR_Only) %>%
  arrange(desc(class), desc(metric)) 

order_vec <- reorder_df %>% pull(ID)

order_binders <- c(order_binder_controls, order_vec[!(order_vec %in% order_binder_controls)])
new_name <- c("UTD", "No_Binder", "m971", paste0("D3.", toupper(substr(order_binders, 10,10))[4:27]))
new_name[order_binders == "819258_2_wiltype"] <- "D3"
new_name_final <- make.unique(new_name, sep = "")
new_name_final <- replace(new_name_final,new_name_final=="D3.N", "D3.N0")
new_name_final <- replace(new_name_final,new_name_final=="D3.I", "D3.I0")
data.frame(
  new_name_final, order_binders
)

p1 <- melt_df %>%
  mutate(variable = factor(as.character(variable), levels = rev(order_cl))) %>% 
  drop_na() %>%
  mutate(binder_name = factor(as.character(ID), levels = (order_binders))) %>% 
  arrange(binder_name) %>% mutate(new_name_final1=rep(new_name_final, each = 4)) %>% 
  mutate(new_name_final1 = factor(as.character(new_name_final1), levels = (new_name_final))) %>% 
  ggplot(.,aes(x = new_name_final1, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") + labs(x = "CD22 MPNN Draw", y = "co-culture line") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

p2 <- raw_in[,c("ID","X0nM", "X0.1nM", "X1nM")] %>%
  reshape2::melt(id.vars = c("ID")) %>% 
  mutate(variable = factor(as.character(variable), levels = rev(c("X0nM", "X0.1nM", "X1nM")))) %>% 
  mutate(binder_name = factor(as.character(ID), levels = (order_binders))) %>% 
  arrange(binder_name) %>% mutate(new_name_final1=rep(new_name_final, each = 3)) %>% 
  mutate(new_name_final1 = factor(as.character(new_name_final1), levels = (new_name_final))) %>%
  ggplot(.,aes(x = new_name_final1, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100))  +
  pretty_plot(fontsize = 8) + L_border() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2
pX <- cowplot::plot_grid(p2, p1, nrow = 2, align = "v", axis = "l", rel_heights = c(0.7, 1.3))
cowplot::ggsave2(pX, file = "../plots/MPNN_CAR_CD22.pdf", width = 4, height = 1.9)


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
pS

cowplot::ggsave2(pS, file = "../plots/concentrations3_CD22.pdf", width = 1.3, height = 1)

